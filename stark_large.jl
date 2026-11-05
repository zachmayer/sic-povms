# Export a companion-form SIC orbit by the Stark-unit method of stark.jl.
# Run: julia --project=. --threads=8 --gcthreads=1 --heap-size-hint=4G stark_large.jl DIMENSION units [BITS=2048]
# Then: julia --project=. --threads=8 --gcthreads=1 --heap-size-hint=4G stark_large.jl DIMENSION tensors [BITS=2048]
# Then run reconstruct_large.py as in README, after Julia exits and releases its heap.
#
# To keep large d tractable, factor the cocycle along the label orbit, evaluate
# the double sine by a full-precision series, and stream units and normal-basis traces.
# Build the ghost with 256 guard bits (the cocycle loses about 40 bits at d = 5779),
# then accept one FFT residual below 2^(-bits+40) at working precision. The series
# coefficients carry 128 more bits, so successive quadrature orders can agree below
# the guarded precision.
# Refine the conjugate basis in BigFloat: Float64 traces lose digits when the
# final tensor contraction cancels. reconstruct_large.py avoids a dense projector;
# refine_large.py supplies the published precision. Their headers give the methods.
# Rows share read-only labels; Julia 1.13 scoped precision and TLS-enabled MPFR/FLINT
# keep their arithmetic independent. Packed tensor files use bounded read buffers;
# threaded slices retain the same H*v. Commit an axis before deleting its input.
# A failed conjugation rereads the committed trace.
# Separate processes release the row workspace before tensor work.
#
# Canonical choices for dimensions without a representatives.toml pin:
# * For (d+1)(d-3)=D*f^2, take the least b>=3 with b^2-4=D*q^2, q|f, and
#   standard-Zauner symmetry. The form <1,-b,1> always exists at b=d-1;
#   its conductor can exceed stark.jl's unrestricted choice.
# * Labels are modulo d (odd d) or 2d (even d). Quotient by the full overlap
#   stabilizer, including a negative-norm unit when present, without conjugation.
# * Axes use increasing prime, decreasing exponent. Choose the first residue
#   pair (a,y) lexicographically for a*I+y*L with that prime-power order both
#   before and after quotienting by previous axes.
# * Try linear, shifted quadratic, then squared quadratic normal bases.
#   Zauner.radix enumerates shifts last index fastest; Fourier moduli exceed 2^-128.
# * As in stark.jl, maximize sigma(sum(u^2)). Choose each multiplication root
#   by its smallest (Re,Im) rounded to 10 decimals; fail on a tie.
# * The first trace-minus-one power K of L and first lexicographic v=(a,b)
#   with det([v,Kv])=1 give the standard-Zauner map. Rank reconstruction fits
#   rounded to 10 decimals, then column-major shift index; take the first
#   converging determinant class. Power iteration starts from e_0.
#
# Pin rationale: representatives.toml; test_large.py tests the unpinned rule at small d.

const cli_args = copy(ARGS)
# Empty ARGS so including stark.jl loads helpers without running its top-level loop.
empty!(ARGS)
# Nemo's threaded mode points GMP at plain malloc, process-wide: GMP/FLINT temporaries
# stop counting toward Julia's GC heap, which otherwise balloons after DFTs and
# defeats the heap hint. FLINT itself stays single-threaded.
ENV["NEMO_THREADED"]="1"
include("stark.jl")
using TOML, SHA
using FileWatching: Pidfile
import Nemo, QuadGK
using Nemo: AcbField, acb_struct
Nemo.__isthreaded[] || error("Nemo loaded before NEMO_THREADED was set")

const REPRESENTATIVES = TOML.parsefile(joinpath(@__DIR__, "representatives.toml"))
const LIVE_BITS = 128

function companion_form(d)
    d >= 4 || error("dimension must be at least four")
    D, f = coredisc((big(d) + 1) * (d - 3))
    for b in 3:d-1
        delta = big(b)^2 - 4
        delta % D == 0 || continue
        q = isqrt(div(delta, D))
        q^2 * D == delta && f % q == 0 || continue
        mod(d,9) == 3 && mod(div(f,q),3) == 0 && continue # skip F_a symmetry
        return AdmissibleTuple(d, QuadBin(1,-b,1))
    end
    error("no supported standard-Zauner companion form for d=$d")
end

function check_configuration(d,bits)
    bits >= 512 || error("construction requires at least 512 bits")
    names = ("stark.jl", "stark_large.jl", "representatives.toml", "Project.toml", "Manifest.toml")
    fingerprint = bytes2hex(sha256(join(read(joinpath(@__DIR__,name),String) for name in names)))
    value = "$d $bits $VERSION $(Sys.MACHINE) $fingerprint\n"
    path = "configuration.txt"
    if isfile(path)
        read(path,String) == value || error("checkpoint source/platform/settings changed; use a fresh LARGE_D_WORK")
    else
        open(path*".part","w") do io
            write_all(io,value)
            sync_file(io)
        end
        Base.Filesystem.rename(path*".part",path)
        sync_directory(".")
    end
end

"""Cyclic axes of (O/dd O)^× modulo the overlap stabilizer, in the header's order.
Represent a+y*beta by a*I+y*L; a negative-norm unit extends the stabilizer.
"""
function arithmetic_generators(F)
    d=Int(F.d); dd=2^iseven(d)*d; b=Int(-F.Q.b)
    if haskey(REPRESENTATIVES,string(d))
        pin=REPRESENTATIVES[string(d)]
        return [permutedims(reshape(Int.(v),2,2)) for v in pin["matrices"]],Int.(pin["orders"])
    end
    beta=(-F.Q.b+F.q*F.a)/2
    O=Order(F.K,[one(F.a),beta])
    R,_=quo(O,ideal(O,dd))
    UR,ur=unit_group(R)
    anti,A=is_antiunitary_with_generator(F)
    # A = a*I + y*L with L[2,2] = 0, so a = A[2,2] and y = A[2,1]; L acts as beta.
    generator = anti ? A[2,2]+A[2,1]*beta : beta
    G,pi=quo(UR,[preimage(ur,R(O(generator)))])
    # Integer pre-check for the search below: a+y*beta can have order n in G only
    # if its matrix to the n-th power lies in the stabilizer it quotients by.
    L=mod.(Int.(F.L),dd)
    stabilizer=[Matrix{Int}(I,2,2)]
    while (P=mod.(stabilizer[end]*(anti ? mod.(Int.(A),dd) : L),dd))!=I
        push!(stabilizer,P)
    end
    chosen=typeof(zero(G))[]
    mats,ords=Matrix{Int}[],Int[]
    remaining,projection=quo(G,chosen)
    while order(remaining)>1
        prime,power=first(sort(collect(factor(exponent(remaining)));by=first))
        n=Int(prime^power)
        function first_axis()
            for a in 0:dd-1,y in 0:dd-1
                gcd(a*a+b*a*y+y*y,dd)==1 || continue
                M=mod.(a*Matrix{Int}(I,2,2)+y*L,dd)
                foldl((P,_)->mod.(P*M,dd),1:n;init=Matrix{Int}(I,2,2)) in stabilizer || continue
                g=pi(preimage(ur,R(O(a+y*beta))))
                order(g)==n && order(projection(g))==n || continue
                return g,M
            end
            error("no direct cyclic complement of order $n")
        end
        g,M=first_axis()
        push!(chosen,g); push!(mats,M); push!(ords,n)
        remaining,projection=quo(G,chosen)
    end
    return mats,ords
end

"""Double-sine series for c=(β+1)/2, x=2(c-w), |x|≤1 (radius of convergence β+1):
log S₂(w;β,1) = -x*log(2)/2 + Σ_{j≥1} dⱼ*(x^(2j+1)-x), where
 dⱼ = ∫₀∞ (t/2)^(2j+1) / ((2j+1)!*2t*sinh(βt/2)*sinh(t/2)) dt
    = 2/(2j+1) * Σ_{m,n≥0} (2(c+mβ+n))^-(2j+1).
S₂ is Zauner.double_sine (Shintani's convention); AFK §8.2 writes its reciprocal.
S₂(w)S₂(β+1-w)=1 and S₂(w+1)=2sin(πw/β)S₂(w) at w=β/2 give S₂(β/2)=1/√2, so the
linear term needs no singular-integral subtraction.
Gauss–Legendre on dyadic intervals supplies the coefficients. Compare successive
orders in ℓ₁ norm; bound both the polynomial tail and the omitted integral.
These are numerical convergence checks, not interval certification. Large b
(b = 299 at 2048 bits) fails the convergence check; the dimensions here have b ≤ 7.
"""
function double_sine_coefficients(b)
    bits=precision(BigFloat)
    coeff=setprecision(BigFloat,bits+128) do
        β=(BigFloat(b)+sqrt(BigFloat(b)^2-4))/2
        c=(β+1)/2
        J=ceil(Int,(bits+96)/(2log2(2c)))
        T=nextpow(2,ceil(Int,2*(bits+96)*log(BigFloat(2))/β))
        # Bound Σ_{m,n≥0}(c+mβ+n)^(-3) by z3 using integrals.
        z3=c^-3+(1+1/β)/(2c^2)+1/(2β*c)
        polynomial_tail=z3/(2*(2J+3)*(2c)^(2J)*(1-(2c)^-2))
        integral_tail=4exp(-β*T/2)/(β*T*(1-exp(-β*T))*(1-exp(-BigFloat(T))))
        tolerance=BigFloat(2)^(-bits-16)
        polynomial_tail+integral_tail<tolerance || error("double-sine tail too large")
        function integrate(order)
            nodes,weights=QuadGK.gauss(BigFloat,order)
            dj=zeros(BigFloat,J)
            for (lo,hi) in [(0,1);[(2^k,2^(k+1)) for k in 0:trailing_zeros(T)-1]]
                mid=BigFloat(lo+hi)/2; half=BigFloat(hi-lo)/2
                for i in eachindex(nodes)
                    t=mid+half*nodes[i]
                    v=half*weights[i]*t^2*exp(-c*t)/(24*expm1(-β*t)*expm1(-t))
                    t2=t^2/4
                    for j in 1:J
                        dj[j]+=v
                        v*=t2/((2j+2)*(2j+3))
                    end
                end
            end
            [-log(BigFloat(2))/2-sum(dj);dj]
        end
        println("DOUBLE_SINE terms=$J cutoff=$T tail_log2=$(Float64(log2(polynomial_tail+integral_tail)))")
        low=integrate(512)
        for order in 768:256:2048
            high=integrate(order)
            delta=sum(abs,high-low)
            println("DOUBLE_SINE order=$order difference_log2=$(Float64(log2(delta)))")
            flush(stdout)
            delta<tolerance && return high
            low=high
        end
        error("double-sine coefficient quadrature did not converge")
    end
    BigFloat.(coeff)
end

"""stark.jl's sigma at full precision: shift to |x|≤1, evaluate the odd series, retain the
q-Pochhammer phase. The integer shift n differs from stark.jl's; the q-Pochhammer and phase
factors compensate, so the product is the same."""
function sigma_series(z,β,coeff)
    m=floor(Int,z/β)
    z-=m*β
    n=round(Int,(β+1)/2-z)-1
    x=β+1-2(z+n+1)
    abs(x)<=1+8eps(BigFloat) || error("double-sine argument outside series interval")
    qpoch(z,β,m)*qpoch(z/β,-1/β,-n)*
        e((6(z+n)^2+6(1-β)*(z+n)+β^2-3β+1)/(24β))*exp(x*evalpoly(x*x,coeff))
end

"""SF cocycle along the reduced label orbit of L=[b -1;1 0].

At each step Lp = pnext + d*(s1,0); the integer shift changes the
first label only, so the exponential q-Pochhammer correction is one.
The Jacobi cocycle telescopes to the modular cocycle for A=L^k.
AFK equations (1.21), (1.22), (1.26), (8.4), and (8.9).
"""
function orbit_shin(F, p1, p2, coeff)
    @assert F.Q.a == F.Q.c == 1
    @assert F.L == [Int(-F.Q.b) -1; 1 0]
    d, β, b = Int(F.d), F.x, Int(-F.Q.b)
    start = (p1,p2)
    val = one(Complex{BigFloat})
    for _ in 1:F.k
        val *= sigma_series((p2*β-p1)/d,β,coeff)
        p1,p2 = mod(b*p1-p2,d),p1
    end
    @assert (p1,p2) == start
    val
end

function ghost_nu(F,p1,p2,coeff)
    d,A = Int(F.d),F.A
    QA = QuadBin(A[2,1],A[2,2]-A[1,1],-A[1,2])
    ζ = -cispi(BigFloat(1)/d)
    c = -e(-BigFloat(rademacher(A))/24)/sqrt(BigFloat(d+1))
    s = iseven(d)*(p1+p2+p1*p2)
    (-1)^s * ζ^(-BigFloat(QA(p1,p2))/(d*(d-2))) * c * orbit_shin(F,p1,p2,coeff)
end

function cocycle_ghost(F,coeff)
    d = Int(F.d)
    χ = zeros(Complex{BigFloat},d,2)
    χ[1,1] = 1
    ζ = -cispi(BigFloat(1)/d)
    maximag = BigFloat(0)
    for (p1,p2) in Iterators.drop(Iterators.product(0:1,0:d-1),1)
        ν = ghost_nu(F,p1,p2,coeff)
        maximag = max(maximag,abs(imag(ν))/abs(ν))
        χ[p2+1,p1+1] = ζ^(p1*p2)*real(ν)
    end
    println("cocycle max relative imaginary part=$maximag")
    # 64 bits cover the phase ζ^(QA/(d(d-2))), whose argument grows like |A| (about 2^47.5 at
    # d = 57966). This checks phases only: the series enters as a real factor.
    maximag<BigFloat(2)^(-precision(BigFloat)+64) || error("cocycle imaginary part too large")
    χ = ifft(χ,1)
    ψ = sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]),1)
    ψ/ψ[1]
end

"Numerical DFT through the installed FLINT C API; return ball midpoints."
function flint_dft(x::Vector{Complex{BigFloat}};inverse=false)
    n=length(x)
    bits=precision(BigFloat)
    C=AcbField(bits+32)
    a=Nemo.acb_vec(C.(x))
    b=Nemo.acb_vec(n)
    if inverse
        ccall((:acb_dft_inverse,Nemo.libflint),Cvoid,
              (Ptr{acb_struct},Ptr{acb_struct},Int,Int),b,a,n,bits+32)
    else
        ccall((:acb_dft,Nemo.libflint),Cvoid,
              (Ptr{acb_struct},Ptr{acb_struct},Int,Int),b,a,n,bits+32)
    end
    vals=Nemo.array(C,b,n)
    Nemo.acb_vec_clear(a,n)
    Nemo.acb_vec_clear(b,n)
    [complex(BigFloat(real(z)),BigFloat(imag(z))) for z in vals]
end


"Quartic ghost residuals by four circular correlations instead of d sums per row."
function ghost_residual_fft(z,ns)
    d=1+length(z)÷2
    ψ=[complex(one(eltype(z)),zero(eltype(z)));complex.(z[1:d-1],z[d:end])]
    N=sum(real(ψ[mod(-k,d)+1]*conj(ψ[k+1])) for k in 0:d-1)
    rows=Dict{Int,Vector{Complex{eltype(z)}}}()
    for p in unique([radix(n,[2,d,d])[2] for n in ns])
        A=[ψ[mod(-k,d)+1]*conj(ψ[mod(k-p,d)+1]) for k in 0:d-1]
        B=[ψ[k+1]*conj(ψ[mod(p-k,d)+1]) for k in 0:d-1]
        rows[p]=flint_dft(flint_dft(A).*flint_dft(B);inverse=true)
    end
    map(ns) do n
        r,p,q=radix(n,[2,d,d])
        v=rows[p][mod(p+q,d)+1]-((p==0)+(q==0))*N^2/(d+1)
        r==0 ? real(v) : imag(v)
    end
end


"Fixed records: Int64 sign, Int64 MPFR exponent, then native UInt64 limbs."
struct PackedTensor{S,N}
    path::String
    shape::NTuple{N,Int}
    bits::Int
end
Base.size(T::PackedTensor)=T.shape
Base.size(T::PackedTensor,j)=T.shape[j]
Base.length(T::PackedTensor)=prod(T.shape)
Base.ndims(T::PackedTensor)=length(T.shape)
record_words(T::PackedTensor{S}) where S=(2+cld(T.bits,64))*(S===BigFloat ? 1 : 2)

function packed(path,shape,bits,::Type{S}=BigFloat;create=false) where S
    T=PackedTensor{S,length(shape)}(relpath(path),Tuple(shape),bits)
    bytes=8record_words(T)*length(T)
    if create
        open(path,"w") do io
            truncate(io,bytes)
        end
    end
    filesize(path)==bytes || error("wrong tensor file size: $path")
    T
end

# Copy limbs into an owned MPFR buffer: no pointer in the file or borrowed limb lifetime.
function unpack!(x::BigFloat,words,k)
    sign=reinterpret(Int64,words[k])
    sign in (-1,1) || error("unwritten record")
    x.sign=sign; x.exp=reinterpret(Int64,words[k+1])
    n=cld(precision(x),64)
    checkbounds(words,k+1+n)
    GC.@preserve x words unsafe_copyto!(Base.unsafe_convert(Ptr{UInt64},x.d),pointer(words,k+2),n)
    x
end
function pack!(words,k,x::BigFloat,bits)
    @assert precision(x)==bits
    words[k]=reinterpret(UInt64,Int64(x.sign)); words[k+1]=reinterpret(UInt64,Int64(x.exp))
    n=cld(precision(x),64)
    checkbounds(words,k+1+n)
    # Special values have no meaningful mantissa; canonicalize their unused bytes.
    if isfinite(x) && !iszero(x)
        GC.@preserve x words unsafe_copyto!(pointer(words,k+2),Base.unsafe_convert(Ptr{UInt64},x.d),n)
    else
        fill!(view(words,k+2:k+1+n),0)
    end
    words
end
function unpack!(z::Complex{BigFloat},words,k)
    unpack!(real(z),words,k); unpack!(imag(z),words,k+2+cld(precision(real(z)),64)); z
end
function pack!(words,k,z::Complex{BigFloat},bits)
    pack!(words,k,real(z),bits); pack!(words,k+2+cld(bits,64),imag(z),bits); words
end
buffer(::Type{BigFloat},bits)=BigFloat(0;precision=bits)
buffer(::Type{Complex{BigFloat}},bits)=complex(buffer(BigFloat,bits),buffer(BigFloat,bits))

# A bounded read cache also gives Base.sum its original Array pairwise reduction tree.
mutable struct PackedReader{S,N} <: AbstractArray{S,N}
    tensor::PackedTensor{S,N}
    io::IOStream
    words::Vector{UInt64}
    first::Int
end
function PackedReader(T::PackedTensor)
    count=max(1,div(1024^2,8record_words(T)))
    PackedReader(T,open(T.path,"r"),Vector{UInt64}(undef,count*record_words(T)),0)
end
Base.size(R::PackedReader)=size(R.tensor)
Base.IndexStyle(::Type{<:PackedReader})=IndexLinear()
Base.close(R::PackedReader)=close(R.io)
function read_record!(x,R::PackedReader,i)
    @assert precision(real(x))==R.tensor.bits
    @assert x isa BigFloat || precision(imag(x))==R.tensor.bits
    w=record_words(R.tensor); count=length(R.words)÷w
    if !(R.first<=i<R.first+count) || R.first==0
        R.first=div(i-1,count)*count+1
        seek(R.io,8w*(R.first-1))
        read!(R.io,view(R.words,1:w*min(count,length(R)-R.first+1)))
    end
    unpack!(x,R.words,1+w*(i-R.first))
end
# Pairwise reduction holds two values at once, so each index needs a fresh buffer.
Base.getindex(R::PackedReader{S},i::Int) where S=read_record!(buffer(S,R.tensor.bits),R,i)
function reading(f,T::PackedTensor)
    R=PackedReader(T)
    try
        f(R)
    finally
        close(R)
    end
end
function write_all(io,x)
    write(io,x)==sizeof(x) || error("short write")
    flush(io)
end
function sync_fd(descriptor)
    # F_FULLFSYNC reaches the macOS drive; fsync elsewhere.
    status=Sys.isapple() ? @ccall(fcntl(descriptor::Cint,51::Cint;0::Cint)::Cint) :
                           @ccall(fsync(descriptor::Cint)::Cint)
    status==0 || error("sync failed")
end
function sync_file(io)
    flush(io)
    sync_fd(fd(io))
end

"Selected real ghost overlaps; disjoint row writes and durable completion markers."
function ghost_units_stream(ψ,P,path)
    d=length(ψ)
    φ=circshift(reverse(ψ),1)
    normfactor=(d+1)/real(φ'ψ)
    ζ=[(-cispi(big(1.0)/d))^(-j) for j in 0:2d-1]
    rows=[Int[] for _ in 1:d]
    for (i,(a,b)) in enumerate(P)
        push!(rows[mod(a,d)+1],i)
    end
    state=path*".rows"
    done=isfile(state) ? deserialize(state) : fill(false,cld(d,32))
    T=packed(path,size(P),precision(BigFloat);create=!isfile(state))
    if !isfile(state)
        open(sync_file,path,"r+")
        save_atomic(state,done)
    end
    commit=ReentrantLock()
    Threads.@threads :dynamic for chunk in eachindex(done)
        done[chunk] && continue
        firstrow=32(chunk-1); lastrow=min(d-1,firstrow+31)
        open(path,"r+") do io
            words=zeros(UInt64,record_words(T))
            for a in firstrow:lastrow
                isempty(rows[a+1]) && continue
                vals=[conj(φ[k+1])*ψ[mod(k-a,d)+1] for k in 0:d-1]
                row=d*flint_dft(vals;inverse=true)
                for i in rows[a+1]
                    p,q=P[i]
                    v=normfactor*real(ζ[mod(p*q,2d)+1]*row[mod(q,d)+1])
                    pack!(words,1,v,T.bits)
                    seek(io,8length(words)*(i-1)); write_all(io,words)
                end
            end
            lock(commit) do
                sync_file(io)
                done[chunk]=true
                save_atomic(state,done)
            end
        end
        println("ghost overlap rows=$firstrow:$lastrow/$d"); flush(stdout)
    end
    T
end

"Marginal normal-basis sums without materializing a shifted Galois tensor."
function basis_marginal(T,j,recipe)
    reading(T) do u
    reading(T) do shifted_u
    a=buffer(BigFloat,T.bits); b=buffer(BigFloat,T.bits)
    shape=size(u)
    result=zeros(BigFloat,shape[j])
    linear=recipe.kind===:linear
    squared=recipe.kind===:shifted_quad_sq
    for ci in CartesianIndices(u)
        read_record!(a,u,LinearIndices(u)[ci])
        v = if linear
            a
        else
            shifted=CartesianIndex(mod.(Tuple(ci).-recipe.shift.-1,shape).+1)
            read_record!(b,shifted_u,LinearIndices(u)[shifted])
            a+a*b
        end
        squared && (v=v^2)
        result[ci[j]]+=v
    end
    result
    end
    end
end

function streamed_normal_recipes(u,d)
    map(1:ndims(u)) do j
        isbasis(recipe)=minimum(abs,flint_dft(complex.(basis_marginal(u,j,recipe));inverse=true))>ldexp(big(1.0),-128)
        linear=Zauner.NormalBasisRecipe(:linear,nothing)
        isbasis(linear) && return linear
        for kind in (:shifted_quad,:shifted_quad_sq), m in 0:length(u)-1
            recipe=Zauner.NormalBasisRecipe(kind,Tuple(Zauner.radix(m,size(u))))
            isbasis(recipe) && return recipe
        end
        error("d=$d no normal basis on axis$j")
    end
end

function normal_bases(u,recipes)
    [basis_marginal(u,j,recipes[j]) for j in 1:ndims(u)]
end

"Transform bounded blocks of independent slices, preserving exactly H*v."
function tensor_axis!(out::PackedTensor,T::PackedTensor{S},H,j) where S
    out.bits==precision(BigFloat) || error("axis output precision differs from working precision")
    size(out)==size(T) || error("axis tensor shapes differ")
    n=size(T,j); stride=prod(size(T)[1:j-1])
    width=min(stride,1024); blocks=cld(stride,width)
    batch=max(1,1024÷stride); outer_count=length(T)÷(n*stride)
    jobs=cld(outer_count,batch)*blocks
    Threads.@threads :dynamic for worker in 1:Threads.nthreads()
        open(T.path,"r") do src
        open(out.path,"r+") do dst
            a=Vector{UInt64}(undef,record_words(T)*n*width*batch)
            b=Vector{UInt64}(undef,record_words(out)*n*width*batch)
            v=[buffer(S,T.bits) for _ in 1:n]
            for job in worker:Threads.nthreads():jobs
                outer,inner=divrem(job-1,blocks)
                first=inner*width; count=min(width,stride-first)
                many=min(batch,outer_count-outer*batch)
                offset=outer*batch*n*stride+first
                # Early axes: read several complete contiguous slabs at once.
                # Later axes: gather one contiguous block along each axis position.
                if blocks==1
                    seek(src,8record_words(T)*offset)
                    read!(src,view(a,1:record_words(T)*n*count*many))
                else
                    for k in 0:n-1
                        seek(src,8record_words(T)*(offset+k*stride))
                        read!(src,view(a,1+k*count*record_words(T):(k+1)*count*record_words(T)))
                    end
                end
                for outer_in in 0:many-1, i in 0:count-1
                    for k in 0:n-1
                        index=outer_in*n*count+k*count+i
                        unpack!(v[k+1],a,1+index*record_words(T))
                    end
                    w=H*v
                    for k in 0:n-1
                        index=outer_in*n*count+k*count+i
                        pack!(b,1+index*record_words(out),w[k+1],out.bits)
                    end
                end
                if blocks==1
                    seek(dst,8record_words(out)*offset)
                    write_all(dst,view(b,1:record_words(out)*n*count*many))
                else
                    for k in 0:n-1
                        seek(dst,8record_words(out)*(offset+k*stride))
                        write_all(dst,view(b,1+k*count*record_words(out):(k+1)*count*record_words(out)))
                    end
                end
            end
        end
        end
    end
    out
end

function trace_product(T,x,prefix,rows_path)
    state=prefix*".jls"
    last=0
    if isfile(state)
        last,T=deserialize(state)
        println("RESUME trace at axis=$(last+1)")
    end
    # A crash after committing state can leave an obsolete predecessor on disk.
    for j in 0:last-1
        old=j==0 ? rows_path : "$(prefix)_axis$j.bin"
        isfile(old) && rm(old)
    end
    for j in last+1:ndims(T)
        out=packed("$(prefix)_axis$j.bin",size(T),T.bits;create=true)
        t=@elapsed tensor_axis!(out,T,hankel(x[j]),j)
        open(sync_file,out.path,"r+")
        save_atomic(state,(j,out))
        rm(T.path)
        T=out
        println("TRACE axis=$j seconds=$t maxrss=$(Sys.maxrss())"); flush(stdout)
    end
    T
end

"Keep guard bits through the ill-conditioned live-basis tensor contraction."
function conjugate_units(σ,trace,forms_x)
    # The committed trace is the retry source; discard interrupted live files.
    for path in ["recognized.bin";["live_axis$j.bin" for j in 1:ndims(trace)]]
        isfile(path) && rm(path)
    end
    y=Vector{Vector{Complex{BigFloat}}}()
    for pair in forms_x
        G,C=(X->σ.(X)).(pair)
        # A failed relation or basis check asks the caller for more precision.
        issuccess(lu(G;check=false)) || return nothing
        A=C/G
        A64=Float64.(A)
        all(isfinite,A64) || return nothing
        λ,V=eigen(A64)
        keys=[(round(real(z);digits=10),round(imag(z);digits=10)) for z in λ]
        allunique(keys) || error("multiplication roots tie at 10 decimals")
        index=argmin(keys)
        # z is an eigenvector of A scaled so that its first entry is its eigenvalue.
        # Solve A*z-z[1]*z=0 with Jacobian A-z[1]*I-z*e₁ᵀ, starting at the selected Float64 root.
        z=Complex{BigFloat}.(V[:,index]*λ[index]/V[1,index])
        for _ in 1:5
            residual=A*z-z[1]*z
            jacobian=complex.(A)-z[1]*I
            jacobian[:,1]-=z
            z-=jacobian\residual
        end
        Gz,Cz=forms(z)
        tolerance=ldexp(big(1.0),-(LIVE_BITS-32))
        max(norm(Gz-G)/norm(G),norm(Cz-C)/norm(C))<tolerance || return nothing
        push!(y,G\z)
    end
    # Recognize at working precision; round only its result to the live precision.
    T=packed("recognized.bin",size(trace),LIVE_BITS;create=true)
    reading(trace) do input
        open(T.path,"r+") do io
            v=buffer(BigFloat,trace.bits); words=zeros(UInt64,record_words(T))
            for i in 1:length(trace)
                read_record!(v,input,i)
                pack!(words,1,BigFloat(σ(v);precision=LIVE_BITS),T.bits)
                write_all(io,words)
            end
        end
    end
    setprecision(BigFloat,LIVE_BITS) do
        for j in 1:ndims(T)
            out=packed("live_axis$j.bin",size(T),LIVE_BITS,Complex{BigFloat};create=true)
            tensor_axis!(out,T,hankel(y[j]),j)
            rm(T.path)
            T=out
        end
    end
    T
end

function generator_orbit(F,mats,ords)
    d=Int(F.d); dd=2^iseven(d)*d
    all(M->all(x->0<=x<dd,M),mats) || throw(ArgumentError("generator matrices must contain canonical residues"))
    powers=[[Matrix{Int}(I,2,2)] for _ in mats]
    for i in eachindex(mats)
        for k in 1:ords[i]-1
            push!(powers[i],mod.(powers[i][end]*mats[i],dd))
        end
    end
    P=Array{NTuple{2,Int}}(undef,Tuple(ords))
    for index in CartesianIndices(P)
        a,b=1,0
        for j in eachindex(mats)
            M=powers[j][index[j]]
            next_a=mod(M[1,1]*a+M[1,2]*b,dd)
            next_b=mod(M[2,1]*a+M[2,2]*b,dd)
            a,b=next_a,next_b
        end
        P[index]=(a,b)
    end
    P
end

function check_generators(F,mats,ords,P)
    d=Int(F.d); dd=2^iseven(d)*d
    stabs=[Int.(mod.(S,dd)) for S in stabilizer_elements(F)]
    for (M,n) in zip(mats,ords)
        gcd(M[1,1]*M[2,2]-M[1,2]*M[2,1],dd)==1 || error("nonunit quotient generator")
        mod.(M*Int.(F.L)-Int.(F.L)*M,dd)==zeros(Int,2,2) || error("noncommuting quotient generator")
        power=Matrix{Int}(I,2,2)
        for _ in 1:n
            power=mod.(power*M,dd)
        end
        @assert power in stabs
    end
    seen=falses(dd,dd)
    for p in P,s in stabs
        a=mod(s[1,1]*p[1]+s[1,2]*p[2],dd)+1
        b=mod(s[2,1]*p[1]+s[2,2]*p[2],dd)+1
        @assert !seen[a,b]
        seen[a,b]=true
    end
    println("group axes=$ords size=$(length(P)) stabilizers=$(length(stabs)) distinct labels=$(count(seen))")
    stabs
end

function standard_map(F)
    d=Int(F.d); dd=2^iseven(d)*d
    K=Matrix{Int}(I,2,2)
    for _ in 1:2F.k
        K=mod.(K*Int.(F.L),dd)
        mod(tr(K)+1,dd)==0 || continue
        for a in 0:dd-1,b in 0:dd-1
            c,e=mod(K[1,1]*a+K[1,2]*b,dd),mod(K[2,1]*a+K[2,2]*b,dd)
            mod(a*e-b*c,dd)==1 && return mod.([e -c;-b a],dd)
        end
    end
    error("no standard-Zauner conjugator")
end


"Write, sync to the drive, then rename: later stages delete the predecessor."
function save_atomic(path,x)
    data=IOBuffer()
    serialize(data,x)
    open(path*".part","w") do io
        write_all(io,take!(data))
        sync_file(io)
    end
    Base.Filesystem.rename(path*".part",path)
    sync_directory(dirname(abspath(path)))
    x
end

"Persist directory entries before deleting a committed file's predecessor."
function sync_directory(path)
    directory=ccall(:open,Cint,(Cstring,Cint),path,0)
    directory>=0 || error("cannot open checkpoint directory")
    try
        sync_fd(directory)
    finally
        ccall(:close,Cint,(Cint,),directory)
    end
end

function stage(f,name)
    path=name*".jls"
    if isfile(path)
        println("RESUME $name"); flush(stdout)
        return deserialize(path)
    end
    GC.gc()
    start=time()
    x=f()
    GC.gc() # discard stage workspace before allocating the checkpoint's identity table
    save_atomic(path,x)
    println("STAGE name=$name seconds=$(time()-start) maxrss=$(Sys.maxrss())")
    flush(stdout)
    x
end

function run_construct(d,bits,phase)
    phase in ("units","tensors") || error("phase must be units or tensors")
    root=pwd()
    work=abspath(get(ENV,"LARGE_D_WORK",joinpath(root,"large_d_work","d$d")))
    mkpath(work)
    cd(work) do
    Pidfile.mkpidlock("lock";wait=false) do
    BLAS.set_num_threads(parse(Int,get(ENV,"SIC_BLAS_THREADS","1")))
    randseed!(0); seed!(0)
    F=setprecision(BigFloat,bits+256) do
        F=companion_form(d)
        β=(BigFloat(-F.Q.b)+sqrt(BigFloat(F.Q.b)^2-4))/2
        @assert precision(F.x)==bits+256 && F.x==β
        F
    end
    check_configuration(d,bits)
    if isfile("export.jls")
        all(name->isfile(joinpath("export",name)),("units.bin","labels.bin","metadata.json")) || error("committed export is missing")
        for path in deserialize("export.jls")
            isfile(path) && rm(path)
        end
        println("RESUME export ",abspath("export")); flush(stdout)
        return
    end
    phase=="tensors" && !any(isfile,("units_$bits.jls","trace_$bits.jls")) &&
        error("run the units phase first")
    @assert !has_fa_symmetry(F)
    println("START d=$d phase=$phase q=$(F.q) f=$(F.f) bits=$bits threads=$(Threads.nthreads()) heap=$(Base.JLOptions().heap_size_hint) blas=$(BLAS.get_num_threads())")
    flush(stdout)
    mats,ords,P,S=stage("group") do
        mats,ords=arithmetic_generators(F)
        P=generator_orbit(F,mats,ords)
        S=check_generators(F,mats,ords,P)
        mats,ords,P,S
    end
    ψ=stage("ghost_$bits") do
        guarded=setprecision(BigFloat,bits+256) do
            coeff=double_sine_coefficients(Int(-F.Q.b))
            cocycle_ghost(F,coeff)
        end
        setprecision(BigFloat,bits) do
            # Round explicitly: Complex{BigFloat}.() would keep the guarded precision.
            ψ=[complex(BigFloat(real(z);precision=bits),BigFloat(imag(z);precision=bits)) for z in guarded]
            residual=maximum(abs,ghost_residual_fft(re_im_proj(ψ),equations(d)))
            println("GHOST residual_log2=$(Float64(log2(residual)))")
            residual<BigFloat(2)^(-bits+40) || error("direct ghost residual too large")
            ψ
        end
    end
    setprecision(BigFloat,bits) do
        trace_path="trace_$bits.jls"
        units_path="units_$bits.jls"
        rows_path="units_$bits.bin"
        u=if isfile(trace_path) || (phase=="units" && isfile(units_path))
            nothing
        else
            stage("units_$bits") do
                ghost_units_stream(ψ,P,rows_path)
            end
        end
        phase=="units" && return
        recipes=stage("recipes_$bits") do
            streamed_normal_recipes(u,d)
        end
        x=stage("bases_$bits") do
            normal_bases(u,recipes)
        end
        s,unit_max=stage("unit_stats_$bits") do
            reading(u) do values
                sum(abs2,values),maximum(abs,values)
            end
        end
        start=time()
        trace=trace_product(u,x,"trace_$bits",rows_path)
        u=nothing
        GC.gc()
        println("STAGE name=trace_$bits seconds=$(time()-start) maxrss=$(Sys.maxrss())")
        flush(stdout)
        forms_x=stage("forms_$bits") do
            forms.(x)
        end
        top=max(reading(v->maximum(abs,v),trace),maximum(maximum(abs,X) for X in Iterators.flatten(forms_x)))
        relation_bits=2Int(F.h)*(ceil(Int,log2(top))+40)
        println("PRECISION relation_bits=$relation_bits working_bits=$bits unit_max_log2=$(Float64(log2(unit_max)))")
        flush(stdout)
        1.5relation_bits+100<bits || error("Need more ghost precision: rerun with at least $(ceil(Int,1.5relation_bits+200)) bits")
        # Hecke's embedding caches contain process-local C state. Rebuild this
        # small object; Julia serialization is safe for the numerical tensors.
        H=classfield(F)
        uu=nothing
        # Relative squared-modulus guard, with room for the final Float64 rounding.
        squared_modulus_tolerance=64eps(Float64)*(d+1)
        for b in unique(min.(round.(Int,relation_bits.*1.5.^(0:3)),bits-100))
            uu=stage("liveunits_$(bits)_$b") do
                live=conjugate_units(conjugator(H,F,b,s),trace,forms_x)
                live===nothing && return nothing
                err=reading(v->maximum(z->abs(abs2(ComplexF64(z))-(d+1)),v),live)
                println("CONJUGATION bits=$b squared_modulus_error=$err"); flush(stdout)
                if err<squared_modulus_tolerance
                    open(sync_file,live.path,"r+")
                    return live
                end
                rm(live.path)
                nothing
            end
            uu!==nothing && break
        end
        uu!==nothing || error("Conjugation squared-modulus check failed")
        export_dir="export"
        mkpath(export_dir)
        reading(uu) do values
            open(joinpath(export_dir,"units.bin.part"),"w") do io
                for z in values
                    write_all(io,ComplexF64(z))
                end
                sync_file(io)
            end
        end
        Base.Filesystem.rename(joinpath(export_dir,"units.bin.part"),joinpath(export_dir,"units.bin"))
        open(io->write_all(io,vec(P)),joinpath(export_dir,"labels.bin"),"w")
        M=standard_map(F)
        matrix_list(A)="["*join([A[1,1],A[1,2],A[2,1],A[2,2]],",")*"]"
        shift=haskey(REPRESENTATIVES,string(d)) ? REPRESENTATIVES[string(d)]["shift"] : -1
        metadata="{\"d\":$d,\"q_b\":$(Int(F.Q.b)),\"preferred_shift\":$shift,\"axes\":[$(join(ords,","))],\"standard_map\":$(matrix_list(M)),\"stabilizers\":[$(join(matrix_list.(S),","))]}"
        open(io->write_all(io,metadata*"\n"),joinpath(export_dir,"metadata.json"),"w")
        for name in ("labels.bin","metadata.json")
            open(sync_file,joinpath(export_dir,name),"r+")
        end
        sync_directory(export_dir)
        save_atomic("export.jls",(trace.path,uu.path))
        rm(trace.path); rm(uu.path)
        println("EXPORT ",abspath(export_dir))
        flush(stdout)
    end
    end
    end
end

# Run only as a script, so test_conjugation.jl can include the helpers.
if abspath(PROGRAM_FILE) == @__FILE__
    length(cli_args) in (2,3) || error("usage: stark_large.jl DIMENSION units|tensors [BITS=2048]")
    run_construct(parse(Int,cli_args[1]),length(cli_args)==3 ? parse(Int,cli_args[3]) : 2048,cli_args[2])
    println("PROCESS maxrss=$(Sys.maxrss())")
end
