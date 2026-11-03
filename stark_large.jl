# Build a companion-form SIC seed by the Stark-unit method of stark.jl.
# Run: julia --project=. --threads=1 --gcthreads=1 stark_large.jl DIMENSION [BITS=2048]
#
# To keep large d tractable, factor the cocycle along the label orbit, use FLINT
# FFT residuals with a frozen Float64 QR, and stream units and normal-basis traces.
# Refine the conjugate basis in BigFloat: Float64 traces lose digits when the
# final tensor contraction cancels. reconstruct_large.py avoids a dense projector;
# refine_large.py supplies the published precision. Their headers give the methods.
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
include("stark.jl")
using TOML, SHA
import Nemo
using Nemo: AcbField, acb_struct

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
        write(path,value)
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
    chosen=typeof(zero(G))[]
    mats,ords=Matrix{Int}[],Int[]
    remaining,projection=quo(G,chosen)
    while order(remaining)>1
        prime,power=first(sort(collect(factor(exponent(remaining)));by=first))
        n=Int(prime^power)
        function first_axis()
            for a in 0:dd-1,y in 0:dd-1
                gcd(a*a+b*a*y+y*y,dd)==1 || continue
                g=pi(preimage(ur,R(O(a+y*beta))))
                order(g)==n && order(projection(g))==n || continue
                return g,mod.(a*Matrix{Int}(I,2,2)+y*Int.(F.L),dd)
            end
            error("no direct cyclic complement of order $n")
        end
        g,M=first_axis()
        push!(chosen,g); push!(mats,M); push!(ords,n)
        remaining,projection=quo(G,chosen)
    end
    return mats,ords
end

"""SF cocycle along the reduced label orbit of L=[b -1;1 0].

At each step Lp = pnext + d*(s1,0); the integer shift changes the
first label only, so the exponential q-Pochhammer correction is one.
The Jacobi cocycle telescopes to the modular cocycle for A=L^k.
AFK equations (1.21), (1.22), (1.26), (8.4), and (8.9).
"""
function orbit_shin(F, p1, p2)
    @assert F.Q.a == F.Q.c == 1
    @assert F.L == [Int(-F.Q.b) -1; 1 0]
    d, β, b = Int(F.d), F.x, Int(-F.Q.b)
    start = (p1,p2)
    val = one(Complex{BigFloat})
    for _ in 1:F.k
        val *= sigma((p2*β-p1)/d,β)
        p1,p2 = mod(b*p1-p2,d),p1
    end
    @assert (p1,p2) == start
    val
end

function ghost_nu(F,p1,p2)
    d,A = Int(F.d),F.A
    QA = QuadBin(A[2,1],A[2,2]-A[1,1],-A[1,2])
    ζ = -cispi(BigFloat(1)/d)
    c = -e(-BigFloat(rademacher(A))/24)/sqrt(BigFloat(d+1))
    s = iseven(d)*(p1+p2+p1*p2)
    (-1)^s * ζ^(-BigFloat(QA(p1,p2))/(d*(d-2))) * c * orbit_shin(F,p1,p2)
end

function cocycle_ghost(F)
    d = Int(F.d)
    χ = zeros(Complex{BigFloat},d,2)
    χ[1,1] = 1
    ζ = -cispi(BigFloat(1)/d)
    maximag = BigFloat(0)
    for (p1,p2) in Iterators.drop(Iterators.product(0:1,0:d-1),1)
        ν = ghost_nu(F,p1,p2)
        maximag = max(maximag,abs(imag(ν))/abs(ν))
        χ[p2+1,p1+1] = ζ^(p1*p2)*real(ν)
    end
    println("cocycle max relative imaginary part=$maximag")
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


"""Full-space high-precision ghost residuals with a frozen Float64 QR solve.
Scale residuals before conversion to avoid underflow. No inexact symmetry basis
constrains the final vector. Abort if the residual fails to halve at a step;
the Jacobian is not refreshed.
"""
function refine_ghost(ψ,bits)
    setprecision(BigFloat,bits) do
        z = BigFloat.(re_im_proj(ψ))
        ns = equations(length(ψ))
        println("building Float64 Jacobian $(length(ns)) x $(length(z))")
        flush(stdout)
        t = @elapsed J=qr!(ghost_jacobian(Float64.(z),ns))
        diagonal=[abs(J.factors[i,i]) for i in 1:length(z)]
        println("QR seconds=$t diagonal_ratio=$(minimum(diagonal)/maximum(diagonal))")
        flush(stdout)
        prev = big(Inf)
        maxsteps=max(100,bits÷20)
        for iter in 0:maxsteps
            t = @elapsed r=ghost_residual_fft(z,ns)
            res = maximum(abs,r)
            println("bump iter=$iter bits=$bits residual=$(Float64(log2(res))) log2 seconds=$t")
            flush(stdout)
            path="data/ghost_$(bits)_progress.jls"
            save_atomic(path,re_im_proj(z))
            if res < big(2.0)^(-bits+40)
                return re_im_proj(z)
            end
            res < prev/2 || error("Float64 inverse stalled at iteration $iter")
            scale = big(2.0)^floor(Int,log2(res))
            delta = J\Float64.(r/scale)
            z -= scale*BigFloat.(delta)
            prev=res
        end
        error("ghost refinement exhausted $maxsteps steps")
    end
end

"Selected real ghost overlaps, one FLINT transform per row; discard unused entries."
function ghost_units_stream(ψ,P,checkpoint_dir)
    d=length(ψ)
    φ=circshift(reverse(ψ),1)
    normfactor=(d+1)/real(φ'ψ)
    ζ=[(-cispi(big(1.0)/d))^(-j) for j in 0:2d-1]
    rows=[Int[] for _ in 1:d]
    for (i,(a,b)) in enumerate(P)
        push!(rows[mod(a,d)+1],i)
    end
    mkpath(checkpoint_dir)
    u=Array{BigFloat}(undef,size(P))
    for firstrow in 0:32:d-1
        lastrow=min(d-1,firstrow+31)
        path=joinpath(checkpoint_dir,"rows_$(firstrow)_$(lastrow).jls")
        if isfile(path)
            indices,values=deserialize(path)
        else
            indices=Int[]
            values=BigFloat[]
            for a in firstrow:lastrow
                isempty(rows[a+1]) && continue
                vals=[conj(φ[k+1])*ψ[mod(k-a,d)+1] for k in 0:d-1]
                row=d*flint_dft(vals;inverse=true)
                for i in rows[a+1]
                    p,q=P[i]
                    push!(values,normfactor*real(ζ[mod(p*q,2d)+1]*row[mod(q,d)+1]))
                end
                append!(indices,rows[a+1])
            end
            save_atomic(path,(indices,values))
            println("ghost overlap rows=$firstrow:$lastrow/$d")
            flush(stdout)
        end
        u[indices]=values
    end
    return u
end


"Marginal normal-basis sums without materializing a shifted Galois tensor."
function basis_marginal(u,j,recipe)
    shape=size(u)
    result=zeros(BigFloat,shape[j])
    linear=recipe.kind===:linear
    squared=recipe.kind===:shifted_quad_sq
    for ci in CartesianIndices(u)
        v = if linear
            u[ci]
        else
            shifted=CartesianIndex(mod.(Tuple(ci).-recipe.shift.-1,shape).+1)
            u[ci]+u[ci]*u[shifted]
        end
        squared && (v=v^2)
        result[ci[j]]+=v
    end
    result
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

function trace_product(u,x,prefix)
    T=u
    start=1
    path="$(prefix)_rolling.jls"
    if isfile(path)
        last,T=deserialize(path)
        start=last+1
        println("RESUME trace at axis=$start")
    end
    for j in start:ndims(u)
        t=@elapsed begin
            H=hankel(x[j])
            T=mapslices(v->H*v,T;dims=j)
        end
        save_atomic(path,(j,T))
        GC.gc()
        println("TRACE axis=$j seconds=$t maxrss=$(Sys.maxrss())")
        flush(stdout)
    end
    T
end

"Keep guard bits through the ill-conditioned live-basis tensor contraction."
function conjugate_units(σ,trace,forms_x)
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
    # Recognition must retain the caller's working precision. Only its result
    # is rounded; lowering precision around σ itself corrupts integer rounding.
    transformed=map(v->BigFloat(σ(v);precision=LIVE_BITS),trace)
    live=setprecision(BigFloat,LIVE_BITS) do
        hankel_prod(transformed,y)
    end
    ComplexF64.(live)
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


function save_atomic(path,x)
    serialize(path*".part",x)
    mv(path*".part",path;force=true)
    x
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
    save_atomic(path,x)
    println("STAGE name=$name seconds=$(time()-start) maxrss=$(Sys.maxrss())")
    flush(stdout)
    x
end

function run_construct(d,bits)
    root=pwd()
    work=abspath(get(ENV,"LARGE_D_WORK",joinpath(root,"large_d_work","d$d")))
    mkpath(joinpath(work,"data"))
    cd(work)
    BLAS.set_num_threads(parse(Int,get(ENV,"SIC_BLAS_THREADS","1")))
    randseed!(0); seed!(0)
    F=companion_form(d)
    check_configuration(d,bits)
    @assert !has_fa_symmetry(F)
    println("START d=$d q=$(F.q) f=$(F.f) bits=$bits blas=$(BLAS.get_num_threads())")
    flush(stdout)
    mats,ords,P,S=stage("group") do
        mats,ords=arithmetic_generators(F)
        P=generator_orbit(F,mats,ords)
        S=check_generators(F,mats,ords,P)
        mats,ords,P,S
    end
    initial=stage("initial") do
        setprecision(BigFloat,128) do
            cocycle_ghost(F)
        end
    end
    ψ=stage("ghost_$bits") do
        progress="data/ghost_$(bits)_progress.jls"
        if isfile(progress)
            start=deserialize(progress)
        else
            start=initial
        end
        refine_ghost(start,bits)
    end
    setprecision(BigFloat,bits) do
        u=stage("units_$bits") do
            ghost_units_stream(ψ,P,"units_$(bits)_rows")
        end
        recipes=stage("recipes_$bits") do
            streamed_normal_recipes(u,d)
        end
        x=stage("bases_$bits") do
            normal_bases(u,recipes)
        end
        trace=stage("trace_$bits") do
            trace_product(u,x,"trace_$bits")
        end
        forms_x=stage("forms_$bits") do
            forms.(x)
        end
        top=max(maximum(abs,trace),maximum(maximum(abs,X) for X in Iterators.flatten(forms_x)))
        relation_bits=2Int(F.h)*(ceil(Int,log2(top))+40)
        println("PRECISION relation_bits=$relation_bits working_bits=$bits unit_max_log2=$(Float64(log2(maximum(abs,u))))")
        flush(stdout)
        1.5relation_bits+100<bits || error("Need more ghost precision: rerun with at least $(ceil(Int,1.5relation_bits+200)) bits")
        s=sum(abs2,u)
        u=nothing
        GC.gc()
        # Hecke's embedding caches contain process-local C state. Rebuild this
        # small object; Julia serialization is safe for the numerical tensors.
        H=classfield(F)
        uu,squared_modulus_error=nothing,Inf
        # Relative squared-modulus guard, with room for the final Float64 rounding.
        squared_modulus_tolerance=64eps(Float64)*(d+1)
        for b in unique(min.(round.(Int,relation_bits.*1.5.^(0:3)),bits-100))
            uu=stage("liveunits_$(bits)_$b") do
                conjugate_units(conjugator(H,F,b,s),trace,forms_x)
            end
            if uu===nothing
                println("CONJUGATION bits=$b insufficient relation precision; advancing precision")
                flush(stdout)
                continue
            end
            squared_modulus_error=maximum(abs,abs2.(uu).-(d+1))
            println("CONJUGATION bits=$b squared_modulus_error=$squared_modulus_error"); flush(stdout)
            squared_modulus_error<squared_modulus_tolerance && break
        end
        uu!==nothing && squared_modulus_error<squared_modulus_tolerance || error("Conjugation squared-modulus check failed")
        export_dir="export"
        mkpath(export_dir)
        write(joinpath(export_dir,"units.bin"),vec(uu))
        write(joinpath(export_dir,"labels.bin"),reinterpret(Int64,vec(P)))
        M=standard_map(F)
        matrix_list(A)="["*join([A[1,1],A[1,2],A[2,1],A[2,2]],",")*"]"
        shift=haskey(REPRESENTATIVES,string(d)) ? REPRESENTATIVES[string(d)]["shift"] : -1
        metadata="{\"d\":$d,\"q_b\":$(Int(F.Q.b)),\"preferred_shift\":$shift,\"axes\":[$(join(ords,","))],\"standard_map\":$(matrix_list(M)),\"stabilizers\":[$(join(matrix_list.(S),","))]}"
        write(joinpath(export_dir,"metadata.json"),metadata*"\n")
        println("EXPORT ",abspath(export_dir))
        flush(stdout)
    end
    cd(root)
    run(`uv run $(joinpath(@__DIR__,"reconstruct_large.py")) $(joinpath(work,"export")) --work $(joinpath(work,"reconstruction")) --output $(joinpath(root,"fiducial_$(d)_stark.txt"))`)
end

# Run only as a script, so test_conjugation.jl can include the helpers.
if abspath(PROGRAM_FILE) == @__FILE__
    length(cli_args) in (1,2) || error("usage: stark_large.jl DIMENSION [BITS=2048]")
    run_construct(parse(Int,cli_args[1]),length(cli_args)==2 ? parse(Int,cli_args[2]) : 2048)
    println("PROCESS maxrss=$(Sys.maxrss())")
end
