# General companion-form Stark/cocycle construction, followed by streamed FFT
# reconstruction. Run: julia --project=. --threads=1 stark_large.jl DIMENSION
# Then: uv run refine_large.py fiducial_DIMENSION_stark.txt --output FILE
#
# Canonical rule for dimensions without a representative pin:
# * Choose the least b>=3 with b^2-4=D*q^2, q|f, and standard-Zauner symmetry,
#   where (d+1)(d-3)=D*f^2. The form is <1,-b,1>; b=d-1 always qualifies
#   for d>=4. Its conductor can exceed stark.jl's choice (e.g. d=206 or 199).
# * Work modulo d for odd d, 2d for even d. Quotient by the FULL overlap
#   stabilizer, including a negative-norm unit when present (no conjugation).
# * Prime-power quotient axes: increasing prime, decreasing exponent; for each
#   axis choose the first (a,y) in lexicographic residue order for a*I+y*L
#   having that order both before and after quotienting by previous axes.
# * Normal-basis recipes: linear, shifted quadratic, squared quadratic, with
#   shifts enumerated by Zauner.radix, last index fastest; Fourier coefficients
#   must exceed 2^-128.
# * Conjugation maximizes sigma(sum(u^2)), as stark.jl. Each multiplication
#   matrix uses the root with smallest (Re,Im), rounded to 10 decimal places;
#   a collision between distinct roots fails rather than choosing an eigenbasis.
# * First trace-minus-one power of L, then first lexicographic (a,b) with
#   det([v,Kv])=1, defines the standard-Zauner map. Reconstruction ranks fits
#   rounded to 10 decimal places, then column-major shift index; the first
#   converging determinant class wins. Power iteration starts from e_0.
# * representatives.toml pins representatives from an unpublished prototype run
#   so its Mac runs and the Linux rebuilds compare byte for byte. Their quotient
#   coordinates came from Hecke, not an arithmetic rule derived from d. The table
#   pins integer generators and determinant shifts tied at roundoff; construction
#   retains the prototype's start rule. It also pins the separate Fibonacci
#   search's random seed, never a saved numerical vector. Every published large-d
#   file uses pinned data; test_large.py tests the canonical rule only at small d.
#
# The cocycle is factored along L's label orbit; ghost residuals and unit rows
# use FLINT FFTs. A frozen Float64 QR corrects the full-precision ghost.
# Normal-basis traces stream along tensor axes. Python reconstruction stores
# owner/sign maps on disk, avoiding complex d-by-d arrays. Final refinement is
# full-space at 640 bits, without a Float64 invariant-space constraint.
# Checkpoints bind source, dimension and precision; use a fresh directory after
# changing precision, source, Julia or platform. Rebuild evidence is in RESULTS.md.

const cli_args = copy(ARGS)
empty!(ARGS)
include("stark.jl")
using TOML, SHA

const REPRESENTATIVES = TOML.parsefile(joinpath(@__DIR__, "representatives.toml"))

function companion_form(d)
    d >= 4 || error("dimension must be at least four")
    D, f = coredisc((big(d) + 1) * (d - 3))
    for b in 3:d-1
        delta = big(b)^2 - 4
        delta % D == 0 || continue
        q = isqrt(div(delta, D))
        q^2 * D == delta && f % q == 0 || continue
        mod(d,9) == 3 && mod(div(f,q),3) == 0 && continue
        return AdmissibleTuple(d, QuadBin(1,-b,1))
    end
    error("no supported standard-Zauner companion form for d=$d")
end

function check_configuration(d,bits)
    bits >= 512 || error("construction requires at least 512 bits")
    names = ("stark.jl", "stark_large.jl", "representatives.toml", "reconstruct_large.py", "search.py", "Project.toml", "Manifest.toml")
    fingerprint = bytes2hex(sha256(join(read(joinpath(@__DIR__,name),String) for name in names)))
    value = "$d $bits $VERSION $(Sys.MACHINE) $fingerprint\n"
    path = "configuration.txt"
    if isfile(path)
        read(path,String) == value || error("checkpoint source/platform/settings changed; use a fresh LARGE_D_WORK")
    else
        write(path,value)
    end
end

function arithmetic_generators(F)
    d,dd,b=Int(F.d),2^iseven(F.d)*Int(F.d),Int(-F.Q.b)
    if haskey(REPRESENTATIVES,string(d))
        pin=REPRESENTATIVES[string(d)]
        return [permutedims(reshape(Int.(v),2,2)) for v in pin["matrices"]],Int.(pin["orders"])
    end
    beta=(-F.Q.b+F.q*F.a)/2
    O=Order(F.K,[one(F.a),beta])
    R,_=quo(O,ideal(O,dd))
    UR,ur=unit_group(R)
    anti,A=is_antiunitary_with_generator(F)
    generator = if anti
        y=A[2,1]
        A[1,1]-y*F.L[1,1]+y*beta
    else
        beta
    end
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

function ghost_nu(F,p1,p2; reduced=true)
    d,A,β = Int(F.d),F.A,F.x
    QA = QuadBin(A[2,1],A[2,2]-A[1,1],-A[1,2])
    ζ = -cispi(BigFloat(1)/d)
    c = -e(-BigFloat(rademacher(A))/24)/sqrt(BigFloat(d+1))
    s = iseven(d)*(p1+p2+p1*p2)
    if reduced
        sh = orbit_shin(F,p1,p2)
    else
        v = Zauner._get_periods(A,β)
        k = Int(divexact(A[1,2]*p2-(A[2,2]-1)*p1,d))
        sh = sigma((p2*β-p1)/d+k,β)
        for i in 2:length(v)-2
            sh *= sigma((v[1]*p2-v[2]*p1)/(d*v[i+1]),v[i]/v[i+1])
        end
    end
    (-1)^s * ζ^(-BigFloat(QA(p1,p2))/(d*(d-2))) * c * sh
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

import Nemo
using Nemo: AcbField, acb_struct

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


using GenericFFT: fft, ifft

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
function bump_float(ψ,bits;name="ghost",maxsteps=100)
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
        for iter in 0:maxsteps
            t = @elapsed r=ghost_residual_fft(z,ns)
            res = maximum(abs,r)
            println("bump iter=$iter bits=$bits residual=$(Float64(log2(res))) log2 seconds=$t")
            flush(stdout)
            path="data/$(name)_$(bits)_progress.jls"
            serialize(path*".part",re_im_proj(z))
            mv(path*".part",path;force=true)
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
function ghost_units_stream(ψ,P;checkpoint_dir="",chunk_rows=32,worker=0,workers=1,assemble=true)
    d=length(ψ)
    φ=circshift(reverse(ψ),1)
    normfactor=(d+1)/real(φ'ψ)
    ζ=[(-cispi(big(1.0)/d))^(-j) for j in 0:2d-1]
    rows=[Int[] for _ in 1:d]
    for (i,(a,b)) in enumerate(P)
        push!(rows[mod(a,d)+1],i)
    end
    isempty(checkpoint_dir) || mkpath(checkpoint_dir)
    0<=worker<workers || error("invalid row worker")
    workers>1 && isempty(checkpoint_dir) && error("parallel row workers need checkpoints")
    u=isempty(checkpoint_dir) ? Array{BigFloat}(undef,size(P)) : nothing
    for firstrow in worker*chunk_rows:workers*chunk_rows:d-1
        lastrow=min(d-1,firstrow+chunk_rows-1)
        path=joinpath(checkpoint_dir,"rows_$(firstrow)_$(lastrow).jls")
        if !isempty(checkpoint_dir) && isfile(path)
            continue
        end
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
        if !isempty(checkpoint_dir)
            serialize(path*".part",(indices,values))
            mv(path*".part",path;force=true)
        else
            u[indices]=values
        end
        println("ghost overlap rows=$firstrow:$lastrow/$d")
        flush(stdout)
    end
    assemble || return nothing
    if !isempty(checkpoint_dir)
        u=Array{BigFloat}(undef,size(P))
        for firstrow in 0:chunk_rows:d-1
            lastrow=min(d-1,firstrow+chunk_rows-1)
            indices,values=deserialize(joinpath(checkpoint_dir,"rows_$(firstrow)_$(lastrow).jls"))
            u[indices]=values
        end
    end
    return u
end


"Marginal normal-basis sums without materializing a shifted Galois tensor."
function basis_marginal(u,j,recipe)
    shape=size(u)
    stride=prod(shape[1:j-1];init=1)
    result=zeros(BigFloat,shape[j])
    linear=recipe.kind===:linear
    zero_shift=linear || all(iszero,recipe.shift)
    squared=recipe.kind===:shifted_quad_sq
    if zero_shift
        for i in eachindex(u)
            v=linear ? u[i] : u[i]+u[i]^2
            squared && (v=v^2)
            result[mod((i-1)÷stride,shape[j])+1]+=v
        end
    else
        LI=LinearIndices(u)
        for ci in CartesianIndices(u)
            shifted=CartesianIndex(mod.(Tuple(ci).-recipe.shift.-1,shape).+1)
            v=u[ci]+u[ci]*u[LI[shifted]]
            squared && (v=v^2)
            result[ci[j]]+=v
        end
    end
    result
end

function normal_recipes(u,d)
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

"Reuse trace forms; conjugate directly without a third large BigFloat tensor."
function conjugate_units(σ,trace,forms_x)
    y=Vector{Vector{ComplexF64}}()
    for pair in forms_x
        G,C=(X->σ.(X)).(pair)
        # The caller retries the integer relation at higher precision. Returning
        # nothing records an unsuccessful precision attempt, never a candidate.
        issuccess(lu(G;check=false)) || return nothing
        A=Float64.(C/G)
        all(isfinite,A) || return nothing
        λ,V=eigen(A)
        keys=[(round(real(z);digits=10),round(imag(z);digits=10)) for z in λ]
        allunique(keys) || error("multiplication roots tie at 10 decimals")
        index=argmin(keys)
        push!(y,ComplexF64.(G)\(V[:,index]*λ[index]/V[1,index]))
    end
    transformed=map(v->ComplexF64(σ(v)),trace)
    hankel_prod(transformed,y)
end

function generator_orbit(F,mats,ords)
    d=Int(F.d)
    multiplier=iseven(d) ? 2 : 1
    d>1 || throw(ArgumentError("dimension must exceed one"))
    d<=typemax(Int)÷multiplier || throw(OverflowError("residue modulus exceeds Int"))
    dd=multiplier*d
    # Both matrix products and scalar updates sum two products of residues.
    # At d=39604, dd=79208 and the largest sum is 12_547_497_698: safe in Int64.
    2big(dd-1)^2<=typemax(Int) || throw(OverflowError("two residue products exceed Int"))
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
    dd=2^iseven(F.d)*Int(F.d)
    stabs=[Int.(mod.(S,dd)) for S in stabilizer_elements(F)]
    for (M,n) in zip(mats,ords)
        gcd(M[1,1]*M[2,2]-M[1,2]*M[2,1],dd)==1 || error("nonunit quotient generator")
        mod.(M*Int.(F.L)-Int.(F.L)*M,dd)==zeros(Int,2,2) || error("noncommuting quotient generator")
        @assert mod.(power_mod_matrix(M,n,dd),dd) in stabs
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

function power_mod_matrix(M,n,dd)
    R=Matrix{Int}(I,2,2)
    while n>0
        isodd(n) && (R=mod.(R*M,dd))
        M=mod.(M*M,dd)
        n>>=1
    end
    R
end

using GenericFFT: fft, ifft

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
        bump_float(start,bits;name="ghost",maxsteps=max(100,bits÷20))
    end
    setprecision(BigFloat,bits) do
        u=stage("units_$bits") do
            ghost_units_stream(ψ,P;checkpoint_dir="units_$(bits)_rows")
        end
        recipes=stage("recipes_$bits") do
            normal_recipes(u,d)
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
        H=stage("classfield") do
            classfield(F)
        end
        uu=nothing
        for b in unique(min.(round.(Int,relation_bits.*1.5.^(0:3)),bits-100))
            uu=stage("liveunits_$(bits)_$b") do
                conjugate_units(conjugator(H,F,b,sum(abs2,u)),trace,forms_x)
            end
            if uu===nothing
                println("CONJUGATION bits=$b insufficient relation precision; advancing precision")
                flush(stdout)
                continue
            end
            err=maximum(abs,abs2.(uu).-(d+1))
            println("CONJUGATION bits=$b modulus_error=$err"); flush(stdout)
            err<1e-3 && break
        end
        uu!==nothing && maximum(abs,abs2.(uu).-(d+1))<1e-3 || error("Conjugation modulus check failed")
        export_dir="export"
        mkpath(export_dir)
        write(joinpath(export_dir,"units.bin"),vec(uu))
        write(joinpath(export_dir,"labels.bin"),reinterpret(Int64,vec(P)))
        M=standard_map(F)
        matrix_list(A)="["*join([A[1,1],A[1,2],A[2,1],A[2,2]],",")*"]"
        legacy=haskey(REPRESENTATIVES,string(d))
        shift=legacy ? REPRESENTATIVES[string(d)]["shift"] : -1
        metadata="{\"d\":$d,\"q_b\":$(Int(F.Q.b)),\"legacy_start\":$legacy,\"preferred_shift\":$shift,\"axes\":[$(join(ords,","))],\"standard_map\":$(matrix_list(M)),\"stabilizers\":[$(join(matrix_list.(S),","))]}"
        write(joinpath(export_dir,"metadata.json"),metadata*"\n")
        println("EXPORT ",abspath(export_dir))
        flush(stdout)
    end
    cd(root)
    python = get(ENV,"UV","uv")
    run(`$python run $(joinpath(@__DIR__,"reconstruct_large.py")) $(joinpath(work,"export")) --work $(joinpath(work,"reconstruction")) --output $(joinpath(root,"fiducial_$(d)_stark.txt"))`)
end

length(cli_args) in (1,2) || error("usage: stark_large.jl DIMENSION [BITS=2048]")
run_construct(parse(Int,cli_args[1]),length(cli_args)==2 ? parse(Int,cli_args[2]) : 2048)
println("PROCESS maxrss=$(Sys.maxrss())")
