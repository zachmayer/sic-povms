# Numerical standard-Zauner SIC in d = 2208, form <1,-7,1>, conductor 3.
# Reuses stark.jl without changing the construction of the published smaller d.
# Departures needed at this size:
#   * Factor the Shintani-Faddeev cocycle along the companion-matrix orbit:
#     the direct q-Pochhammer recurrence has 4.57 trillion terms.
#   * FLINT FFT ghost residuals and a frozen Float64 QR solve replace dense
#     high-precision solves; scale corrections to avoid Float64 underflow.
#     Refresh QR after the first 512-bit correction, as the original run did.
#   * Stream selected ghost units rather than storing d^2 BigFloat overlaps.
#   * Search determinant classes of shifts and complete by FFT projections,
#     avoiding the original shift-by-d-by-d arrays and dense Fourier products.
# The seven quotient generators and coordinate map are the choices of the
# original 2208 run, fixed as integer data rather than chosen by Hecke's RNG.
# These are construction parameters, not a saved ghost or fiducial seed.
#
# julia --project=. --threads=1 stark_2208.jl
# Writes fiducial_2208_stark.txt; refine_2208.py writes the 150-digit result.

isempty(ARGS) || error("usage: julia --project=. stark_2208.jl")
include("stark.jl")                 # empty ARGS: definitions only
import Nemo
using Nemo: AcbField, acb_struct

"Cocycle product for L=[7 -1;1 0]; first-label integer shifts have correction 1."
function orbit_shin(F, p1, p2)
    d, β, b = Int(F.d), F.x, Int(-F.Q.b)
    start, val = (p1, p2), one(Complex{BigFloat})
    for _ in 1:F.k
        val *= sigma((p2 * β - p1) / d, β)
        p1, p2 = mod(b * p1 - p2, d), p1
    end
    (p1, p2) == start || error("cocycle orbit did not close")
    return val
end

function cocycle_ghost(F)
    d, A = Int(F.d), F.A
    QA = QuadBin(A[2,1], A[2,2] - A[1,1], -A[1,2])
    ζ = -cispi(BigFloat(1) / d)
    c = -e(-BigFloat(rademacher(A)) / 24) / sqrt(BigFloat(d + 1))
    χ = zeros(Complex{BigFloat}, d, 2)
    χ[1,1] = 1
    for (p1, p2) in Iterators.drop(Iterators.product(0:1, 0:d-1), 1)
        ν = (-1)^(p1 + p2 + p1 * p2) * ζ^(-BigFloat(QA(p1,p2)) / (d * (d - 2))) * c * orbit_shin(F, p1, p2)
        χ[p2+1,p1+1] = ζ^(p1 * p2) * real(ν)
    end
    χ = ifft(χ, 1)
    ψ = sqrt(abs(χ[1,1])) * circshift(cumprod(χ[:,2] ./ χ[:,1]), 1)
    return ψ / ψ[1]
end

"DFT through the pinned Nemo/FLINT C API; retain numerical midpoints."
function flint_dft(x::Vector{Complex{BigFloat}}; inverse=false)
    n, bits = length(x), precision(BigFloat)
    C = AcbField(bits + 32)
    a, b = Nemo.acb_vec(C.(x)), Nemo.acb_vec(n)
    if inverse
        ccall((:acb_dft_inverse, Nemo.libflint), Cvoid,
              (Ptr{acb_struct}, Ptr{acb_struct}, Int, Int), b, a, n, bits + 32)
    else
        ccall((:acb_dft, Nemo.libflint), Cvoid,
              (Ptr{acb_struct}, Ptr{acb_struct}, Int, Int), b, a, n, bits + 32)
    end
    vals = Nemo.array(C, b, n)
    Nemo.acb_vec_clear(a, n)
    Nemo.acb_vec_clear(b, n)
    return [complex(BigFloat(real(z)), BigFloat(imag(z))) for z in vals]
end

"Quartic ghost equations as four circular correlations."
function ghost_residual_fft(z, ns)
    d = 1 + length(z) ÷ 2
    ψ = [complex(one(eltype(z)), zero(eltype(z))); complex.(z[1:d-1], z[d:end])]
    N = sum(real(ψ[mod(-k,d)+1] * conj(ψ[k+1])) for k in 0:d-1)
    rows = Dict{Int,Vector{Complex{eltype(z)}}}()
    for p in unique([radix(n, [2,d,d])[2] for n in ns])
        A = [ψ[mod(-k,d)+1] * conj(ψ[mod(k-p,d)+1]) for k in 0:d-1]
        B = [ψ[k+1] * conj(ψ[mod(p-k,d)+1]) for k in 0:d-1]
        rows[p] = flint_dft(flint_dft(A) .* flint_dft(B); inverse=true)
    end
    return map(ns) do n
        r, p, q = radix(n, [2,d,d])
        v = rows[p][mod(p+q,d)+1] - ((p == 0) + (q == 0)) * N^2 / (d + 1)
        r == 0 ? real(v) : imag(v)
    end
end

function bump_fft(ψ, bits)
    return setprecision(BigFloat, bits) do
        z, ns = BigFloat.(re_im_proj(ψ)), equations(length(ψ))
        J = qr(ghost_jacobian(Float64.(z), ns))
        prev = big(Inf)
        for iter in 0:100
            r = ghost_residual_fft(z, ns)
            res = maximum(abs, r)
            println("ghost $bits bits, step $iter: log2 residual = $(Float64(log2(res)))")
            flush(stdout)
            res < big(2.0)^(-bits + 40) && return re_im_proj(z)
            res < prev / 2 || error("ghost correction stalled")
            scale = big(2.0)^floor(Int, log2(res))
            z -= scale * BigFloat.(J \ Float64.(r / scale))
            # The first correction removes the cocycle's Float64 integration
            # error. Fix this refresh explicitly instead of relying on a restart.
            bits == 512 && iter == 0 && (J = qr(ghost_jacobian(Float64.(z), ns)))
            prev = res
        end
        error("ghost refinement exhausted 100 steps")
    end
end

"Fixed prime-power quotient axes, independently checked modulo the stabilizer."
function arithmetic_orbit(F)
    mats = [[3679 0; 0 3679], [1207 3725; 691 786], [1093 4032; 384 2821],
            [2991 23; 4393 3152], [4111 1614; 2802 2161], [3471 1847; 2569 3152], [673 0; 0 673]]
    ords, dd = [2, 8, 16, 3, 16, 3, 11], 4416
    S = [Matrix{Int}(I, 2, 2)]
    L = Int.(F.L)
    while mod.(S[end] * L, dd) != S[1]
        push!(S, mod.(S[end] * L, dd))
    end
    length(S) == 24 || error("wrong stabilizer order")
    powers = [[Matrix{Int}(I, 2, 2)] for _ in mats]
    for j in eachindex(mats)
        mod.(mats[j] * L - L * mats[j], dd) == zeros(Int, 2, 2) || error("generator does not commute with L")
        for _ in 1:ords[j]-1
            push!(powers[j], mod.(powers[j][end] * mats[j], dd))
        end
        mod.(powers[j][end] * mats[j], dd) in S || error("wrong quotient order")
    end
    P = Array{NTuple{2,Int}}(undef, Tuple(ords))
    seen = falses(dd, dd)
    for index in CartesianIndices(P)
        p = [1, 0]
        for j in eachindex(mats)
            p = mod.(powers[j][index[j]] * p, dd)
        end
        P[index] = Tuple(p)
        for s in S
            a, b = mod.(s * p, dd) .+ 1
            seen[a,b] && error("quotient orbit is not free")
            seen[a,b] = true
        end
    end
    return P, S
end

"Selected real ghost units: one row transform at a time, discarding other labels."
function ghost_units_stream(ψ, P)
    d, φ = length(ψ), circshift(reverse(ψ), 1)
    normfactor = (d + 1) / real(φ'ψ)
    ζ = [(-cispi(big(1.0) / d))^(-j) for j in 0:2d-1]
    rows = [Int[] for _ in 1:d]
    for (i, (a, b)) in enumerate(P)
        push!(rows[mod(a,d)+1], i)
    end
    u = Array{BigFloat}(undef, size(P))
    for a in 0:d-1
        isempty(rows[a+1]) && continue
        vals = [conj(φ[k+1]) * ψ[mod(k-a,d)+1] for k in 0:d-1]
        row = d * flint_dft(vals; inverse=true)
        for i in rows[a+1]
            p, q = P[i]
            u[i] = normfactor * real(ζ[mod(p*q,2d)+1] * row[mod(q,d)+1])
        end
    end
    return u
end

function fft_projector(W, tau)
    d = size(W, 1)
    R = ifft(W .* tau, 2) / (d + 1)
    B = [R[mod(j-k,d)+1,k+1] + (j == k) / d for j in 0:d-1, k in 0:d-1]
    return Hermitian((B + B') / 2)
end

function fft_overlaps(ψ, tau)
    d = length(ψ)
    V = [conj(ψ[mod(k+a,d)+1]) * ψ[k+1] for a in 0:d-1, k in 0:d-1]
    return conj.((d + 1) * d / norm(ψ)^2 .* tau .* ifft(V, 2))
end

function complete_sic(uu, P, S)
    d, dd = 2208, 4416
    M = [338 4023; 3951 5]          # first determinant-one map in the original lexicographic search
    owner, sign = zeros(Int32, d, d), zeros(Int8, d, d)
    for (g, p) in enumerate(P), s in S
        n = M * s
        a = mod(n[1,1] * p[1] + n[1,2] * p[2], dd)
        b = mod(n[2,1] * p[1] + n[2,2] * p[2], dd)
        x, y = mod(a,d), mod(b,d)
        owner[x+1,y+1] = g
        sign[x+1,y+1] = (-1)^div(a*b - x*y, d)
    end
    # Equal determinants differ by a unitary Clifford in the centralizer:
    # rank-one existence is constant on each such class of shifts.
    detseen, shifts = Set{Int}(), CartesianIndex[]
    for (i, p) in zip(CartesianIndices(P), P)
        a, b = mod(p[1] - 7p[2], dd), p[2]
        determinant = mod(a*a + 7a*b + b*b, dd)
        if !(determinant in detseen)
            push!(detseen, determinant)
            push!(shifts, i)
        end
    end
    CI = CartesianIndices(uu)
    at(g, s) = conj(uu[CartesianIndex(mod.(Tuple(CI[g]) .- Tuple(s), size(uu)) .+ 1)])
    A = [zeros(ComplexF64, length(shifts), d) for _ in 1:2]
    for row in 1:2, b in 1:d
        g = owner[row,b]
        g == 0 && continue
        c = sign[row,b] * τ((row-1) * (b-1), d) / (d + 1)
        for (i, s) in enumerate(shifts)
            A[row][i,b] = c * at(g, s)
        end
    end
    D, T = 1 / d .+ real.(ifft(A[1], 2)), ifft(A[2], 2)
    fit = d * vec(sum(abs.(abs2.(T) .- D .* circshift(D, (0,-1))), dims=2))
    tau = τ.(pq(d), d)
    known = owner .!= 0
    known[1,1] = true
    ranked = sortperm(round.(fit; digits=10))  # fixed index breaks numerical ties
    shifts[first(ranked)] == CartesianIndex(1,1,2,1,2,1,4) || error("the original shift is no longer uniquely preferred")
    for index in ranked
        s = shifts[index]
        W = zeros(ComplexF64, d, d)
        for i in eachindex(W)
            owner[i] != 0 && (W[i] = sign[i] * at(owner[i], s))
        end
        λ, V = eigen(fft_projector(W, tau))
        λ[end] > count(known) / d^2 / 2 || continue
        ψ = V[:,end]
        for iter in 1:250
            C = fft_overlaps(ψ, tau)
            err = maximum(abs, abs2.(C)[2:end] .- (d + 1)) / (d + 1)^2
            if err < 1e-13
                println("shift $s, completion $iter: SIC error $err")
                return ψ
            end
            ψ = normalize(fft_projector(ifelse.(known .& (err > 1e-9), W, sqrt(d + 1) * C ./ abs.(C)), tau) * ψ)
        end
    end
    error("no determinant representative produces a SIC")
end

function construct_2208()
    randseed!(0); seed!(0)
    F = AdmissibleTuple(2208, QuadBin(1,-7,1))
    P, S = arithmetic_orbit(F)
    lap("fixed arithmetic quotient")
    ψ = setprecision(() -> cocycle_ghost(F), BigFloat, 128)
    ψ = bump_fft(bump_fft(ψ, 512), 2048)
    lap("2048-bit ghost")
    return setprecision(BigFloat, 2048) do
        u = ghost_units_stream(ψ, P)
        lap("streamed ghost units")
        recipes = normal_recipes(u, 2208)
        x = ghost_basis(u; recipes)[1]
        hprod = hankel_prod(u, x)
        top = max(maximum(abs, hprod), maximum(maximum(abs, X) for X in Iterators.flatten(forms.(x))))
        bits = 2F.h * (ceil(Int, log2(top)) + 40)
        1.5bits + 100 < 2048 || error("insufficient ghost precision")
        lap("normal bases and trace forms")
        H = classfield(F)
        lap("quadratic class field")
        for b in unique(min.(round.(Int, bits .* 1.5 .^ (0:3)), 1948))
            uu = sic_units(conjugator(H, F, b, sum(u .^ 2)), u, recipes)
            maximum(abs, abs2.(uu) .- 2209) < 1e-3 || continue
            lap("conjugated SIC units")
            return complete_sic(uu, P, S)
        end
        error("conjugated units failed modulus test")
    end
end

t = @elapsed ψ = construct_2208()
ψ *= conj(ψ[1]) / abs(ψ[1]) / norm(ψ)
write("fiducial_2208_stark.txt", join(scott(real(z)) * scott(k == 1 ? 0.0 : imag(z)) * "i\n" for (k,z) in enumerate(ψ)))
println("wrote fiducial_2208_stark.txt in $t s")
