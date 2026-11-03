# The production conjugator against a 512-bit reference. A two-element basis whose
# 1e14-sized terms cancel to 200 (factor 1.1e11; d = 39604's worst unit had 1.8e11).
include(joinpath(@__DIR__,"stark_large.jl"))

setprecision(BigFloat,512) do
    B=big(1000.0)
    z=[B-sqrt(big(5.0)),B+sqrt(big(5.0))]
    G,C=forms(z)
    y=G\z
    trace=[big(10.0)^14,(200-big(10.0)^14*y[1])/y[2]]
    expected=ComplexF64.(hankel_prod(trace,[y]))
    actual=conjugate_units(identity,trace,[(G,C)])
    @assert actual!==nothing
    @assert actual==expected
    println("cancellation factor=",Float64(sum(abs,trace.*y)/abs(expected[1])))
    println("PASS: production conjugator matches the 512-bit reference")
end
