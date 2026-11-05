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
    source=Ref{Any}(copy(trace))
    actual=conjugate_units(identity,source,[(G,C)])
    @assert actual!==nothing
    @assert actual==expected
    @assert source[]===nothing
    mktempdir() do dir
        path=joinpath(dir,"trace.jls")
        save_atomic(path,(1,trace))
        # A rejected squared-modulus result must be reproducible at the next rung.
        source[]=last(deserialize(path))
        @assert conjugate_units(identity,source,[(G,C)])==expected
        @assert source[]===nothing
    end
    println("cancellation factor=",Float64(sum(abs,trace.*y)/abs(expected[1])))
    println("PASS: production conjugator matches the 512-bit reference")
end

setprecision(BigFloat,512) do
    for shape in ((2,3,4),(1,2),(4,))
        T=reshape(sin.(BigFloat.(1:prod(shape))),shape)
        x=[cos.(BigFloat.(1:n)) for n in shape]
        expected=hankel_prod(T,x)
        mktempdir() do dir
            units=joinpath(dir,"units.jls")
            prefix=joinpath(dir,"trace")
            save_atomic(units,T)
            actual=trace_product(copy(T),x,prefix,units)
            @assert actual==expected
            @assert last(deserialize(prefix*".jls"))==expected
            @assert !isfile(units) && !isfile(prefix*"_rolling.jls")
            # Resume after an axis commit, with no raw units left.
            partial=tensor_axis!(copy(T),hankel(x[1]),1)
            save_atomic(prefix*"_rolling.jls",(1,partial))
            @assert trace_product(nothing,x,prefix,units)==expected
        end
        y=[cis.(BigFloat.(1:n)) for n in shape]
        setprecision(BigFloat,LIVE_BITS) do
            low=BigFloat.(T)
            expected=hankel_prod(low,y)
            actual=Array{Union{BigFloat,Complex{BigFloat}}}(low)
            for j in 1:ndims(T)
                tensor_axis!(actual,hankel(y[j]),j,j==1 ? BigFloat : Complex{BigFloat})
            end
            @assert actual==expected
        end
    end
    println("PASS: tensor axes match hankel_prod exactly; rolling checkpoints resume")
end
