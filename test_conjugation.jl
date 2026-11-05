# Exact record storage, ordered reductions, threaded axes and conjugation precision.
include(joinpath(@__DIR__,"stark_large.jl"))

function packed_copy(path,A)
    bits=precision(real(first(A)))
    T=packed(path,size(A),bits,eltype(A);create=true)
    open(path,"w") do io
        words=zeros(UInt64,record_words(T))
        for value in A
            pack!(words,1,value)
            write(io,words)
        end
        sync_file(io)
    end
    T
end
read_all(T)=reading(collect,T)

for bits in (512,2048,2304)
    setprecision(BigFloat,bits) do
        values=[BigFloat(0),-BigFloat(0),BigFloat(Inf),BigFloat(-Inf),BigFloat(NaN)]
        append!(values,[ldexp((-1)^i*rand(BigFloat),17i-1100) for i in 1:128])
        x=buffer(BigFloat,bits)
        words=zeros(UInt64,2+cld(bits,64))
        for value in values
            pack!(words,1,value)
            unpack!(x,words,1)
            @assert isequal(x,value) && precision(x)==bits
            @assert x.sign==value.sign && x.exp==value.exp
            again=similar(words)
            pack!(again,1,x)
            @assert again==words
        end
        mktempdir() do dir
            for A in (values,complex.(values,reverse(values)))
                T=packed_copy(joinpath(dir,"records.bin"),A)
                @assert isequal(read_all(T),A)
            end
        end
    end
end
println("PASS: real and complex records round-trip at 512, 2048 and 2304 bits")

setprecision(BigFloat,512) do
    mktempdir() do dir
        # Cross both Julia's pairwise-reduction boundaries and the file read cache.
        A=sin.(BigFloat.(1:20001))
        T=packed_copy(joinpath(dir,"stats.bin"),A)
        reading(T) do R
            @assert sum(abs2,R)==sum(abs2,A)
            @assert maximum(abs,R)==maximum(abs,A)
        end
        A=reshape(A[1:30],2,3,5)
        T=packed_copy(joinpath(dir,"marginals.bin"),A)
        for kind in (:linear,:shifted_quad,:shifted_quad_sq), j in 1:ndims(A)
            recipe=Zauner.NormalBasisRecipe(kind,kind===:linear ? nothing : (1,2,4))
            expected=zeros(BigFloat,size(A,j))
            for ci in CartesianIndices(A)
                value=if kind===:linear
                    A[ci]
                else
                    shifted=CartesianIndex(mod.(Tuple(ci).-recipe.shift.-1,size(A)).+1)
                    A[ci]+A[ci]*A[shifted]
                end
                kind===:shifted_quad_sq && (value=value^2)
                expected[ci[j]]+=value
            end
            @assert basis_marginal(T,j,recipe)==expected
        end
    end
end
println("PASS: pairwise statistics and all marginal recipes retain their arithmetic")

setprecision(BigFloat,512) do
    # Last shape crosses the 1024-slice tile and has a partial final tile.
    for shape in ((2,3,4),(1,2),(4,),(3,5,7,11,2,7))
        A=reshape(sin.(BigFloat.(1:prod(shape))),shape)
        x=[cos.(BigFloat.(1:n)) for n in shape]
        expected=hankel_prod(A,x)
        mktempdir() do dir
            units=joinpath(dir,"units.bin")
            prefix=joinpath(dir,"trace")
            T=packed_copy(units,A)
            # An incomplete successor has no committed state and must be overwritten.
            write(prefix*"_axis1.bin",UInt8[0,1,2])
            actual=trace_product(T,x,prefix,units)
            @assert read_all(actual)==expected
            @assert read_all(last(deserialize(prefix*".jls")))==expected
            @assert !isfile(units)
            @assert read_all(trace_product(nothing,x,prefix,units))==expected

            # A committed first axis is sufficient after the original units disappear.
            prefix=joinpath(dir,"resume")
            T=packed_copy(units,A)
            partial=packed(prefix*"_axis1.bin",shape,512;create=true)
            tensor_axis!(partial,T,hankel(x[1]),1)
            save_atomic(prefix*".jls",(1,partial))
            rm(units)
            write(prefix*"_axis2.bin",UInt8[3,4])
            @assert read_all(trace_product(nothing,x,prefix,units))==expected

            # First axis is real input; every subsequent axis is complex input.
            y=[cis.(BigFloat.(1:n)) for n in shape]
            setprecision(BigFloat,LIVE_BITS) do
                low=map(v->BigFloat(v;precision=LIVE_BITS),A)
                expected_live=hankel_prod(low,y)
                T=packed_copy(units,low)
                for j in 1:ndims(A)
                    out=packed(joinpath(dir,"complex$j.bin"),shape,LIVE_BITS,Complex{BigFloat};create=true)
                    tensor_axis!(out,T,hankel(y[j]),j)
                    T=out
                end
                @assert read_all(T)==expected_live
            end
        end
    end
end
println("PASS: threaded real/complex axes match hankel_prod exactly; checkpoints resume")

# Terms near 1e14 cancel to 200 (factor 1.1e11; d=39604 needed 1.8e11).
setprecision(BigFloat,512) do
    B=big(1000.0)
    z=[B-sqrt(big(5.0)),B+sqrt(big(5.0))]
    G,C=forms(z)
    y=G\z
    trace=[big(10.0)^14,(200-big(10.0)^14*y[1])/y[2]]
    expected=ComplexF64.(hankel_prod(trace,[y]))
    mktempdir() do dir
        cd(dir) do
            source=packed_copy("trace.bin",trace)
            original_hash=sha256(read(source.path))
            save_atomic("trace.jls",(1,source))
            actual=conjugate_units(identity,source,[(G,C)])
            @assert actual!==nothing && ComplexF64.(read_all(actual))==expected
            # A failed modulus guard retries from the untouched committed trace.
            source=last(deserialize("trace.jls"))
            again=conjugate_units(identity,source,[(G,C)])
            @assert ComplexF64.(read_all(again))==expected
            @assert sha256(read(source.path))==original_hash
        end
    end
    println("cancellation factor=",Float64(sum(abs,trace.*y)/abs(expected[1])))
end
println("PASS: production conjugator matches the 512-bit reference and preserves trace")

setprecision(BigFloat,512) do
    mktempdir() do dir
        d=70 # Three independently committed chunks, including a partial final chunk.
        ψ=[complex(BigFloat(k+1)/d,BigFloat(mod(k,3))/d) for k in 0:d-1]
        P=reshape([(a,b) for a in 0:d-1 for b in 0:1],d,2)
        baseline=ghost_units_stream(ψ,P,joinpath(dir,"baseline.bin"))
        expected=read(baseline.path)
        path=joinpath(dir,"resumed.bin")
        cp(baseline.path,path)
        save_atomic(path*".rows",[true,false,false])
        words=zeros(UInt64,record_words(baseline))
        open(path,"r+") do io
            # Committed first-chunk records survive; unfinished chunks may contain anything.
            for (i,(a,b)) in enumerate(P)
                a<32 && continue
                seek(io,8length(words)*(i-1))
                write(io,words)
            end
            sync_file(io)
        end
        @assert read(path)!=expected
        ghost_units_stream(ψ,P,path)
        @assert read(path)==expected
        @assert all(deserialize(path*".rows"))
        # Completed markers skip every row even when the supplied ghost would differ.
        ghost_units_stream(fill(one(eltype(ψ)),d),P,path)
        @assert read(path)==expected
        # A crash before the initial completion map committed restarts the raw file.
        initial=joinpath(dir,"uncommitted.bin")
        write(initial,UInt8[1,2,3])
        ghost_units_stream(ψ,P,initial)
        @assert read(initial)==expected
    end
end
println("PASS: unit row markers resume after incomplete writes and skip committed rows")

