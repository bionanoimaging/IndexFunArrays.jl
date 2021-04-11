apply_covariance = (x,M) -> sum(dot.(x,M*x))

# List of functions and names we want to predefine
function generate_functions_expr()
    # offset and scale is already wrapped in the generator function
    x_expr = :(scale .* (x .- offset))
    x_expr1 = :(scale[1] .* (x[1] .- offset[1]))
    x_expr2 = :(scale[2] .* (x[2] .- offset[2]))
    x_expr3 = :(scale[3] .* (x[3] .- offset[3]))
    x_expr4 = :(scale[4] .* (x[4] .- offset[4]))
    x_expr5 = :(scale[5] .* (x[5] .- offset[5]))
    x_expr6 = :(prod(x .== offset))
    x_expr7 = :(cis(dot((x .- offset), scale)))
    x_expr8 = :(exp(.- norm((x .- offset).*scale))) # scale is 1/(2 sigma^2)
    x_expr9 = :(exp(.- norm((x .- offset).*scale)) ./ prod(sqrt.(pi ./ abs2.(scale))))
    x_expr10 = :(exp(.- apply_covariance((x .- offset), scale) ))
    x_expr11 = :(exp(.- apply_covariance((x .- offset), scale)) ./ prod(sqrt.(pi ./ abs2(scale))))

    functions = [
        (:(rr2), :(x -> T(sum(abs2.($x_expr))))),
        (:(rr),  :(x -> T(sqrt.(sum(abs2.($x_expr)))))),
        (:(xx),  :(x -> T($x_expr1))),
        (:(yy),  :(x -> T($x_expr2))),
        (:(zz),  :(x -> T($x_expr3))),
        (:(ee),  :(x -> T($x_expr4))),
        (:(tt),  :(x -> T($x_expr5))),
        (:(delta),  :(x -> T($x_expr6))),
        (:(phiphi), :(x -> T(atan.($x_expr2, $x_expr1)))),  # this is the arcus tangens of y/x yielding a spiral phase ramp
        (:(exp_is),  :(x -> T($x_expr7))),  # exp(2pi i s (x-o)) # by modifying s, this becomes exp(i kx)
        (:(exp_sqr),  :(x -> T($x_expr8))),  # maximum-normalized Gaussian
        (:(exp_sqr_norm),  :(x -> T($x_expr9))),  # integral-normalized Gaussian (over infinite ROI)
        (:(exp_sqr_cov),  :(x -> T($x_expr10))),  # Gaussian based on the covariance matrix
        (:(exp_sqr_cov_norm),  :(x -> T($x_expr11))),  # normalized Gaussian based on the covariance matrix
    ]
    return functions
end

function wrap_fkt(f)  # encapsulates the expression in a proper function with index, offset and scale argument
    return (idx, offset, scale) -> f(idx)  # Not that f contains offset and scale
end

function wrap_zipped(g) # encapsulates a wrapped function as called with a sinle tuple of offset and scale
    return (idx, mytuple) -> g(idx,mytuple...)  
end

curry(f, x) = (xs...) -> f(x, xs...)   # just a shorthand to remove x

mat_to_tvec = (v) -> [Tuple(v[:,n]) for n in 1:size(v,2)] # converts a 2d matrix to a Vector of Tuples

# we automatically generate the functions for rr2, rr, ...
# We set the types for the arguments correctly in the default cases
for F in generate_functions_expr() 

    # default functions with offset and scaling behavior
    @eval function $(Symbol(:_, F[1]))(::Type{T}, size::NTuple{N, Int},
                           offset,
                           scale,
                           dims, 
                           accumulator) where{N, T} 
        offset = get_offset(size, offset)
        scale_init = get_scale(size, scale)
        scale = ntuple(i -> i ∈ dims ? scale_init[i] : zero(scale_init[1]), N)

        #print("Dispatch T0\n")
        IndexFunArray(T, $(F[2]), size) 
    end

    #=
    # a version (by Felix) that allows a Vector as input
    @eval function $(Symbol(F[1], :(base)))(::Type{T}, size::NTuple{N, Int},
                           offset::Vector,
                           scale::Vector,
                           dims) where{N, T}
        offsets_a = get_offset.(Ref(size), offset)
        scales_a = get_scale.(Ref(size), scale)
        
        g(x) = begin
            res = zero(T)
            for (offset, scale) in zip(offsets_a, scales_a)
                res += $(F[2])(x)
            end
            return res
        end
        IndexFunArray(T, g, size) 
    end
    =#
    @eval function $(Symbol(:_, F[1]))(::Type{T}, size::NTuple{N, Int},
        offset::IterType,   # a version that supports an iterable collection of tuples for offset and scale 
        scale::IterType,
        dims=ntuple(+, N),
        accumulator = sum) where{N, T} 
        #print("Dispatch TupleTuple\n")
        w = (x, offset, scale) -> $(F[2])(x) # adds offset and scale as requires parameters to the expression
        # the line below makes a function only depending on the index position but iterating over both, offsets and scales
        fkt = (idx) -> accumulator(curry(wrap_zipped(w),idx),zip(offset,scale))
        IndexFunArray(T, fkt, size) 
    end

    @eval function  $(Symbol(:_, F[1]))(::Type{T}, size::NTuple{N, Int},
        offset::IterType,   # a version that supports an iterable collection of tuples for offset and scale 
        scale=ScaUnit,
        dims=ntuple(+, N),
        accumulator = sum) where{N, M, T} 

       # print("Dispatch T1\n")
        $(Symbol(:_, F[1]))(T, size, offset, repeated(scale), dims, accumulator)
    end

    @eval function  $(Symbol(:_, F[1]))(::Type{T}, size::NTuple{N, Int},
        offset,   # a version that supports an iterable collection of tuples for offset and scale 
        scale::IterType,
        dims=ntuple(+, N);
        accumulator = sum) where{N, M, T} 

        #print("Dispatch T2\n")
        $(Symbol(:_, F[1]))(T, size, repeated(offset), scale, dims, accumulator)
    end

    # default functions with offset and scaling behavior
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int};
                           offset=CtrFT,
                           scale=ScaUnit,
                           dims=ntuple(+, N), accumulator=sum) where{N, T} 
        # $(F[1])(T, size, offset, scale, dims) 
        #print("Dispatch T3\n")
        $(Symbol(:_, F[1]))(T, size, offset, scale, dims, accumulator) 
    end
    
    # change order of offset and scale
    @eval function $(F[1])(size::NTuple{N, Int}; 
                           offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum) where {N}
        T = $DEFAULT_T
        $(F[1])(T, size, scale=scale, offset=offset, dims=dims, accumulator=accumulator) 
    end

    # convenient wrapper to provide an array as input
    @eval function $(F[1])(arr::AbstractArray{T, N}; 
                           offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum) where {T, N}
        T2 = default_type(T, $DEFAULT_T)                   
        $(F[1])(T2, size(arr), scale=scale, offset=offset, accumulator=accumulator)
    end

    @eval export $(F[1])
end 


