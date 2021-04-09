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
    x_expr8 = :(exp(.- sum(abs2.(x .- offset).*scale))) # scale is 1/(2 sigma^2)
    x_expr9 = :(exp(.- sum(abs2.(x .- offset).*scale)) ./ prod(sqrt.(pi ./ scale)))

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
    ]
    return functions
end


# we automatically generate the functions for rr2, rr, ...
# We set the types for the arguments correctly in the default cases
for F in generate_functions_expr() 
    # default functions with offset and scaling behavior
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int};
                           offset=CtrFT,
                           scale=ScaUnit,
                           dims=ntuple(+, N)) where{N, T} 
        offset = get_offset(size, offset)
        scale_init = get_scale(size, scale)
        scale = ntuple(i -> i ∈ dims ? scale_init[i] : zero(scale_init[1]), N)
        IndexFunArray(T, $(F[2]), size) 
    end
    
    # change order of offset and scale
    @eval function $(F[1])(size::NTuple{N, Int}; 
                           offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N}
        T = $DEFAULT_T
        $(F[1])(T, size, scale=scale, offset=offset, dims=dims) 
    end

    # convenient wrapper to provide an array as input
    @eval function $(F[1])(arr::AbstractArray{T, N}; 
                           offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {T, N}
        T2 = default_type(T, $DEFAULT_T)                   
        $(F[1])(T2, size(arr), scale=scale, offset=offset)
    end



    @eval export $(F[1])
end 


