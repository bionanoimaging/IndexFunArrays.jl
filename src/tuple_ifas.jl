function generate_tuple_functions_expr()
    x_expr = :(scale .* (x .- offset))

    functions = [
        (:(idx),  :(x -> T.($x_expr))), # returns a tuple
    ]
    return functions
end

# we automatically generate tuple functions
# We set the types for the arguments correctly in the default cases
for F in generate_tuple_functions_expr() 
    # default functions with offset and scaling behavior
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int};
                           offset=CtrFT,
                           scale=ScaUnit, dims=ntuple(+, N)) where{N, T} 
        T2 = typeof(T.(size))
        offset = get_offset(size, offset)
        scale_init = get_scale(size, scale)
        scale = ntuple(i -> i ∈ dims ? scale_init[i] : zero(scale_init[1]), N)
        IndexFunArray(T2, $(F[2]), size) 
    end
    
    # change order of offset and scale
    @eval function $(F[1])(size::NTuple{N, Int}; 
                           offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N}
        T2 = $DEFAULT_T
        $(F[1])(T2, size, scale=scale, offset=offset, dims=dims) 
    end

    # convenient wrapper to provide an array as input
    @eval function $(F[1])(arr::AbstractArray{T, N}; 
                           offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {T, N}
        T2 = default_type(T, $DEFAULT_T)                   
        $(F[1])(T2, size(arr), scale=scale, offset=offset, dims=dims)
    end

    @eval export $(F[1])
end 


