apply_covariance = (x,M) -> sum(dot.(x,M*x))

optional_posZ(x::NTuple{1,T}, offset::NTuple{1,T}) where {T} = one(T)
optional_posZ(x::NTuple{2,T}, offset::NTuple{2,T}) where {T} = one(T)
optional_posZ(x::NTuple{N,T}, offset::NTuple{N,T}) where {T,N} = x[3]-offset[3]

abs2_scale(xo, scale) = ifelse(isinf(scale), abs2(xo) / eps(typeof(scale)), abs2(xo)*scale)

# List of functions and names we want to predefine
function generate_functions_expr()
    # offset and scale is already wrapped in the generator function
    x_expr = :(scale .* (x .- offset))
    x_expr1 = :(scale[1] .* (x[1] .- offset[1]))
    x_expr2 = :(scale[2] .* (x[2] .- offset[2]))
    x_expr3 = :(scale[3] .* (x[3] .- offset[3]))
    x_expr4 = :(scale[4] .* (x[4] .- offset[4]))
    x_expr5 = :(scale[5] .* (x[5] .- offset[5]))
    x_expr14 = :(maximum(abs.(scale .* (x .- offset)))) # useful for rectangular things
    x_expr15 = :(minimum(abs.(scale .* (x .- offset))))
    x_expr6 = :(prod(x .== offset))
    x_expr7 = :(cis(dot((x .- offset), scale)))
    # x_expr8 = :(exp(.- sum(abs2.(x .- offset).*scale))) # scale is 1/(2 sigma^2)
    # x_expr9 = :(exp(.- sum(abs2.(x .- offset).*scale)) ./ prod(sqrt.(pi ./ scale)))
    x_expr8 = :(exp(.- sum(abs2_scale.(x .- offset, scale)))) # scale is 1/(2 sigma^2)
    x_expr9 = :(exp(.- sum(abs2_scale.(x .- offset, scale))) ./ prod(ifelse.(isinf.(scale), 1, sqrt.(pi ./ scale))))
    x_expr10 = :(exp(.- apply_covariance((x .- offset), scale) ))
    x_expr11 = :(exp(.- apply_covariance((x .- offset), scale)) ./ prod(sqrt.(pi ./ abs2(scale))))
    x_expr12 = :(abs2.(scale[1] .* (x[1] .- offset[1])) .+ abs2.(scale[2] .* (x[2] .- offset[2])))
    x_expr13 = :((x[1] .- offset[1]).* scale[1] .+ (x[2] .- offset[2]) .* scale[2])

    functions = [
        (:(rr2), :(x -> T(sum(abs2.($x_expr))))),
        (:(rr),  :(x -> T(sqrt.(sum(abs2.($x_expr)))))),
        (:(xx),  :(x -> T($x_expr1))),
        (:(yy),  :(x -> T($x_expr2))),
        (:(zz),  :(x -> T($x_expr3))),
        (:(ee),  :(x -> T($x_expr4))),
        (:(tt),  :(x -> T($x_expr5))),
        (:(box1),  :(x -> T(maximum(abs.($x_expr)) .<= 1) ) ),
        (:(disc1),  :(x -> T(sum(abs2.($x_expr)) .<= 1))),
        (:(idx_min),  :(x -> T($x_expr15))),
        (:(idx_max),  :(x -> T($x_expr14))),
        (:(delta),  :(x -> T($x_expr6))),
        (:(phiphi), :(x -> T(atan.($x_expr2, $x_expr1)))),  # this is the arcus tangens of y/x yielding a spiral phase ramp
        (:(phase_kz), :(x -> T(optional_posZ(x,offset) .* sqrt.(max.(one(T) .- $x_expr12, zero(T)))))),  # can be used for constucting a free-space propagator in optics
        (:(phase_kxy), :(x -> T($x_expr13))),  # useful for xy shifting
        (:(exp_is),  :(x -> T($x_expr7))),  # exp(2pi i s (x-o)) # by modifying s, this becomes exp(i kx)
        (:(exp_sqr),  :(x -> T($x_expr8))),  # maximum-normalized Gaussian
        (:(exp_sqr_norm),  :(x -> T($x_expr9))),  # integral-normalized Gaussian (over infinite ROI)
        (:(exp_sqr_cov),  :(x -> T($x_expr10))),  # Gaussian based on the covariance matrix
        (:(exp_sqr_cov_norm),  :(x -> T($x_expr11))),  # normalized Gaussian based on the covariance matrix
    ]
    return functions
end

function wrap_fkt(f)  # encapsulates the expression in a proper function with index, offset and scale argument
    return (idx, offset, scale, weight) -> f(idx)  # Not that f contains offset and scale
end

function wrap_zipped(g) # encapsulates a wrapped function as called with a sinle tuple of offset and scale
    return (idx, mytuple) -> g(idx,mytuple...)  
end

curry(f, x) = (xs...) -> f(x, xs...)   # just a shorthand to remove x

## These functions ensure that also numbers can be iterated and zipped
function cast_iter(vals::Matrix)
    Tuple(Tuple(vals[:,n]) for n in 1:size(vals,2))
    # Tuple(vals[:,axes(vals, 2)])
end

function cast_iter(vals::Vector{<:Number})
    Tuple(vals)   # makes the matrix iterable, interpreting it as a series of vectors
end

function cast_iter(vals::IterType)
    vals
end

function cast_iter(vals)
    repeated(vals) # during the zip operation this type always yields the same value
end

function cast_number_iter(vals::Vector{<:Number})
    Tuple(vals)   # makes the matrix iterable, interpreting it as a series of vectors
end

function cast_number_iter(vals::NTuple{N,<:Number} where N)
    vals
end

function cast_number_iter(vals::Number)
    repeated(vals) # during the zip operation this type always yields the same value
end


function optional_mat_to_iter(vals)  # only for matrices
    vals
end

function optional_mat_to_iter(vals::Matrix)
    cast_iter(vals)
end

function apply_dims(scale2, dims, N)  # replaces scale entries not in dims with zeros
    ntuple(i -> i ∈ dims ? scale2[i] : zero(eltype(scale2)), N)  # zero(scale2[1])
end

function apply_dims(scales::IterType, dims, N)  # replaces scale entries not in dims with zeros
    Tuple(ntuple(i -> i ∈ dims ? scale3[i] : zero(eltype(scale3)), N)  for scale3 in  scales)
end

# we automatically generate the functions for rr2, rr, ...
# We set the types for the arguments correctly in the default cases
for F in generate_functions_expr() 

    # default functions with offset and scaling behavior. This version allows no list of numbers or tuples
    @eval function $(Symbol(:_, F[1]))(::Type{T}, size::NTuple{N, Int},
                           offset::Union{Type{<:Ctr}, Number, NTuple{N,Number}, Vector{<:Number}},
                           scale::Union{Type{<:Sca}, Number, NTuple{N,Number}, Vector{<:Number}},
                           dims, 
                           accumulator,
                           weight::Number) where{N, T} 
        offset = get_offset(size, offset)
        scale_init = get_scale(size, scale)
        scale = apply_dims(scale_init, dims, N) # replaces scale entries not in dims with zeros

        g(x) = weight .* ($(F[2])(x)) # 
        IndexFunArray(T,g, size) 
    end

    @eval function $(Symbol(:_, F[1]))(::Type{T}, size::NTuple{N, Int},
        offset,   # a version that supports an iterable collection of tuples for offset and scale 
        scale,
        dims=ntuple(+, N),
        accumulator = sum,
        weight=1) where{N, T} 

        # its important to use different variable names due to a julia bug
        offset_ = get_offset(size, offset)
        scale_init = get_scale(size, scale)
        scale_ = apply_dims(scale_init, dims, N)

        offset2 = cast_iter(offset_)
        # @show offset2
        scale2 = cast_iter(scale_)
        weight2 = cast_number_iter(weight)
        w = (x, offset, scale, weight3) -> weight3 .* $(F[2])(x) # adds offset and scale as requires parameters to the expression
        # the line below makes a function only depending on the index position but iterating over both, offsets and scales
        fkt = (idx) -> accumulator(curry(wrap_zipped(w),idx),zip(offset2,scale2,weight2))

        IndexFunArray(T, fkt, size) 
    end
    # default functions with offset and scaling behavior
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int};
                           offset=CtrFT,
                           scale=ScaUnit,
                           dims=ntuple(+, N), accumulator=sum, weight=1) where{N, T} 
        # $(F[1])(T, size, offset, scale, dims) 
        # print("Dispatch T3\n")
        $(Symbol(:_, F[1]))(T, size, offset, scale, dims, accumulator, weight) 
    end
    
    # change order of offset and scale
    @eval function $(F[1])(size::NTuple{N, Int}; 
                           offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum, weight=1) where {N}
        T = $DEFAULT_T
        $(F[1])(T, size, scale=scale, offset=offset, dims=dims, accumulator=accumulator, weight=weight) 
    end

    # convenient wrapper to provide an array as input
    @eval function $(F[1])(arr::AbstractArray{T, N}; 
                           offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum, weight=1) where {T, N}
        T2 = default_type(T, $DEFAULT_T)                   
        $(F[1])(T2, size(arr), scale=scale, offset=offset, dims=dims, accumulator=accumulator, weight=weight)
    end

    @eval export $(F[1])
end 

export axes1d

function axes1d(::Type{T}, size::NTuple{N, Int}; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), keepdims=false) where {T,N}
    offset = get_offset(size, offset)
    scale_init = get_scale(size, scale)
    scale = apply_dims(scale_init, dims, N) # replaces scale entries not in dims with zeros
    if keepdims
        (ramp(T, d, size[d]; offset=offset[d],scale=scale[d]) for d in dims)
    else
        (xx(T, (size[d],); offset=offset[d], scale=scale[d]) for d in dims)
    end
end

function axes1d(size::NTuple{N, Int}; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), keepdims=false) where {N}
    T = DEFAULT_T
    axes1d(T, size; offset=offset, scale=scale, dims=dims, keepdims=keepdims)
end

function axes1d(arr::AbstractArray{T, N}; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), keepdims=false) where {T, N}
    axes1d(size(arr); offset=offset, scale=scale, dims=dims, keepdims=keepdims)
end
