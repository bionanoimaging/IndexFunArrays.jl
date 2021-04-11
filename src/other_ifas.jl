export idx, cpx, exp_ikx
export gaussian, normal
export ramp

function ramp(::Type{T}, dim::Int, dim_size::Int;
    offset=CtrFT, scale=ScaUnit) where {T}
    size = single_dim_size(dim,dim_size)
    offset = get_offset(size, offset)
    scale_n = get_scale(size, scale)
    f = ((x) -> scale_n[dim] .* (x[dim] .- offset[dim]))
    IndexFunArray(T, f, size) 
end

function ramp(dim::Int, dim_size::Int; offset=CtrFT, scale=ScaUnit)
    ramp(DEFAULT_T, dim, dim_size; offset=offset, scale=scale)
end


# values in the complex plane
function cpx(::Type{T}, size::NTuple{N, Int}; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T}
    ig = idx(T, size, offset=offset, scale=scale, dims=dims).generator
    f(x) = complex(ig(x)...)
    return IndexFunArray(Complex{T}, f, size)
end
# values in the complex plane
function cpx(size::NTuple{N, Int}; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T}
    cpx(DEFAULT_T, size, offset=offset, scale=scale, dims=dims)
end
function cpx(arr::AbstractArray{T, N}; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T}
    cpx(T, size(arr), offset=offset, scale=scale, dims=dims)
end

function exp_ikx(::Type{T}, size::NTuple{N, Int}; shift_by=size.÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N), accumulator=sum) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(-2pi .* x .* y), get_scale(size, scale), shift_by)
    return exp_is(T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator)
end

function exp_ikx(size::NTuple{N, Int}; shift_by=size.÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N), accumulator=sum) where {N,T}
    myscale = apply_tuple_list((x,y)-> DEFAULT_T.(-2pi .* x .* y), get_scale(size, scale), shift_by)
    return exp_is(complex(DEFAULT_T), size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator)
end

function exp_ikx(arr::AbstractArray{T, N}; shift_by=size(arr).÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N), accumulator=sum) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(-2pi .* x .* y), get_scale(size, scale), shift_by)
    return exp_is(complex(typeof(arr[1])),size(arr), scale = myscale,  offset=offset, dims = dims, accumulator=accumulator)
end

function gaussian(::Type{T}, size::NTuple{N, Int}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(x ./( 2 .* y .*y)), get_scale(size, scale), sigma)
    return exp_sqr(T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator)
end

function gaussian(size::NTuple{N, Int}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum) where {N,T}
    myscale = apply_tuple_list((x,y)-> DEFAULT_T.(x ./( 2 .* y .*y)), get_scale(size, scale), sigma)
    return exp_sqr(DEFAULT_T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator)
end

function gaussian(arr::AbstractArray{T, N}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(x ./( 2 .*y.*y)), get_scale(size, scale), sigma)
    return exp_sqr(arr, scale = myscale,  offset=offset, dims = dims, accumulator=accumulator)
end

function normal(::Type{T}, size::NTuple{N, Int}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(x ./( 2 .*y.*y)), get_scale(size, scale), sigma)
    return exp_sqr_norm(T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator)
end

function normal(size::NTuple{N, Int}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum) where {N,T}
    myscale = apply_tuple_list((x,y)-> DEFAULT_T.(x ./( 2 .*y.*y)), get_scale(size, scale), sigma)
    return exp_sqr_norm(DEFAULT_T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator)
end

function normal(arr::AbstractArray{T, N}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(x ./( 2 .*y.*y)), get_scale(size, scale), sigma)
    return exp_sqr_norm(arr, scale = myscale,  offset=offset, dims = dims, accumulator=accumulator)
end
