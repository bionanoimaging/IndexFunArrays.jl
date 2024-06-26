export idx, cpx, exp_ikx, propagator
export gaussian, normal
export ramp, box, disc

function ramp(::Type{T}, dim::Int, dim_size::Int;
    offset=CtrFT, scale=ScaUnit) where {T}
    size = single_dim_size(dim,dim_size)
    offset = get_offset(size, offset)
    scale_n = get_scale(size, scale)
    f = ((x) -> T.(scale_n[dim] .* (x[dim] .- offset[dim])))
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
function cpx(size::NTuple{N, Int}; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N}
    cpx(DEFAULT_T, size, offset=offset, scale=scale, dims=dims)
end
function cpx(arr::AbstractArray{T, N}; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T}
    cpx(T, size(arr), offset=offset, scale=scale, dims=dims)
end

function exp_ikx(::Type{T}, size::NTuple{N, Int}; shift_by=size.÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N),
                accumulator=sum, weight=1) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(-2pi .* x .* y), get_scale(size, scale), optional_mat_to_iter(shift_by))
    return exp_is(T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator, weight=weight)
end

function exp_ikx(size::NTuple{N, Int}; shift_by=size.÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N), accumulator=sum, weight=1) where {N}
    myscale = apply_tuple_list((x,y)-> DEFAULT_T.(-2pi .* x .* y), get_scale(size, scale), optional_mat_to_iter(shift_by))
    return exp_is(complex(DEFAULT_T), size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator, weight=weight)
end

function exp_ikx(arr::AbstractArray{T, N}; shift_by=size(arr).÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N), accumulator=sum, weight=1) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(-2pi .* x .* y), get_scale(size(arr), scale), optional_mat_to_iter(shift_by))
    return exp_is(complex(eltype(arr)),size(arr), scale = myscale,  offset=offset, dims = dims, accumulator=accumulator)
end

function propagator(::Type{T}, size::NTuple{N, Int}; Δz=one(T), shift_by=0, k_max=T(0.5), offset=CtrFT, scale=ScaFT, dims=ntuple(+, N),
    accumulator=sum, weight=1) where {N,T}
myscale = apply_tuple_list((x,y)-> T.(x ./ y), get_scale(size, scale), optional_mat_to_iter(k_max))
myscale2 = apply_tuple_list((x,y)-> T.(-2pi .* x .* y), get_scale(size, scale), optional_mat_to_iter(shift_by))
return cispi.((2 .*Δz) .* 
    phase_kz(T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator, weight=weight) .+
    2 .* phase_kxy(T, size, scale = myscale2, offset=offset, dims = dims, accumulator=accumulator, weight=weight))
end

function propagator(size::NTuple{N, Int}; Δz=one(DEFAULT_T),  k_max=DEFAULT_T(0.5), shift_by=0, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N), accumulator=sum, weight=1) where {N}
myscale = apply_tuple_list((x,y)-> DEFAULT_T.(x ./ y), get_scale(size, scale), optional_mat_to_iter(k_max))
myscale2 = apply_tuple_list((x,y)-> DEFAULT_T.(-2pi .* x .* y), get_scale(size, scale), optional_mat_to_iter(shift_by))
return cispi.((2 .*Δz) .* 
        phase_kz(DEFAULT_T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator, weight=weight) .+
        2 .* phase_kxy(DEFAULT_T, size, scale = myscale2, offset=offset, dims = dims, accumulator=accumulator, weight=weight))
    end

function propagator(arr::AbstractArray{T, N}; Δz=one(T), k_max=T(0.5), shift_by=0, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N), accumulator=sum, weight=1) where {N,T}
myscale = apply_tuple_list((x,y)-> T.(x ./ y), get_scale(size(arr), scale), optional_mat_to_iter(k_max))
myscale2 = apply_tuple_list((x,y)-> T.(-2pi .* x .* y), get_scale(size(arr), scale), optional_mat_to_iter(shift_by))
return cispi.((2 .*Δz) .* 
        phase_kz(eltype(arr),size(arr), scale = myscale,  offset=offset, dims = dims, accumulator=accumulator) .+
        2 .* phase_kxy(eltype(arr), size, scale = myscale2, offset=offset, dims = dims, accumulator=accumulator, weight=weight))
end

# These versions of Gaussians interpret sigma as the standard-deviation along each axis separately
function gaussian(::Type{T}, size::NTuple{N, Int}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum, weight=1) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(x.^2 ./( 2 .* y .*y)), get_scale(size, scale), optional_mat_to_iter(sigma))
    return exp_sqr(T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator, weight=weight)
end

function gaussian(size::NTuple{N, Int}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum, weight=1) where {N}
    myscale = apply_tuple_list((x,y)-> DEFAULT_T.(x.^2 ./( 2 .* y .*y)), get_scale(size, scale), optional_mat_to_iter(sigma))
    return exp_sqr(DEFAULT_T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator, weight=weight)
end

function gaussian(arr::AbstractArray{T, N}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum, weight=1) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(x.^2 ./( 2 .*y.*y)), get_scale(size(arr), scale), optional_mat_to_iter(sigma))
    return exp_sqr(arr, scale = myscale,  offset=offset, dims = dims, accumulator=accumulator, weight=weight)
end

function normal(::Type{T}, size::NTuple{N, Int}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum, weight=1) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(x.^2 ./( 2 .*y.*y)), get_scale(size, scale), optional_mat_to_iter(sigma))
    return exp_sqr_norm(T, size, scale = myscale, offset=offset, dims = dims, accumulator=accumulator, weight=weight)
end

function normal(sz::NTuple{N, Int}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum, weight=1) where {N}
    myscale = apply_tuple_list((x,y)-> DEFAULT_T.(x.^2 ./( 2 .*y.*y)), get_scale(sz, scale), optional_mat_to_iter(sigma))
    return exp_sqr_norm(DEFAULT_T, sz, scale = myscale, offset=offset, dims = dims, accumulator=accumulator, weight=weight)
end

function normal(arr::AbstractArray{T, N}; sigma=1.0, offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N), accumulator=sum, weight=1) where {N,T}
    myscale = apply_tuple_list((x,y)-> T.(x.^2 ./( 2 .*y.*y)), get_scale(size(arr), scale), optional_mat_to_iter(sigma))
    return exp_sqr_norm(arr, scale = myscale,  offset=offset, dims = dims, accumulator=accumulator, weight=weight)
end

function box(::Type{T}, sz::NTuple{N, Int}, boxsize=sz./2; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T}
    myscale = get_scale(sz, scale)
    # with the quarterpixel offset we achieve the correct integer behaviour
    offset = get_offset(sz, offset) .- 0.25 
    return box1(T, sz, scale = myscale .* 2 ./ boxsize, offset=offset, dims = dims)
end

function box(sz::NTuple{N, Int}, boxsize=sz./2; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N}
    box(DEFAULT_T, sz, boxsize; offset=offset, scale=scale, dims=dims)
end

function box(::Type{TR}, arr::AbstractArray{T, N}, boxsize=size(arr)./2; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T,TR}
    return box(TR, size(arr), boxsize; offset=offset, scale=scale, dims=dims)
end

function box(arr::AbstractArray{T, N}, boxsize=size(arr)./2; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T}
    box(eltype(arr), size(arr), boxsize; offset=offset, scale=scale, dims=dims)
end

function disc(::Type{T}, sz::NTuple{N, Int}, disc_radius=sz./2; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T}
    myscale = get_scale(sz, scale)
    return disc1(T, sz, scale = myscale ./ disc_radius, offset=offset, dims = dims)
end

function disc(sz::NTuple{N, Int}, disc_radius=sz./2; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N}
    return disc(DEFAULT_T, sz, disc_radius; offset=offset, scale=scale, dims=dims)
end

function disc(::Type{TR}, arr::AbstractArray{T, N}, disc_radius=size(arr)./2; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T,TR}
    return disc(TR, size(arr), disc_radius; offset=offset, scale=scale, dims=dims)
end

function disc(arr::AbstractArray{T, N}, disc_radius=size(arr)./2; offset=CtrFT, scale=ScaUnit, dims=ntuple(+, N)) where {N,T}
    disc(eltype(arr), arr, disc_radius; offset=offset, scale=scale, dims=dims)
end
