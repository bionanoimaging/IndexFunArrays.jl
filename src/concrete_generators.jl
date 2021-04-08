export Sca,ScaUnit,ScaNorm,ScaFT,ScaFTEdge
# export ScaRFT,ScaRFTEdge
export Ctr,CtrCorner,CtrFFT,CtrFT,CtrMid,CtrEnd, CtrRFFT,CtrRFT  # These do not work, as they need the size information about the real-space array
export rr, rr2
export xx, yy, zz, ee, tt, ramp
export phiphi
export idx, cpx, exp_ikx

using LinearAlgebra


const DEFAULT_T = Float64

# These are the type promotion rules, taken from float.jl but written in terms of types
# see also 
promote_type()
default_type(::Type{Bool}, def_T)    = def_T
default_type(::Type{Int8}, def_T)    = def_T
default_type(::Type{Int16}, def_T)   = def_T
default_type(::Type{Int32}, def_T)   = def_T
default_type(::Type{Int64}, def_T)   = def_T # LOSSY
default_type(::Type{Int128}, def_T)  = def_T # LOSSY
default_type(::Type{UInt8}, def_T)   = def_T
default_type(::Type{UInt16}, def_T)  = def_T
default_type(::Type{UInt32}, def_T)  = def_T
default_type(::Type{UInt64}, def_T)  = def_T # LOSSY
default_type(::Type{UInt128}, def_T) = def_T # LOSSY
default_type(::Type{T}, def_T) where{T} = T # all other types remain to be the same


# define types to specify where the center point is
abstract type Ctr end  # Center of the array
struct CtrCorner <: Ctr end  # corner voxel is zero
struct CtrFFT <: Ctr end # corresponding to FFTs
struct CtrRFFT <: Ctr end # corresponding to RFFTs
struct CtrFT <: Ctr end # corresponding to FTs  (meaning shifted FFTs)
struct CtrRFT <: Ctr end # corresponding to RFTs
struct CtrMid <: Ctr end # middle of the array
struct CtrEnd <: Ctr end # other corner voxel is zero

"""
    Ctr

Abstract type to specify the reference position
from which several other types subtype.

# Possible subtypes
* `CtrCorner`: Set the reference pixel in the corner
* `CtrFFT`: Set the reference pixel to the FFT center.
* `CtrFT`: Set the reference pixel to the FT center. FT means that the zero frequency is at the FFT convention center (`size ÷ 2 + 1`).
* `CtrRFFT`: Set the reference pixel to the RFFT center. Same as `CtrFFT` but the first dimension has center at 1. 
* `CtrRFT`: Set the reference pixel to the RFT center. FT means that the zero frequency is at the FFT convention center (`size ÷ 2 + 1`). 
            Same as `CtrFT` but the first dimension has center at 1.
* `CtrMid`: Set the reference pixel to real mid. For uneven arrays it is the center pixel, for even arrays it is the centered around a half pixel.
* `CtrEnd` Set the reference to the end corner (last pixel)
"""
Ctr

get_offset(size, ::Type{CtrCorner}) = size.*0 .+ 1.0
get_offset(size, ::Type{CtrFT}) = size.÷2 .+ 1.0
get_offset(size, ::Type{CtrRFT}) = Base.setindex(size.÷2,0,1) .+ 1.0
get_offset(size, ::Type{CtrFFT}) = size.*0 .+ 1.0
get_offset(size, ::Type{CtrRFFT}) = size.*0 .+ 1.0
get_offset(size, ::Type{CtrMid}) = (size.+1)./2.0
get_offset(size, ::Type{CtrEnd}) = size.+0.0
get_offset(size, t::Number) = ntuple(i -> t, length(size))
get_offset(size, t::NTuple) = t



abstract type Sca end # scaling of the array
struct ScaUnit <: Sca end # pixel distance is one
struct ScaNorm <: Sca end # total size along each dimension normalized to 1.0
struct ScaFT <: Sca end # reciprocal Fourier coordinates
struct ScaRFT <: Sca end # reciprocal Fourier coordinates for rFTs. 
struct ScaFTEdge <: Sca end # such that the edge of the Fourier space is 1.0
struct ScaRFTEdge <: Sca end # such that the edge of the Fourier space is 1.0

"""
    Sca 

Abstract type to indicate a scaling from which several other types subtype.

# Possible subtypes
* `ScaUnit`: No scaling of the indices 
* `ScaNorm`: Total length along each dimension is normalized to 1
* `ScaFT`: Reciprocal Fourier coordinates
* `ScaFTEdge`: Such that the edge (in FFT sense) of the pixel is 1.0
"""
Sca

get_scale(size, ::Type{ScaUnit}) = ntuple(_ -> one(Int), length(size))
get_scale(size, ::Type{ScaNorm}) = 1 ./ (max.(size .- 1, 1)) 
get_scale(size, ::Type{ScaFT}) = 0.5 ./ (max.(size ./ 2, 1))
# get_scale(size, ::Type{ScaRFT}) = 0.5 ./ (max.(Base.setindex(size./ 2,size[1]-1,1), 1))  # These scales are wrong! They need the information on the real-space size!
get_scale(size, ::Type{ScaFTEdge}) = 1 ./ (max.(size ./ 2, 1))  
# get_scale(size, ::Type{ScaRFTEdge}) = 1 ./ (max.(Base.setindex(size./ 2,size[1]-1,1), 1))
get_scale(size, t::Number) = ntuple(i -> t, length(size)) 
get_scale(size, t::NTuple) = t 

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
    ]
    return functions
end


function generate_window_functions_expr()
    # x_exprW = :(clamp.(- border_in .+ abs.(scale .* (x .- offset))./(border_out-border_in),0,1))
    x_exprW = :(clamp.(1 .-(abs.(scale .* (x .- offset)).-border_in)./(border_out .- border_in),0,1)) 
    x_exprRW = :(clamp.(1 .-(sqrt.(sum((scale .* (x .- offset)).^2)).-border_in)./(border_out .- border_in),0,1)) 

    functions = [  # see https://en.wikipedia.org/wiki/Window_function 
        (:(window_linear),  :(x -> T(prod(($x_exprW))))),
        (:(window_edge),  :(x -> T(prod(($x_exprW).>0.5)))),
        (:(window_hanning),  :(x -> T(prod(sinpi.(0.5 .* ($x_exprW)).^2)))),
        (:(window_half_cos),  :(x -> T(prod(sinpi.(0.5 .* ($x_exprW)))))),
        (:(window_hamming),  :(x -> T(prod(0.54 .-0.46.*cospi.(($x_exprW)))))),
        (:(window_blackman_harris),  :(x -> T(prod(0.35875 .- 0.48829.*cospi.($x_exprW).+0.14128.*cospi.(2 .*$x_exprW).-0.01168.*cospi.(3 .*$x_exprW))))),
        (:(window_radial_linear),  :(x -> T($x_exprRW))),
        (:(window_radial_edge),  :(x -> T(($x_exprRW).>0.5))),
        (:(window_radial_hanning),  :(x -> T(sinpi.(0.5 .* ($x_exprRW)).^2))),
        (:(window_radial_hamming),  :(x -> T((0.54 .-0.46.*cospi.(($x_exprRW)))))),
        (:(window_radial_blackman_harris),  :(x -> T((0.35875 .-0.48829.*cospi.($x_exprRW).+0.14128.*cospi.(2 .*$x_exprRW).-0.01168.*cospi.(3 .*$x_exprRW))))),
    ]
    return functions
end


function generate_tuple_functions_expr()
    x_expr = :(scale .* (x .- offset))

    functions = [
        (:(idx),  :(x -> T.($x_expr))), # returns a tuple
    ]
    return functions
end


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

function exp_ikx(::Type{T}, size::NTuple{N, Int}; shift_by=size.÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N)) where {N,T}
    return exp_is(T, size, scale = -2pi .* get_scale(size, scale) .* shift_by, offset=offset, dims = dims)
end

function exp_ikx(size::NTuple{N, Int}; shift_by=size.÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N)) where {N,T}
    return exp_is(complex(DEFAULT_T), size, scale = -2pi .* get_scale(size, scale) .* shift_by, offset=offset, dims = dims)
end

function exp_ikx(arr::AbstractArray{T, N}; shift_by=size(arr).÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N)) where {N,T}
    return exp_is(complex(typeof(arr[1])),size(arr), scale = -2pi .* get_scale(size(arr), scale) .* shift_by,  offset=offset, dims = dims)
end

# complex exponential for shifting
#=
function exp_ikx(::Type{T}, size::NTuple{N, Int}; shift_by=size.÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N)) where {N,T}
    scale_n = T.(get_scale(size, scale))
    offset = T.(get_offset(size, offset))
    finds(x) = (scale_n .* (x .- offset))
    f(x) = cis(dot((2π .* T.(.-shift_by)), finds(x)))
    return IndexFunArray(typeof(f(size)), f, size)
end

function exp_ikx(size::NTuple{N, Int}; shift_by=size.÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N)) where {N,T}
    exp_ikx(DEFAULT_T, size, shift_by=shift_by, offset=offset, scale=scale, dims=dims)
end
function exp_ikx(arr::AbstractArray{T, N}; shift_by=size(arr).÷2, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N)) where {N,T}
    exp_ikx(T, size(arr), shift_by=shift_by, offset=offset, scale=scale, dims=dims)
end
=#

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


# we automatically generate the functions for different windows like hanning 
for F in generate_window_functions_expr() 
    # default functions with certain offset and scaling behavior
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int};
                           offset=CtrFT,
                           scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N)) where{N, T} 
        scale_init = get_scale(size, scale)
        scale = ntuple(i -> i ∈ dims ? scale_init[i] : zero(scale_init[1]), N)
        offset = get_offset(size, offset)
        IndexFunArray(T, $(F[2]), size) 
    end
    
    # change order of offset and scale
    @eval function $(F[1])(size::NTuple{N, Int}; 
                           offset=CtrFT, 
                           scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N)) where{N} 
        T = $DEFAULT_T 
        $(F[1])(T, size, scale=scale, offset=offset, border_in=border_in, border_out=border_out, dims=dims) 
    end

    # convenient wrapper to provide an array as input
    @eval function $(F[1])(arr::AbstractArray{T, N}; 
                           offset=CtrFT, 
                           scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N)) where{N, T} 
        $(F[1])(default_type(T, $DEFAULT_T), size(arr), scale=scale, offset=offset, border_in=border_in, border_out=border_out, dims=dims)
    end

    @eval export $(F[1])
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



# include docstrings
include("concrete_generators_docstrings.jl")
