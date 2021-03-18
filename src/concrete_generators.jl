export rr, rr2
export xx, yy, zz
export phiphi

# define types to specify where the center point is
abstract type Ctr end  # Center of the array
struct Ctr_Corner <: Ctr end  # corner voxel is zero
struct Ctr_FFT <: Ctr end # corresponding to FFTs
struct Ctr_FT <: Ctr end # corresponding to FTs  (meaning shifted FFTs)
struct Ctr_Mid <: Ctr end # middle of the array
struct Ctr_End <: Ctr end # other corner voxel is zero


get_offset(size, ::Type{Ctr_Corner}) = size.*0 .+ 1.0
get_offset(size, ::Type{Ctr_FT}) = size.÷2 .+ 1.0
get_offset(size, ::Type{Ctr_FFT}) = size.*0 .+ 1.0
get_offset(size, ::Type{Ctr_Mid}) = (size.+1)./2.0
get_offset(size, ::Type{Ctr_End}) = size.+0.0
get_offset(size, t::NTuple) =  t


export Ctr,Ctr_Corner,Ctr_FFT,Ctr_FT,Ctr_Mid,Ctr_End

abstract type Sca end # scaling of the array
struct Sca_Unit <: Sca end # pixel distance is one
struct Sca_Norm <: Sca end # total size along each dimension normalized to 1.0
struct Sca_FT <: Sca end # reciprocal Fourier coordinates

export Sca,Sca_Unit,Sca_Norm,Sca_FT

get_scale(size, ::Type{Sca_Unit}) = size .* 0 .+ 1.0
get_scale(size, ::Type{Sca_Norm}) = 1.0 ./ size
get_scale(size, ::Type{Sca_FT}) = 1.0 ./ size  # needs revision!
get_scale(size, t::NTuple) = t  # needs revision!

# List of functions and names we want to predefine
function generate_functions_expr()
    # offset and scale is already wrapped in the generator function
    x_expr = :(scale .* (x .- offset))
    x_expr1 = :(scale[1] .* (x[1] .- offset[1]))
    x_expr2 = :(scale[2] .* (x[2] .- offset[2]))
    x_expr3 = :(scale[3] .* (x[3] .- offset[3]))

    functions = [
        (:(rr2), :(x -> T(sum(abs2.($x_expr))))),
        (:(rr),  :(x -> T(sqrt.(sum(abs2.($x_expr)))))),
        (:(xx),  :(x -> T($x_expr1))),
        (:(yy),  :(x -> T($x_expr2))),
        (:(zz),  :(x -> T($x_expr3))),
        (:(phiphi), :(x -> T(atan.($x_expr2, $x_expr3)))),
    ]
    return functions
end

# we automatically generate the functions for rr2, rr, ...
# We set the types for the arguments correctly in the default cases
for F in generate_functions_expr() 
    # fallback type of GeneratorArray
    default_T = Float64

    # default functions with certain offset and scaling behavior
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int},
                           offset::Type{CT}=Ctr_FT,
                           scale::Type{SC}=Sca_Unit) where{SC, CT, N, T} 
        offset = get_offset(size, CT)
        scale = get_scale(size, SC)
        GeneratorArray(T, $(F[2]), size) 
    end

    @eval function $(F[1])(size::NTuple{N, Int},
                           offset::Type{CT}=Ctr_FT,
                           scale::Type{SC}=Sca_Unit) where{SC, CT, N} 
        T = $default_T 
        $(F[1])(T, size, offset, scale) 
    end

    # tuple based offset and scale
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int},
                           offset::NTuple,
                           scale::NTuple) where{SC, CT, N, T} 
        GeneratorArray(T, $(F[2]), size) 
    end

    @eval function $(F[1])(size::NTuple{N, Int},
                           offset::NTuple,
                           scale::NTuple) where{SC, CT, N} 
        T = $default_T 
        $(F[1])(T, size, offset, scale) 
    end

    # only offset provided 
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int},
                           offset::NTuple) where{SC, CT, N, T} 
        scale = ntuple(_ -> one(T), length(size))
        GeneratorArray(T, $(F[2]), size) 
    end

    @eval function $(F[1])(size::NTuple{N, Int},
                           offset::NTuple) where{SC, CT, N} 
        T = $default_T 
        $(F[1])(T, size, offset) 
    end

    @eval export $(F[1])
end 

export getPropagator  # This has actually does evaluate. It should be programmed as a point-wise operation
function getPropagator(k0, dZ) 
    return x -> exp.(2im.*π.* sqrt.(max.(0.0,k0^2 .- rr2(x,Ctr_FT))))
end

#rr(size::NTuple{N,Int},::Type{T}=Float64,::Type{CT}=Ctr_Corner) where{T,N,CT} = GeneratorArray(x->sum(abs2.(x)), get_offset(size,CT), T, size) 
#rr(size::NTuple{N,Int}, ::Type{CT}=Ctr_Corner) where{N,CT} = rr(size, Float64, CT)
#rr(arr::AbstractArray{T,N},::Type{T}) where{T,N,CT} = rr(size(arr),T,Ctr_Corner) 
#rr(arr::AbstractArray{T,N},::Type{CT}=Ctr_Corner) where{T,N,CT} = rr(arr, Float64,CT)

