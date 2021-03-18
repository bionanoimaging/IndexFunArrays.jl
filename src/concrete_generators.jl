export rr, rr2
export xx, yy, zz
export phiphi

# define types to specify where the center point is
abstract type Ctr end  # Center of the array
struct CtrCorner <: Ctr end  # corner voxel is zero
struct CtrFFT <: Ctr end # corresponding to FFTs
struct CtrFT <: Ctr end # corresponding to FTs  (meaning shifted FFTs)
struct CtrMid <: Ctr end # middle of the array
struct CtrEnd <: Ctr end # other corner voxel is zero


get_offset(size, ::Type{CtrCorner}) = size.*0 .+ 1.0
get_offset(size, ::Type{CtrFT}) = size.÷2 .+ 1.0
get_offset(size, ::Type{CtrFFT}) = size.*0 .+ 1.0
get_offset(size, ::Type{CtrMid}) = (size.+1)./2.0
get_offset(size, ::Type{CtrEnd}) = size.+0.0
get_offset(size, t::NTuple) =  t


export Ctr,CtrCorner,CtrFFT,CtrFT,CtrMid,CtrEnd

abstract type Sca end # scaling of the array
struct ScaUnit <: Sca end # pixel distance is one
struct ScaNorm <: Sca end # total size along each dimension normalized to 1.0
struct ScaFT <: Sca end # reciprocal Fourier coordinates
struct ScaFTEdge <: Sca end # such that the edge of the Fourier space is 1.0

export Sca,ScaUnit,ScaNorm,ScaFT,ScaFTEdge

get_scale(size, ::Type{ScaUnit}) = ntuple(_ -> one(Int), length(size))
get_scale(size, ::Type{ScaNorm}) = 1 ./ (size .- 1)
get_scale(size, ::Type{ScaFT}) = 0.5 ./ (size .÷ 2)
get_scale(size, ::Type{ScaFTEdge}) = 1 ./ (size .÷ 2)  
get_scale(size, t::NTuple) = t 

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

function generate_window_functions_expr()
    x_exprW = :(min.(1.0, max.(0, - border_in .+ abs.(scale .* (x .- offset)))./(border_out-border_in)))

    functions = [
        (:(hanning),  :(x -> T(prod(cospi.(0.5 .* ($x_exprW)).^2)))),
    ]
    return functions
end

# we automatically generate the functions for rr2, rr, ...
# We set the types for the arguments correctly in the default cases
for F in generate_functions_expr() 
    # fallback type of GeneratorArray
    default_T = Float64

    # default functions with certain offset and scaling behavior
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int};
                           offset=CtrFT,
                           scale=ScaUnit) where{N, T} 
        offset = get_offset(size, offset)
        scale = get_scale(size, scale)
        GeneratorArray(T, $(F[2]), size) 
    end
    
    # change order of offset and scale
    @eval function $(F[1])(size::NTuple{N, Int}; 
                           offset=CtrFT, scale=ScaUnit) where {N}
        T = $default_T 
        $(F[1])(T, size, scale=scale, offset=offset) 
    end

    # convenient wrapper to provide an array as input
    @eval function $(F[1])(arr::AbstractArray{T, N}; 
                           offset=CtrFT, scale=ScaUnit) where {T, N}
        $(F[1])(T, size(arr), scale=scale, offset=offset)
    end

    @eval export $(F[1])
end 

# we automatically generate the functions for rr2, rr, ...
# We set the types for the arguments correctly in the default cases
for F in generate_window_functions_expr() 
    # fallback type of GeneratorArray
    default_T = Float64

    # default functions with certain offset and scaling behavior
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int};
                           offset=CtrFT,
                           scale=ScaFTEdge, border_in=0.8, border_out=1.0) where{N, T} 
        offset = get_offset(size, offset)
        scale = get_scale(size, scale)
        GeneratorArray(T, $(F[2]), size) 
    end
    
    # change order of offset and scale
    @eval function $(F[1])(size::NTuple{N, Int}; 
                           offset=CtrFT, 
                           scale=ScaFTEdge, border_in=0.8, border_out=1.0) where{N} 
        T = $default_T 
        $(F[1])(T, size, scale=scale, offset=offset, border_in=border_in, border_out=border_out) 
    end

    # convenient wrapper to provide an array as input
    @eval function $(F[1])(arr::AbstractArray{T, N}; 
                           offset=CtrFT, 
                           scale=ScaFTEdge, border_in=0.8, border_out=1.0) where{N, T} 
        $(F[1])(T, size(arr), scale=scale, offset=offset, border_in=border_in, border_out=border_out)
    end

    @eval export $(F[1])
end 
