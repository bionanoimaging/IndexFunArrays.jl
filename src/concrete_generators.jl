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

export Ctr,Ctr_Corner,Ctr_FFT,Ctr_FT,Ctr_Mid,Ctr_End

abstract type Sca end # scaling of the array
struct Sca_Unit <: Sca end # pixel distance is one
struct Sca_Norm <: Sca end # total size along each dimension normalized to 1.0
struct Sca_FT <: Sca end # reciprocal Fourier coordinates

export Sca,Sca_Unit,Sca_Norm,Sca_FT

get_scale(size, ::Type{Sca_Unit}) = size.*0 .+ 1.0
get_scale(size, ::Type{Sca_Norm}) = 1.0./size
get_scale(size, ::Type{Sca_FT}) = 1.0./size  # needs revision!

# List of functions and names
Fkts = [
    (:(rr2), :(x->sum(abs2.(x)))),
    (:(rr), :(x->sqrt.(sum(abs2.(x))))),
    (:(xx), :(x->x[1])),
    (:(yy), :(x->x[2])),
    (:(zz), :(x->x[3])),
    (:(phiphi), :(x->atan.(x[2],x[1]))),
]
for F in Fkts
    @eval $(F[1])(size::NTuple{N,Int}, offset::NTuple{N,Int}, scale::NTuple{N,Int}, ::Type{T}=Float64,) where{T,N,CT} = GeneratorArray(T, $(F[2]), convert(NTuple{N,T},offset), convert(NTuple{N,T},scale), size) 
    @eval $(F[1])(size::NTuple{N,Int},::Type{CT}=Ctr_FT,::Type{SC}=Sca_Unit, ::Type{T}=Float64,) where{T,N,CT} = GeneratorArray(T, $(F[2]), get_offset(size,CT), get_scale(size,SC), size) 
    @eval $(F[1])(anArray::AbstractArray{T,N}, offset::NTuple{N,Int},scale::NTuple{N,Int}, ::Type{T}=Float64,) where{T,N,CT} = GeneratorArray(T, $(F[2]), convert(NTuple{N,T},offset), convert(NTuple{N,T},scale), size(anArray)) 
    @eval $(F[1])(anArray::AbstractArray{T,N},::Type{CT}=Ctr_FT,::Type{SC}=Sca_Unit,::Type{T}=Float64,) where{T,N,CT} = GeneratorArray(T, $(F[2]), get_offset(size(anArray),CT), get_scale(size,SC), size(anArray)) 
    @eval export $(F[1])
end 

export getPropagator  # This has actually does evaluate. It should be programmed as a point-wise operation
function getPropagator(k0, dZ) 
    return x -> exp.(2im.*π.* sqrt.(max.(0.0,k0^2 .- rr2(x,Ctr_FT))))
end

#rr(size::NTuple{N,Int},::Type{T}=Float64,::Type{CT}=Ctr_Corner) where{T,N,CT} = GeneratorArray(x->sum(abs2.(x)), get_offset(size,CT), T, size) 
#rr(size::NTuple{N,Int}, ::Type{CT}=Ctr_Corner) where{N,CT} = rr(size, Float64, CT)
#rr(anArray::AbstractArray{T,N},::Type{T}) where{T,N,CT} = rr(size(anArray),T,Ctr_Corner) 
#rr(anArray::AbstractArray{T,N},::Type{CT}=Ctr_Corner) where{T,N,CT} = rr(anArray, Float64,CT)

