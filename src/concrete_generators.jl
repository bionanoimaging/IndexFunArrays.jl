abstract type Ctr end  # Center of the array
struct Ctr_Corner <: Ctr end  # corner voxel is zero
struct Ctr_FFT <: Ctr end # corresponding to FFTs
struct Ctr_FT <: Ctr end # corresponding to FTs  (meaning shifted FFTs)
struct Ctr_Mid <: Ctr end # middle of the array
struct Ctr_End <: Ctr end # other corner voxel is zero

export Ctr,Ctr_Corner,Ctr_FFT,Ctr_FT,Ctr_Mid,Ctr_End

abstract type Sca end # scaling of the array
struct Sca_Unit <: Sca end # pixel distance is one
struct Sca_Norm <: Sca end # total size along each dimension normalized to 1.0
struct Sca_FT <: Sca end # reciprocal Fourier coordinates

export Sca,Sca_Unit,Sca_Norm,Sca_FT

get_offset(size, ::Type{Ctr_Corner}) = size.*0 .+ 1.0
get_offset(size, ::Type{Ctr_FT}) = size.÷2 .+ 1.0
get_offset(size, ::Type{Ctr_FFT}) = size.*0 .+ 1.0
get_offset(size, ::Type{Ctr_Mid}) = (size.+1)./2.0
get_offset(size, ::Type{Ctr_End}) = size.+0.0

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
    @eval $(F[1])(size::NTuple{N,Int}, offset::NTuple{N,Int}, ::Type{T}=Float64,) where{T,N,CT} = GeneratorArray($(F[2]), convert(NTuple{N,T},offset), T, size) 
    @eval $(F[1])(size::NTuple{N,Int},::Type{CT}=Ctr_FT, ::Type{T}=Float64,) where{T,N,CT} = GeneratorArray($(F[2]), get_offset(size,CT), T, size) 
    @eval $(F[1])(anArray::AbstractArray{T,N}, offset::NTuple{N,Int},::Type{T}=Float64,) where{T,N,CT} = GeneratorArray($(F[2]), convert(NTuple{N,T},offset), T, size(anArray)) 
    @eval $(F[1])(anArray::AbstractArray{T,N},::Type{CT}=Ctr_FT,::Type{T}=Float64,) where{T,N,CT} = GeneratorArray($(F[2]), get_offset(size(anArray),CT), T, size(anArray)) 
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

