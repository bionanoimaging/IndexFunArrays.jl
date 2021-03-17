module GeneratorArrays

struct GeneratorArray{T,N,F} <: AbstractArray{T,N} where {F} # T refers to the result type
    generator::F      # stores the generator function to be applied to the indices. This NEEDs to be a defined type to avoid allocations
    offset::NTuple{N,T}
    scale::NTuple{N,T}
    dims::NTuple{N,Int}   # stores the sizes
end

GeneratorArray(gen::F, ofs::NTuple{N,T}, sca::NTuple{N,T}, ::Type{T}, dims::Int...) where {T,N,F} = GeneratorArray{T,N,F}(gen, ofs, sca, dims);
GeneratorArray(gen::F, ofs::NTuple{N,T}, sca::NTuple{N,T}, ::Type{T}, dims::NTuple{N,Int}) where {T,N,F} = GeneratorArray{T,N,F}(gen, ofs, sca, dims);
GeneratorArray(gen::F, ofs::NTuple{N,T}, ::Type{T}, dims::NTuple{N,Int}) where {T,N,F} = GeneratorArray{T,N,F}(gen, ofs, dims.*0 .+1.0, dims);
GeneratorArray(gen::F, ::Type{T}, dims::NTuple{N,Int}) where {T,N,F} = GeneratorArray{T,N,F}(gen, dims.*0 .+1.0, dims.*0 .+1.0, dims);
GeneratorArray(gen::F, dims::NTuple{N,Int}) where {T,N,F} = GeneratorArray{Float64,N,F}(gen, dims.*0 .+1.0, dims.*0 .+1.0, dims);

Base.size(A::GeneratorArray) = A.dims
Base.similar(A::GeneratorArray, ::Type{T}, dims::Dims) where {T} = GeneratorArray(A.generator, A.offset, A.scale, dims)
Base.getindex(A::GeneratorArray{T,N}, I::Vararg{Int,N}) where {T,N} = A.generator((I.-A.offset).*A.scale) # A.generator(convert(NTuple{N,T},I)) 
Base.setindex!(A::GeneratorArray{T,N}, v, I::Vararg{Int,N}) where {T,N} = (print("Warning: Attempt to assign to GeneratorArray\n"))

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

get_offset(dims, ::Type{Ctr_Corner}) = dims.*0 .+ 1.0
get_offset(dims, ::Type{Ctr_FT}) = dims.÷2 .+ 1.0
get_offset(dims, ::Type{Ctr_FFT}) = dims.*0 .+ 1.0
get_offset(dims, ::Type{Ctr_Mid}) = (dims.+1)./2.0
get_offset(dims, ::Type{Ctr_End}) = dims.+0.0

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
    @eval $(F[1])(dims::NTuple{N,Int}, offset::NTuple{N,Int}, ::Type{T}=Float64,) where{T,N,CT} = GeneratorArray($(F[2]), convert(NTuple{N,T},offset), T, dims) 
    @eval $(F[1])(dims::NTuple{N,Int},::Type{CT}=Ctr_FT, ::Type{T}=Float64,) where{T,N,CT} = GeneratorArray($(F[2]), get_offset(dims,CT), T, dims) 
    @eval $(F[1])(anArray::AbstractArray{T,N}, offset::NTuple{N,Int},::Type{T}=Float64,) where{T,N,CT} = GeneratorArray($(F[2]), convert(NTuple{N,T},offset), T, size(anArray)) 
    @eval $(F[1])(anArray::AbstractArray{T,N},::Type{CT}=Ctr_FT,::Type{T}=Float64,) where{T,N,CT} = GeneratorArray($(F[2]), get_offset(size(anArray),CT), T, size(anArray)) 
    @eval export $(F[1])
end 

export getPropagator  # This has actually does evaluate. It should be programmed as a point-wise operation
function getPropagator(k0, dZ) 
    return x -> exp.(2im.*π.* sqrt.(max.(0.0,k0^2 .- rr2(x,Ctr_FT))))
end

#rr(dims::NTuple{N,Int},::Type{T}=Float64,::Type{CT}=Ctr_Corner) where{T,N,CT} = GeneratorArray(x->sum(abs2.(x)), get_offset(dims,CT), T, dims) 
#rr(dims::NTuple{N,Int}, ::Type{CT}=Ctr_Corner) where{N,CT} = rr(dims, Float64, CT)
#rr(anArray::AbstractArray{T,N},::Type{T}) where{T,N,CT} = rr(size(anArray),T,Ctr_Corner) 
#rr(anArray::AbstractArray{T,N},::Type{CT}=Ctr_Corner) where{T,N,CT} = rr(anArray, Float64,CT)

end # module
