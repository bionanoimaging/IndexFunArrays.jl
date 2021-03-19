module GeneratorArrays

export GeneratorArray
import Base.size
export size

include("concrete_generators.jl")

 # T refers to the result type
struct GeneratorArray{T, N, F} <: AbstractArray{T, N} where {F}
    # stores the generator function to be applied to the indices. 
    generator::F     
    #a offset which is subtracted (.-) from the indices
    # output size of the array 
    size::NTuple{N, Int}

    # Constructor function
    function GeneratorArray(::Type{T}, gen::F, size::NTuple{N,Int}) where {T,N,F,G,H}
        if typeof(gen(size)) != T
            error("The generator function does not have Type T as indicated")
        end
        return new{T, N, F}(gen, size) 
    end

end


"""
    GeneratorArray([T], gen::F, size::NTuple{N,Int}) where {N,F}

Generate a GeneratorArray object which behaves like an array
but does not allocate the full array. Instead it calculates
the elements when needed. This is useful to prevent array allocations.
`gen` is a function which takes the array indices wrapped as tuple as input.
The output of `gen` determines the element type of the resulting array.
`size` is the output size of the resulting array.
`T` can be the optional element type of the arrays. 
`gen` needs to have `T` as return type, otherwise the GeneratorArray
might be type unstable.

# Examples
```julia-repl
julia> GeneratorArray(x -> sum(x), (3, 3))
3×3 GeneratorArray{Int64, 2, var"#182#183"}:
 2  3  4
 3  4  5
 4  5  6

julia> GeneratorArray(x -> sum(abs2.(x)), (3, 3))
3×3 GeneratorArray{Int64, 2, var"#184#185"}:
  2   5  10
  5   8  13
 10  13  18

julia> GeneratorArray(x -> (x[1], x[2], "Julia"), (3,3))
3×3 GeneratorArray{Tuple{Int64, Int64, String}, 2, var"#18#19"}:
 (1, 1, "Julia")  (1, 2, "Julia")  (1, 3, "Julia")
 (2, 1, "Julia")  (2, 2, "Julia")  (2, 3, "Julia")
 (3, 1, "Julia")  (3, 2, "Julia")  (3, 3, "Julia")
```
"""
function GeneratorArray(gen::F, size::NTuple{N,Int}) where {N,F}
    # type T of the generator array is not provided
    # evaluate it once to get a type
    # but can be wrong if the generator function is not type stable
    T = typeof(gen(size))
    return GeneratorArray(T, gen, size)
end

# define AbstractArray function to allow to treat the generator as an array
# See https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
Base.size(A::GeneratorArray) = A.size
Base.similar(A::GeneratorArray, ::Type{T}, size::Dims) where {T} = GeneratorArray(A.generator, size)

# calculate the entry according to the index
Base.getindex(A::GeneratorArray{T,N}, I::Vararg{Int, N}) where {T,N} = return A.generator(I)

# not supported
Base.setindex!(A::GeneratorArray{T,N}, v, I::Vararg{Int,N}) where {T,N} = begin 
    error("Attempt to assign entries to GeneratorArray which is immutable.")
end



"""
    size(x::AbstractArray,dim::NTuple; keep_dims=true)

Overloads `Base.size` and allows to access the size at several dimensions
in one call.

# Examples
```jldoctest
julia> x = ones((2,4,6,8, 10));

julia> size(x, (2,3))
(1, 4, 6, 1, 1)

julia> size(x, (2,3,4), keep_dims=false)
(4, 6, 8)

julia> size(x, (4,3,2), keep_dims=false)
(8, 6, 4)
```

"""
function size(x::AbstractArray{T},dim::NTuple{N,Int}; keep_dims=true) where{T,N}
    if ~keep_dims
        return map(n->size(x,n),dim)
    end
    sz=ones(Int, ndims(x))
    for n in dim
        sz[n]=size(x,n) 
    end
    return Tuple(sz)
end 

end # module
