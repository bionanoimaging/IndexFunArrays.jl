module GeneratorArrays

export GeneratorArray

include("concrete_generators.jl")

 # T refers to the result type
struct GeneratorArray{T,N,F,H} <: AbstractArray{T, N} where {F}
    # stores the generator function to be applied to the indices. 
    generator::F     
    #a offset which is subtracted (.-) from the indices
    offset::NTuple{N,H}
    #a scaling which is applied after subtraction (.*) s
    scale::NTuple{N, H}
    # output size of the array 
    size::NTuple{N, Int}

    function GeneratorArray(::Type{T}, gen::F, 
                            ofs::NTuple{N, G}, sca::NTuple{N, H}, 
                            size::NTuple{N,Int}) where {T,N,F,G,H}
        if typeof(gen(size)) != T
            error("The generator function does not have Type T as indicated")
        end
        ofs = convert(NTuple{N, eltype(T)}, ofs)
        sca = convert(NTuple{N, eltype(T)}, sca)
        gen2(x) = gen(sca .* (x .- ofs))
        F_new = typeof(gen2)
        return new{T, N, F_new, eltype(T)}(gen2, ofs, sca, size) 
    end

end


# We define three easy constructors
# * generator function, size
# * generator, offset, scaling, size
# * generator, offset, size
# 
# And all of the above are available with optional type T

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

This is the simplest GeneratorArray constructor and only suited for arrays
containing `<:Numbers`.

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
```
"""
function GeneratorArray(gen::F, size::NTuple{N,Int}) where {N,F}
    # type T of the generator array is not provided
    # evaluate it once to get a type
    # but can be wrong if the generator function is not type stable
    T = typeof(gen(size))
    return GeneratorArray(T, gen, size)
end

function GeneratorArray(::Type{T}, gen::F, size::NTuple{N,Int}) where {T, N,F}
    if typeof(gen(size)) != T
        error("The generator function does not have Type T as indicated")
    end
    GeneratorArray(T, gen, 
                      ntuple(x -> zero(T), length(size)), 
                      ntuple(x -> one(T), length(size)),
                      size);
end

"""
    GeneratorArray(gen::F, 
                   ofs::NTuple{N, G},
                   size::NTuple{N,Int}) where {N,F,G,H}

`ofs` is a offset subtracted from the array indices.

# Examples
```julia-repl
julia> GeneratorArray(x -> sum(x), (1,1), (3, 3))
3×3 GeneratorArray{Int64, 2, var"#19#20"}:
 0  1  2
 1  2  3
 2  3  4
```
"""
function GeneratorArray(gen::F, 
                        ofs::NTuple{N, G},
                        size::NTuple{N,Int}) where {N,F,G,H}
    T = typeof(gen(size))
    return GeneratorArray(T, gen, ofs, size) 
end

function GeneratorArray(::Type{T}, gen::F, 
                        ofs::NTuple{N, G},
                        size::NTuple{N,Int}) where {T,N,F,G,H}
    if typeof(gen(size)) != T
        error("The generator function does not have Type T as indicated")
    end
    ofs = convert(NTuple{N, T}, ofs)
    sca = ntuple(x -> one(T), length(size))
    return GeneratorArray(gen, ofs, sca, size) 
end

"""
    GeneratorArray(gen::F, 
                   ofs::NTuple{N, G}, sca::NTuple{N, H}, 
                   size::NTuple{N,Int}) where {N,F,G,H}

`ofs` is a offset subtracted from the array indices.
`sca` is a scaling parameter subtracted from the array indices.

This functions is generic and is able to generate non-primitive 
element array types.

# Examples
```julia-repl
julia> GeneratorArray(x -> sum(x), (1,1), (5,5), (3, 3))
3×3 GeneratorArray{Int64, 2, var"#23#24"}:
  0   5  10
  5  10  15
 10  15  20

julia> GeneratorArray(x -> x, (1,0), (1,1), (3, 3))
3×3 GeneratorArray{Int64, 2, var"#39#40"}:
 (0, 1)  (0, 2)  (0, 3)
 (1, 1)  (1, 2)  (1, 3)
 (2, 1)  (2, 2)  (2, 3)
```
"""
function GeneratorArray(gen::F, 
                        ofs::NTuple{N, G}, sca::NTuple{N, H}, 
                        size::NTuple{N,Int}) where {N,F,G,H}
    T = typeof(gen(size))
    return GeneratorArray(T, gen, ofs, sca, size)
end



# define AbstractArray function to allow to treat the generator as an array
# See https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
Base.size(A::GeneratorArray) = A.size
Base.similar(A::GeneratorArray, ::Type{T}, size::Dims) where {T} = begin
    GeneratorArray(A.generator, A.offset, A.scale, size)
end

@inline ind_manip(I, scale, offset) = scale .* (I .- offset)

# subtract offset and apply scaling afterwards
Base.getindex(A::GeneratorArray{T,N}, I::Vararg{Int, N}) where {T,N} = begin
    # moving the index maninpulation into the generator
    # increases performance
    #return A.generator((I .- A.offset))
    return A.generator(I)
end

# not possible
Base.setindex!(A::GeneratorArray{T,N}, v, I::Vararg{Int,N}) where {T,N} = begin 
    error("Attempt to assign entries to GeneratorArray which is immutable.")
end

end # module
