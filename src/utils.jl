"""
    single_dim_size(dim::Int,dim_size::Int)

Returns a tuple (length `dim`) of singleton sizes except at the final position `dim`, which contains `dim_size`

# Example
```jldoctest
julia> IndexFunArrays.single_dim_size(4, 3)
(1, 1, 1, 3)

julia> IndexFunArrays.single_dim_size(4, 5)
(1, 1, 1, 5)

julia> IndexFunArrays.single_dim_size(2, 5)
(1, 5)
```
"""
function single_dim_size(dim::Int,dim_size::Int)
    Base.setindex(Tuple(ones(Int, dim)),dim_size,dim)
end


"""
    selectsizes(x::AbstractArray, dim; keep_dims=true)

Additional size method to access the size at several dimensions
in one call.
`keep_dims` allows to return the other dimensions as singletons.

# Examples
```jldoctest
julia> x = ones((2,4,6,8, 10));

julia> selectsizes(x, (2,3))
(1, 4, 6, 1, 1)

julia> selectsizes(x, 5)
(1, 1, 1, 1, 10)

julia> selectsizes(x, (5,))
(1, 1, 1, 1, 10)

julia> selectsizes(x, (2,3,4), keep_dims=false)
(4, 6, 8)

julia> selectsizes(x, (4,3,2), keep_dims=false)
(8, 6, 4)
```

"""
function selectsizes(x::AbstractArray{T},dim::NTuple{N,Int};
                    keep_dims=true) where{T,N}
    if ~keep_dims
        return map(n->size(x,n),dim)
    end
    sz = ones(Int, ndims(x))
    for n in dim
        sz[n] = size(x,n) 
    end
    return Tuple(sz)
end 

function selectsizes(x::AbstractArray, dim::Integer; keep_dims=true)
    selectsizes(x, Tuple(dim), keep_dims=keep_dims)
end
