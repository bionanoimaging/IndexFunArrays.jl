

"""
    rr2([T=Float64], size::size::NTuple{N, Int};
        offset=CtrFT,
        scale=ScaUnit)

Calculates the squared radius to a reference pixel.
In this case `CtrFT` is the center defined by the FFT convention.
`ScaUnit` leaves the values unscaled.
`offset` and `scale` can be either of `< :Ctr`, `Sca` respectively
or simply tuples with the same shape as `size`.
Look at `?Ctr` and `Sca` for all options.

Note that this function is based on a `GeneratorArray` and therefore does
not allocate the full memory needed to represent the array.

# Examples
```jldoctest
julia> rr2((4, 4))
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 8.0  5.0  4.0  5.0
 5.0  2.0  1.0  2.0
 4.0  1.0  0.0  1.0
 5.0  2.0  1.0  2.0

julia> rr2((3, 3))
3×3 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 2.0  1.0  2.0
 1.0  0.0  1.0
 2.0  1.0  2.0

julia> rr2((3, 3), offset=(1, 1))
3×3 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}:
 0.0  1.0  4.0
 1.0  2.0  5.0
 4.0  5.0  8.0

julia> rr2(Int, (3, 3), offset=(1, 1), scale=(10, 10))
3×3 GeneratorArray{Int64, 2, GeneratorArrays.var"#4#5"{Int64, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}:
   0  100  400
 100  200  500
 400  500  800
```
"""
rr2
