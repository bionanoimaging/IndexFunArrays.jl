"""
    rr2([T=Float64], size::size::NTuple{N, Int};
        offset=CtrFT,
        scale=ScaUnit)

Calculates the squared radius to a reference pixel.
In this case `CtrFT` is the center defined by the FFT convention.
`ScaUnit` leaves the values unscaled.
`offset` and `scale` can be either of `<:Ctr`, `<:Sca` respectively
or simply tuples with the same shape as `size`.
Look at `?Ctr` and `?Sca` for all options.

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
```

## Change Reference Position
```jldoctest
julia> rr2((3,3), offset=CtrCorner)
3×3 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 0.0  1.0  4.0
 1.0  2.0  5.0
 4.0  5.0  8.0

julia> rr2((4,4), offset=CtrFT)
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 8.0  5.0  4.0  5.0
 5.0  2.0  1.0  2.0
 4.0  1.0  0.0  1.0
 5.0  2.0  1.0  2.0

julia> rr2((4,4), offset=CtrMid)
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 4.5  2.5  2.5  4.5
 2.5  0.5  0.5  2.5
 2.5  0.5  0.5  2.5
 4.5  2.5  2.5  4.5

julia> rr2((4,4), offset=CtrEnd)
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 18.0  13.0  10.0  9.0
 13.0   8.0   5.0  4.0
 10.0   5.0   2.0  1.0
  9.0   4.0   1.0  0.0

julia> rr2((3, 3), offset=(1, 1))
3×3 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}:
 0.0  1.0  4.0
 1.0  2.0  5.0
 4.0  5.0  8.0
```

## Change Scaling
```jldoctest
julia> rr((4,4), scale=ScaUnit)
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 2.82843  2.23607  2.0  2.23607
 2.23607  1.41421  1.0  1.41421
 2.0      1.0      0.0  1.0
 2.23607  1.41421  1.0  1.41421

julia> rr((4,4), scale=ScaNorm)
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 0.942809  0.745356  0.666667  0.745356
 0.745356  0.471405  0.333333  0.471405
 0.666667  0.333333  0.0       0.333333
 0.745356  0.471405  0.333333  0.471405

julia> rr((4,4), scale=ScaFT)
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 0.707107  0.559017  0.5   0.559017
 0.559017  0.353553  0.25  0.353553
 0.5       0.25      0.0   0.25
 0.559017  0.353553  0.25  0.353553

julia> rr((4,4), scale=ScaFTEdge)
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 1.41421  1.11803   1.0  1.11803
 1.11803  0.707107  0.5  0.707107
 1.0      0.5       0.0  0.5
 1.11803  0.707107  0.5  0.707107

julia> rr2(Int, (3, 3), offset=(1, 1), scale=(10, 10))
3×3 GeneratorArray{Int64, 2, GeneratorArrays.var"#4#5"{Int64, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}:
   0  100  400
 100  200  500
 400  500  800
```
"""
rr2

"""
    rr([T=Float64], size::size::NTuple{N, Int};
       offset=CtrFT,
       scale=ScaUnit)

See `rr2` for all options.

# Examples
```jldoctest
julia> rr((3, 3))
3×3 GeneratorArray{Float64, 2, GeneratorArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 1.41421  1.0  1.41421
 1.0      0.0  1.0
 1.41421  1.0  1.41421

julia> rr((3, 3), offset=CtrCorner)
3×3 GeneratorArray{Float64, 2, GeneratorArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 0.0  1.0      2.0
 1.0  1.41421  2.23607
 2.0  2.23607  2.82843
```
"""
rr


"""
    xx([T=Float64], size::size::NTuple{N, Int};
       offset=CtrFT,
       scale=ScaUnit)

A distance ramp along first dimension.
```jldoctest
julia> xx((4,4))
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#14#15"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 -2.0  -2.0  -2.0  -2.0
 -1.0  -1.0  -1.0  -1.0
  0.0   0.0   0.0   0.0
  1.0   1.0   1.0   1.0
```
"""
xx


"""
    yy([T=Float64], size::size::NTuple{N, Int};
       offset=CtrFT,
       scale=ScaUnit)

A distance ramp along second dimension.
```jldoctest
julia> yy((4,4))
4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#19#20"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 -2.0  -1.0  0.0  1.0
 -2.0  -1.0  0.0  1.0
 -2.0  -1.0  0.0  1.0
 -2.0  -1.0  0.0  1.0
```
"""
yy



"""
    zz([T=Float64], size::size::NTuple{N, Int};
       offset=CtrFT,
       scale=ScaUnit)

A distance ramp along third dimension.
```jldoctest
julia> zz((1, 1, 4))
1×1×4 GeneratorArray{Float64, 3, GeneratorArrays.var"#24#25"{Float64, Tuple{Float64, Float64, Float64}, Tuple{Int64, Int64, Int64}}}:
[:, :, 1] =
 -2.0

[:, :, 2] =
 -1.0

[:, :, 3] =
 0.0

[:, :, 4] =
 1.0
```
"""
zz
