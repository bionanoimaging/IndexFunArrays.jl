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
## Application to selected dimensions
Note that the code below yields a 3D array but with a one-sized trailing dimension. This can then be used for broadcasting.
```jldoctest
julia> x = ones(5,6,5);

julia> y=rr2(size(x,(1,2)))
5×6×1 GeneratorArray{Float64, 3, GeneratorArrays.var"#4#5"{Float64, Tuple{Float64, Float64, Float64}, Tuple{Int64, Int64, Int64}}}:
[:, :, 1] =
 13.0  8.0  5.0  4.0  5.0  8.0
 10.0  5.0  2.0  1.0  2.0  5.0
  9.0  4.0  1.0  0.0  1.0  4.0
 10.0  5.0  2.0  1.0  2.0  5.0
 13.0  8.0  5.0  4.0  5.0  8.0
```
Similarly you can also use dimensions 2 and 3 yielding an array of `size(y) == (1,6,5)`. 
Note that the necessary modification to the `Base.size` function is currently provided by this package.
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
    xx([T=Float64], size::NTuple{N, Int};
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
    yy([T=Float64], size::NTuple{N, Int};
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
    zz([T=Float64], size::NTuple{N, Int};
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

"""
    window_linear([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional (separable) window with a linear transition from zero at the borders (`border_out`) to one (`border_in`).
```jldoctest
julia> window_linear((8,9),border_in=0.0)
8×9 GeneratorArray{Float64, 2, GeneratorArrays.var"#34#35"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Float64, Float64}}:
 0.0  0.0     0.0    0.0     0.0   0.0     0.0    0.0     0.0
 0.0  0.0625  0.125  0.1875  0.25  0.1875  0.125  0.0625  0.0
 0.0  0.125   0.25   0.375   0.5   0.375   0.25   0.125   0.0
 0.0  0.1875  0.375  0.5625  0.75  0.5625  0.375  0.1875  0.0
 0.0  0.25    0.5    0.75    1.0   0.75    0.5    0.25    0.0
 0.0  0.1875  0.375  0.5625  0.75  0.5625  0.375  0.1875  0.0
 0.0  0.125   0.25   0.375   0.5   0.375   0.25   0.125   0.0
 0.0  0.0625  0.125  0.1875  0.25  0.1875  0.125  0.0625  0.0
```
"""
window_linear

"""
    window_radial_linear([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional radial window with a linear transition from zero at the borders (`border_out`) to one (`border_in`).
With the default offset and scale the borders are specified relative to the edge.
```jldoctest
julia> window_radial_linear((8,9),border_in=0.0)
8×9 GeneratorArray{Float64, 2, GeneratorArrays.var"#59#60"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Float64, Float64}}:
 0.0  0.0        0.0        0.0       0.0   0.0       0.0        0.0        0.0
 0.0  0.0        0.0986122  0.209431  0.25  0.209431  0.0986122  0.0        0.0
 0.0  0.0986122  0.292893   0.440983  0.5   0.440983  0.292893   0.0986122  0.0
 0.0  0.209431   0.440983   0.646447  0.75  0.646447  0.440983   0.209431   0.0
 0.0  0.25       0.5        0.75      1.0   0.75      0.5        0.25       0.0
 0.0  0.209431   0.440983   0.646447  0.75  0.646447  0.440983   0.209431   0.0
 0.0  0.0986122  0.292893   0.440983  0.5   0.440983  0.292893   0.0986122  0.0
 0.0  0.0        0.0986122  0.209431  0.25  0.209431  0.0986122  0.0        0.0
```
"""
window_radial_linear


"""
    window_edge([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional (separable) window with a sudden transition half way between the borders (`border_out`) to one (`border_in`).
See `?window_linear` for more details on the arguments.
"""
window_edge

"""
    window_radial_edge([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional radial window (disk) with a sudden transition half way between the borders (`border_out`) to one (`border_in`).
See `?window_radial_linear` for more details on the arguments.
"""
window_radial_edge

"""
    window_hanning([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional (separable) window with a von Hann transition between the borders (`border_out`) to one (`border_in`).
See `?window_linear` for more details on the arguments.
"""
window_hanning

"""
    window_radial_hanning([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional radial window with a von Hann transition between the borders (`border_out`) to one (`border_in`).
See `?window_radial_linear` for more details on the arguments.
"""
window_radial_hanning

"""
    window_hamming([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional (separable) window with a Hamming transition between the borders (`border_out`) to one (`border_in`).
See `?window_linear` for more details on the arguments.
"""
window_hamming

"""
    window_radial_hamming([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional radial window with a Hamming transition between the borders (`border_out`) to one (`border_in`).
See `?window_radial_linear` for more details on the arguments.
"""
window_radial_hamming

"""
    window_blackman_harris([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional (separable) window with a  transition according to Blackman/Harris between the borders (`border_out`) to one (`border_in`).
See `?window_linear` for more details on the arguments.
"""
window_blackman_harris

"""
    window_radial_blackman_harris([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0)  

A multidimensional radial window with a Hamming transition according to Blackman/Harris between the borders (`border_out`) to one (`border_in`).
See `?window_radial_linear` for more details on the arguments.
"""
window_radial_blackman_harris

