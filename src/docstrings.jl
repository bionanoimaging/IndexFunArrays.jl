"""
    idx([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit)

See `rr2` for a description of all options.

Returns basically the `CartesianIndices` but as a tuple and accounting for optional
`scale`, `offset` and data type.
Note that `T` is enforced element-wise for the return tuple elements.


```jldoctest
julia> idx(Int, (3,3), offset=CtrCorner)
3×3 IndexFunArray{Tuple{Int64, Int64}, 2, IndexFunArrays.var"#39#40"{Int64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 (0, 0)  (0, 1)  (0, 2)
 (1, 0)  (1, 1)  (1, 2)
 (2, 0)  (2, 1)  (2, 2)

julia> idx(Int, (3,3), offset=(0,0))
3×3 IndexFunArray{Tuple{Int64, Int64}, 2, IndexFunArrays.var"#39#40"{Int64, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}:
 (1, 1)  (1, 2)  (1, 3)
 (2, 1)  (2, 2)  (2, 3)
 (3, 1)  (3, 2)  (3, 3)

julia> idx(Float32, (3,3), offset=(0,0))
3×3 IndexFunArray{Tuple{Float32, Float32}, 2, IndexFunArrays.var"#39#40"{Float32, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}:
 (1.0, 1.0)  (1.0, 2.0)  (1.0, 3.0)
 (2.0, 1.0)  (2.0, 2.0)  (2.0, 3.0)
 (3.0, 1.0)  (3.0, 2.0)  (3.0, 3.0)
```
"""
idx

"""
    cpx([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit)

See `rr2` for a description of all options.

Returns an IndexFunArray where each positon corresponds to a complex value
according to its position. The parameters `offset` and `scale` can be used
accordingly (see `rr2`).
Note that `T` is enforced element-wise for the return tuple elements.


```jldoctest
julia> cpx(Int, (3,3), offset=CtrCorner)
3×3 Matrix{Complex{Int64}}:
 0+0im  0+1im  0+2im
 1+0im  1+1im  1+2im
 2+0im  2+1im  2+2im

julia> cpx(Int, (3,3), offset=(0,0))
3×3 Matrix{Complex{Int64}}:
 1+1im  1+2im  1+3im
 2+1im  2+2im  2+3im
 3+1im  3+2im  3+3im

julia> cpx((3,3), offset=(0,0))
3×3 Matrix{ComplexF64}:
 1.0+1.0im  1.0+2.0im  1.0+3.0im
 2.0+1.0im  2.0+2.0im  2.0+3.0im
 3.0+1.0im  3.0+2.0im  3.0+3.0im
```
"""
cpx


"""
    rr2([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

Calculates the squared radius to a reference pixel.
In this case `CtrFT` is the center defined by the FFT convention.
`ScaUnit` leaves the values unscaled.
`offset` and `scale` can be either of `<:Ctr`, `<:Sca` respectively
or simply tuples with the same shape as `size`.
Look at `?Ctr` and `?Sca` for all options.
`dims` is a keyword argument to specifiy over which dimensions the
operation will effectively happen.

The arguments `offset`, and `scale` support list-mode, which means that
supplying a tuple of tuples or a vector of tuples or a matrix causes the
function to automatically generate a superposition of multiple versions of itself.
The type of superposition is controlled by the `accumulator` argument. The relative strength
of the individual superposition is controlled via the `weight` argument, which can be
a tuple or vector. Have a look at the Voronoi-example below.

Note that this function is based on a `IndexFunArray` and therefore does
not allocate the full memory needed to represent the array.


# Examples
```jldoctest
julia> rr2((4, 4))
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 8.0  5.0  4.0  5.0
 5.0  2.0  1.0  2.0
 4.0  1.0  0.0  1.0
 5.0  2.0  1.0  2.0
```

## Change Reference Position
```jldoctest
julia> rr2((3,3), offset=CtrCorner)
3×3 IndexFunArray{Float64, 2, IndexFunArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 0.0  1.0  4.0
 1.0  2.0  5.0
 4.0  5.0  8.0

julia> rr2((4,4), offset=CtrFT)
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 8.0  5.0  4.0  5.0
 5.0  2.0  1.0  2.0
 4.0  1.0  0.0  1.0
 5.0  2.0  1.0  2.0

julia> rr2((4,4), offset=CtrMid)
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 4.5  2.5  2.5  4.5
 2.5  0.5  0.5  2.5
 2.5  0.5  0.5  2.5
 4.5  2.5  2.5  4.5

julia> rr2((4,4), offset=CtrEnd)
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 18.0  13.0  10.0  9.0
 13.0   8.0   5.0  4.0
 10.0   5.0   2.0  1.0
  9.0   4.0   1.0  0.0

julia> rr2((3, 3), offset=(1, 1))
3×3 IndexFunArray{Float64, 2, IndexFunArrays.var"#4#5"{Float64, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}:
 0.0  1.0  4.0
 1.0  2.0  5.0
 4.0  5.0  8.0
```

## Change Scaling
```jldoctest
julia> rr((4,4), scale=ScaUnit)
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 2.82843  2.23607  2.0  2.23607
 2.23607  1.41421  1.0  1.41421
 2.0      1.0      0.0  1.0
 2.23607  1.41421  1.0  1.41421

julia> rr((4,4), scale=ScaNorm)
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 0.942809  0.745356  0.666667  0.745356
 0.745356  0.471405  0.333333  0.471405
 0.666667  0.333333  0.0       0.333333
 0.745356  0.471405  0.333333  0.471405

julia> rr((4,4), scale=ScaFT)
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 0.707107  0.559017  0.5   0.559017
 0.559017  0.353553  0.25  0.353553
 0.5       0.25      0.0   0.25
 0.559017  0.353553  0.25  0.353553

julia> rr((4,4), scale=ScaFTEdge)
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 1.41421  1.11803   1.0  1.11803
 1.11803  0.707107  0.5  0.707107
 1.0      0.5       0.0  0.5
 1.11803  0.707107  0.5  0.707107

julia> rr2(Int, (3, 3), offset=(1, 1), scale=(10, 10))
3×3 IndexFunArray{Int64, 2, IndexFunArrays.var"#4#5"{Int64, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}:
   0  100  400
 100  200  500
 400  500  800
```
## Application to selected dimensions
Note that the code below yields a 3D array but with a one-sized trailing dimension. This can then be used for broadcasting.
```jldoctest
julia> x = ones(5,6,5);

julia> y=rr2(selectsizes(x,(1,2)))
5×6×1 IndexFunArray{Float64, 3, IndexFunArrays.var"#4#5"{Float64, Tuple{Float64, Float64, Float64}, Tuple{Int64, Int64, Int64}}}:
[:, :, 1] =
 13.0  8.0  5.0  4.0  5.0  8.0
 10.0  5.0  2.0  1.0  2.0  5.0
  9.0  4.0  1.0  0.0  1.0  4.0
 10.0  5.0  2.0  1.0  2.0  5.0
 13.0  8.0  5.0  4.0  5.0  8.0
```
Similarly you can also use dimensions 2 and 3 yielding an array of `size(y) == (1,6,5)`. 
Note that the necessary modification to the `Base.size` function is currently provided by this package.

## Using List-Mode Arguments
The code below generates 160 Voronoi cells at random positions. The `accumulator` was set to  mimimum
yielding in each pixel the square distance to the closest Voronoi center. See `gaussian` for another example
of using list-mode arguments.
```jldoctest
julia> y = rr2((1000,1000),offset = (1000.0,1000.0) .* rand(2,160), accumulator=minimum);

```
    
---
    rr2(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`rr2(eltype(arr), size(arr), scaling=scaling, offset=offset)`.
"""
rr2


"""
    rr([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

See `rr2` for a description of all options.

# Examples
```jldoctest
julia> rr((3, 3))
3×3 IndexFunArray{Float64, 2, IndexFunArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 1.41421  1.0  1.41421
 1.0      0.0  1.0
 1.41421  1.0  1.41421

julia> rr((3, 3), offset=CtrCorner)
3×3 IndexFunArray{Float64, 2, IndexFunArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 0.0  1.0      2.0
 1.0  1.41421  2.23607
 2.0  2.23607  2.82843
```
---
    rr(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`rr(eltype(arr), size(arr), scaling=scaling, offset=offset)`.
"""
rr


"""
    xx([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

A distance ramp along first dimension.
See `rr2` for a description of all options.
    ```jldoctest
julia> xx((4,4))
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#14#15"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 -2.0  -2.0  -2.0  -2.0
 -1.0  -1.0  -1.0  -1.0
  0.0   0.0   0.0   0.0
  1.0   1.0   1.0   1.0
```
---
    xx(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`xx(eltype(arr), size(arr), scaling=scaling, offset=offset)`.
"""
xx


"""
    yy([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

A distance ramp along second dimension.
See `rr2` for a description of all options.
    ```jldoctest
julia> yy((4,4))
4×4 IndexFunArray{Float64, 2, IndexFunArrays.var"#19#20"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 -2.0  -1.0  0.0  1.0
 -2.0  -1.0  0.0  1.0
 -2.0  -1.0  0.0  1.0
 -2.0  -1.0  0.0  1.0
```
---
    yy(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`yy(eltype(arr), size(arr), scaling=scaling, offset=offset)`.
"""
yy



"""
    zz([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

A distance ramp along third dimension.
See `rr2` for a description of all options.
    ```jldoctest
julia> zz((1, 1, 4))
1×1×4 IndexFunArray{Float64, 3, IndexFunArrays.var"#24#25"{Float64, Tuple{Float64, Float64, Float64}, Tuple{Int64, Int64, Int64}}}:
[:, :, 1] =
 -2.0

[:, :, 2] =
 -1.0

[:, :, 3] =
 0.0

[:, :, 4] =
 1.0
```

---
    zz(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`zz(eltype(arr), size(arr), scaling=scaling, offset=offset)`.
"""
zz

"""
    ee([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

A distance ramp along forth (element) dimension. This dimension is often used as a color or wavelength channel.
See `rr2` for a description of all options.
    ```jldoctest
julia> ee((1, 1, 1, 4))
1×1×1×4 IndexFunArray{Float64, 4, IndexFunArrays.var"#60#63"{Float64, NTuple{4, Float64}}}:
[:, :, 1, 1] =
 -2.0

[:, :, 1, 2] =
 -1.0

[:, :, 1, 3] =
 0.0

[:, :, 1, 4] =
 1.0
```

---
    ee(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`ee(eltype(arr), size(arr), scaling=scaling, offset=offset)`.
"""
ee

"""
    tt([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

A distance ramp along fifth (time) dimension.
See `rr2` for a description of all options.
    ```jldoctest
julia> tt((1, 1, 1, 1, 4))
1×1×1×1×4 IndexFunArray{Float64, 5, IndexFunArrays.var"#69#72"{Float64, NTuple{5, Float64}}}:
[:, :, 1, 1, 1] =
 -2.0

[:, :, 1, 1, 2] =
 -1.0

[:, :, 1, 1, 3] =
 0.0

[:, :, 1, 1, 4] =
 1.0
```

---
    tt(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`tt(eltype(arr), size(arr), scaling=scaling, offset=offset)`.
"""
tt

"""
    phiphi([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

An azimutal spiral phase ramp using atan(). The azimuthal phase spans dimensions 1 and 2.
See `rr2` for a description of all options.
    ```jldoctest
julia> phiphi((5, 5))
5×5 IndexFunArray{Float64, 2, IndexFunArrays.var"#29#30"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 -2.35619   -2.67795   3.14159  2.67795   2.35619
 -2.03444   -2.35619   3.14159  2.35619   2.03444
 -1.5708    -1.5708    0.0      1.5708    1.5708
 -1.10715   -0.785398  0.0      0.785398  1.10715
 -0.785398  -0.463648  0.0      0.463648  0.785398
```

---
    phiphi(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`phiphi(eltype(arr), size(arr), scaling=scaling, offset=offset)`.
"""
phiphi 



"""
    ramp(::Type{T}, dim::Int, dim_size::Int;
        offset=CtrFT, scale=ScaUnit,
        weight=1,
        accumulator=sum) where {T}

Generates a dim-dimensional ramp of size `dim_size` to be used for broadcasting through multidimensional expressions.
`dim` highest dimension of the oriented array to be generated. This is also the ramp direction.
`dim_size` size along this dimension.

For details about offset and scale and dims see rr2.
```jldoctest
julia> ramp(Float32, 1, 7; offset=(2,))
7-element IndexFunArray{Float32, 1, IndexFunArrays.var"#434#435"{Tuple{Int64}, Tuple{Int64}, Int64}}:
 -1
  0
  1
  2
  3
  4
  5
```
ramp(dim::Int, dim_size::Int; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`ramp(Float64, dim::Int, dim_size::Int; offset=CtrFt, scaling=ScaUnit)`.

"""
ramp

"""
    delta([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

A `delta` peak positioned at offset. See `rr2()` for a description of all options. 
Note that `scale` does not influence the result. Also note that this function operates on
a comparison for equality, which means a sub-pixel offset of the delta results into zero.
```jldoctest
julia> delta((5,5))
5×5 IndexFunArray{Float64, 2, IndexFunArrays.var"#79#81"{Float64, Tuple{Float64, Float64}}}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> delta((6,6))
6×6 IndexFunArray{Float64, 2, IndexFunArrays.var"#79#81"{Float64, Tuple{Float64, Float64}}}:
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0

 julia> delta((5,5),offset=CtrCorner)
 5×5 IndexFunArray{Float64, 2, IndexFunArrays.var"#79#81"{Float64, Tuple{Float64, Float64}}}:
  1.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0
  
```

---
    delta(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`delta(eltype(arr), size(arr), scaling=scaling, offset=offset)`.
"""
delta

"""
    gaussian([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        sigma=1.0,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

A gaussian peak positioned at offset. Note that the gaussian is NOT normalized by its integral, but by its maximum.
For a version with a normalized integral, see `normal`.
List-mode is supported also for the argument sigma. 
See the final example below which generates 60 Gaussians at random positions with random strengths and width along X and Y.

# Arguments:
* `sigma`: the (standard deviation) width of the Gaussian. If a tuple is supplied, each entry is interpreted as the width along the correspondin dimension. 
* `offset`: the center position of the Gaussian. You can use a tuple or the indicators `CtrCorner`, `CtrEnd`, `CtrFT`, `CtrRFT` etc.
* `scale`: the scale of the pixel. By default `ScaUnit` is assumed
* `dims`: the dimensions over which to apply this function to.
* `weight`: the strength of the result. Supports list-mode (see rr2 for documentation)
* `accumulator`: the method used for superimposing list-mode data. Only applies in list-mode
```jldoctest
julia> gaussian((5,5))
5×5 IndexFunArray{Float64, 2, IndexFunArrays.var"#100#102"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 0.0183156  0.082085  0.135335  0.082085  0.0183156
 0.082085   0.367879  0.606531  0.367879  0.082085
 0.135335   0.606531  1.0       0.606531  0.135335
 0.082085   0.367879  0.606531  0.367879  0.082085
 0.0183156  0.082085  0.135335  0.082085  0.0183156

julia> gaussian((6,6), sigma=5.0)
 6×6 IndexFunArray{Float64, 2, IndexFunArrays.var"#100#102"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
  0.165299  0.272532  0.367879  0.40657   0.367879  0.272532
  0.272532  0.449329  0.606531  0.67032   0.606531  0.449329
  0.367879  0.606531  0.818731  0.904837  0.818731  0.606531
  0.40657   0.67032   0.904837  1.0       0.904837  0.67032
  0.367879  0.606531  0.818731  0.904837  0.818731  0.606531
  0.272532  0.449329  0.606531  0.67032   0.606531  0.449329

julia> gaussian(Float32,(5,5),offset=CtrCorner)
  5×5 IndexFunArray{Float32, 2, IndexFunArrays.var"#100#102"{Float32, Tuple{Float64, Float64}, Tuple{Float32, Float32}}}:
   1.0          0.606531     0.135335    0.011109    0.000335463
   0.606531     0.367879     0.082085    0.00673795  0.000203468
   0.135335     0.082085     0.0183156   0.00150344  4.53999f-5
   0.011109     0.00673795   0.00150344  0.00012341  3.72665f-6
   0.000335463  0.000203468  4.53999f-5  3.72665f-6  1.12535f-7

julia> y = gaussian((100,100),offset = (100,100) .* rand(2,60), weight=rand(60), sigma=2.0 .*(0.3 .+rand(2,60)));

```

---
gaussian(arr::AbstractArray; offset=CtrFt, sigma=1.0, scaling=ScaUnit)

This is a wrapper for 
`gaussian(eltype(arr), size(arr), sigma=sigma, scaling=scaling, offset=offset)`.
"""
gaussian

"""
    normal([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        sigma=1.0,
        dims=ntuple(+, N),
        scale=ScaUnit,
        weight=1,
        accumulator=sum)

A gaussian peak positioned at offset. Note that this normal distribution (Gaussian) is normalized by its integral.
For a version with normalized to the maximum, see `gaussian`.

# Arguments:
* `sigma`: the (standard deviation) width of the Gaussian
* `offset`: the center position of the Gaussian. You can use a tuple or the indicators `CtrCorner`, `CtrEnd`, `CtrFT`, `CtrRFT` etc.
* `scale`: the scale of the pixel. By default `ScaUnit` is assumed
* `dims`: the dimensions over which to apply this function to.
* `weight`: the strength of the result. Supports list-mode (see rr2 for documentation)
* `accumulator`: the method used for superimposing list-mode data. Only applies in list-mode
```jldoctest
julia> normal((5,5))
5×5 IndexFunArray{Float64, 2, IndexFunArrays.var"#107#109"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 0.00291502  0.0130642  0.0215393  0.0130642  0.00291502
 0.0130642   0.0585498  0.0965324  0.0585498  0.0130642
 0.0215393   0.0965324  0.159155   0.0965324  0.0215393
 0.0130642   0.0585498  0.0965324  0.0585498  0.0130642
 0.00291502  0.0130642  0.0215393  0.0130642  0.00291502

julia> sum(normal((5,5)))
 0.9818147610543744

julia> sum(normal((25,25)))
 1.000000010701151

julia> sum(normal((55,55), sigma=(5.0,2.0)))
0.9999999639340563
```

---
normal(arr::AbstractArray; offset=CtrFt, sigma=1.0, scaling=ScaUnit)

This is a wrapper for 
`normal(eltype(arr), size(arr), sigma=sigma, scaling=scaling, offset=offset)`.
"""
normal

"""
    exp_ikx([T=Float64], size::NTuple{N, Int};
        offset=CtrFT,
        shift_by=size.÷2
        dims=ntuple(+, N),
        scale=ScaFT,
        weight=1,
        accumulator=sum)

A complex-valued phase ramp according to `exp(-2pi i <k,x>)`. If applied as a multiplicative factor in Fourier space,
it will lead to a shift of `x` pixels in real space. Note that this effect is actually realized by a change to the scaling parameter.
The default shift is `size.÷2` which corresponds to `CtrFT`, however, the Ctr... arguments cannot be used for `shift_by`.

The argument `shift_by` supports list-mode, which can be used to conveniently perform multiple shifts simulatneously.
See the final example below, which generates delta peaks at randomized subpixel positions.

# Arguments:
* `offset`: the center position of the Gaussian. You can use a tuple or the indicators `CtrCorner`, `CtrEnd`, `CtrFT`, `CtrRFT` etc.
* `shift_by`: the amount to shift by in real space.
* `scale`: the scale of the pixel. By default `ScaUnit` is assumed
* `dims`: the dimensions over which to apply this function to.
* `weight`: the strength of the result. Supports list-mode (see rr2 for documentation)
* `accumulator`: the method used for superimposing list-mode data. Only applies in list-mode
```jldoctest
julia> a = rr((4,3),offset=CtrCorner)
4×3 IndexFunArray{Float64, 2, IndexFunArrays.var"#37#39"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
 0.0  1.0      2.0
 1.0  1.41421  2.23607
 2.0  2.23607  2.82843
 3.0  3.16228  3.60555

julia> using FourierTools; real(iffts(ffts(a).*exp_ikx(a, shift_by=(2.0,1.0))))
4×3 Matrix{Float64}:
 2.82843   2.0          2.23607
 3.60555   3.0          3.16228
 2.0      -2.22045e-16  1.0
 2.23607   1.0          1.41421

 julia> using FourierTools; y = real(ift(exp_ikx((101,101),weight=rand(60), shift_by=101.0 .*rand(2,60))));

```

---
    exp_ikx(arr::AbstractArray; offset=CtrFt, shift_by==size(arr).÷2, scaling=ScaUnit)

This is a wrapper for 
`exp_ikx(eltype(arr), size(arr), shift_by=shift_by, scaling=scaling, offset=offset)`.
"""
exp_ikx

"""
    propagator(size::NTuple{N, Int}; Δz=1.0, shift_by=(0,0), k_max=0.5, offset=CtrFT, scale=ScaFT, dims=ntuple(+, N), accumulator=sum, weight=1) where {N,T}

    A complex-valued phase which describes the `k_z= Δz * sqrt(1-k_x^2-k_y^2)` phase change upon free space propagation. 
    
    The argument `shift_by` and `Δz` arguments support list-mode, which can be used to conveniently perform multiple shifts simulatneously.
    
    # Arguments:
    * `offset`: the center position of the Gaussian. You can use a tuple or the indicators `CtrCorner`, `CtrEnd`, `CtrFT`, `CtrRFT` etc.
    * `Δz`: the amount to propagate along the third dimension z by in real space.
    * `shift_by`: the amount to shift by in x and y spatial direct in real space.
    * `k_max`: indicates the sampling relative to the Nyquist frequency. In optics this should be `pixelpitch./λ` with the `pixelpitch` as a tuple and the wavelength `λ = n*λ₀` in the medium.
    * `scale`: the scale of the pixel. By default `ScaUnit` is assumed
    * `dims`: the dimensions over which to apply this function to.
    * `weight`: the strength of the result. Supports list-mode (see rr2 for documentation)
    * `accumulator`: the method used for superimposing list-mode data. Only applies in list-mode
---
propagator(arr::AbstractArray; Δz=1.0, shift_by=(0,0), k_max=0.5, offset=CtrFt, scaling=ScaUnit)

This is a wrapper for 
`propagator(eltype(arr), size(arr), Δz=Δz, shift_by=shift_by, k_max=k_max, scaling=scaling, offset=offset)`.

"""
propagator

"""
    window_linear([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional (separable) window with a linear transition from zero at the borders (`border_out`) to one (`border_in`).
```jldoctest
julia> window_linear((8,9),border_in=0.0)
8×9 IndexFunArray{Float64, 2, IndexFunArrays.var"#34#35"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Float64, Float64}}:
 0.0  0.0     0.0    0.0     0.0   0.0     0.0    0.0     0.0
 0.0  0.0625  0.125  0.1875  0.25  0.1875  0.125  0.0625  0.0
 0.0  0.125   0.25   0.375   0.5   0.375   0.25   0.125   0.0
 0.0  0.1875  0.375  0.5625  0.75  0.5625  0.375  0.1875  0.0
 0.0  0.25    0.5    0.75    1.0   0.75    0.5    0.25    0.0
 0.0  0.1875  0.375  0.5625  0.75  0.5625  0.375  0.1875  0.0
 0.0  0.125   0.25   0.375   0.5   0.375   0.25   0.125   0.0
 0.0  0.0625  0.125  0.1875  0.25  0.1875  0.125  0.0625  0.0
```

---
    window_linear(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit,
                                      border_in=0.8, border_out=1.0)

This is a wrapper for 
`window_linear(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_linear

"""
    window_radial_linear([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional radial window with a linear transition from zero at the borders (`border_out`) to one (`border_in`).
Note that `border_in` and `border_out` need to be scalar values for `radial` type windows. The elypticity can be adjusted via the `scale` parameter.
With the default offset and scale the borders are specified relative to the edge.
```jldoctest
julia> window_radial_linear((4,5),border_in=0.0)
4×5 IndexFunArray{Float64, 2, IndexFunArrays.var"#59#60"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Float64, Float64}}:
 0.0  0.0       0.0  0.0       0.0
 0.0  0.292893  0.5  0.292893  0.0
 0.0  0.5       1.0  0.5       0.0
 0.0  0.292893  0.5  0.292893  0.0
```

---
    window_radial_linear(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit,
                         border_in=0.8, border_out=1.0, dims=ntuple(+, N))

This is a wrapper for 
`window_radial_linear(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_radial_linear


"""
    window_edge([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional (separable) window with a sudden transition half way between the borders (`border_out`) to one (`border_in`).
See `?window_linear` for more details on the arguments.

---
    window_edge(arr::AbstractArray; offset=CtrFt, scaling=ScaUnit,
                                      border_in=0.8, border_out=1.0)

This is a wrapper for 
`window_edge(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_edge

"""
    window_radial_edge([T=Float64], size::NTuple; 
                offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional radial window (disk) with a sudden transition half way between the borders (`border_out`) to one (`border_in`).
Note that `border_in` and `border_out` need to be scalar values for `radial` type windows. The elypticity can be adjusted via the `scale` parameter.
See `?window_radial_linear` for more details on the arguments.

---
    window_radial_edge(arr::AbstractArray; 
                       offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

This is a wrapper for 
`window_radial_edge(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_radial_edge

"""
    window_hanning([T=Float64], size::NTuple; 
                       offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional (separable) window with a von Hann transition between the borders (`border_out`) to one (`border_in`).
See `?window_linear` for more details on the arguments.

---
    window_hanning(arr::AbstractArray;
                       offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

This is a wrapper for 
`window_hanning(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_hanning

"""
    window_radial_hanning([T=Float64], size::NTuple; 
                       offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional radial window with a von Hann transition between the borders (`border_out`) to one (`border_in`).
Note that `border_in` and `border_out` need to be scalar values for `radial` type windows. The elypticity can be adjusted via the `scale` parameter.
See `?window_radial_linear` for more details on the arguments.

---
    window_radial_hanning(arr::AbstractArray;
                          offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  
                                      

This is a wrapper for 
`window_radial_hanning(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_radial_hanning

"""
    window_hamming([T=Float64], size::NTuple; 
                   offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional (separable) window with a Hamming transition between the borders (`border_out`) to one (`border_in`).
See `?window_linear` for more details on the arguments.

---
    window_hamming(arr::AbstractArray;
                   offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

This is a wrapper for 
`window_hamming(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_hamming

"""
    window_radial_hamming([T=Float64], size::NTuple; 
                          offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional radial window with a Hamming transition between the borders (`border_out`) to one (`border_in`).
Note that `border_in` and `border_out` need to be scalar values for `radial` type windows. The elypticity can be adjusted via the `scale` parameter.
See `?window_radial_linear` for more details on the arguments.

---
    window_radial_hamming(arr::AbstractArray;
                          offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

This is a wrapper for 
`window_radial_hamming(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_radial_hamming

"""
    window_blackman_harris([T=Float64], size::NTuple; 
                          offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional (separable) window with a  transition according to Blackman/Harris between the borders (`border_out`) to one (`border_in`).
See `?window_linear` for more details on the arguments.

---
    window_blackman_harris(arr::AbstractArray;
                          offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

This is a wrapper for 
`window_blackman_harris(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_blackman_harris

"""
    window_radial_blackman_harris([T=Float64], size::NTuple; 
                          offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional radial window with a Hamming transition according to Blackman/Harris between the borders (`border_out`) to one (`border_in`).
Note that `border_in` and `border_out` need to be scalar values for `radial` type windows. The elypticity can be adjusted via the `scale` parameter.
See `?window_radial_linear` for more details on the arguments.

---
    window_radial_blackman_harris(arr::AbstractArray;
                                  offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))
This is a wrapper for 
`window_radial_blackman_harris(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_radial_blackman_harris

"""
    window_gaussian([T=Float64], size::NTuple; 
                          offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional (separable) window with a  transition according to a decaying Gaussian between the borders (`border_out`) to one (`border_in`).
By default the standard deviation `sigma` of the Gaussian is adjusted such that the `2 sigma` level is reached at `border_out`.
However, this window is not clipped at the outer border, thus allowing the sigma to be adjusted by placing `border_out` closer to `border_in`.
See `?window_linear` for more details on the arguments.
```jldoctest
julia> w1 = window_gaussian((9,9), border_in=(0.3,0.3), border_out=(0.6,1))
9×9 IndexFunArray{Float64, 2, IndexFunArrays.var"#286#288"{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Tuple{Float64, Int64}}}:
 0.000109245  0.000259904  0.000413188  0.000449917  0.000449917  0.000449917  0.000413188  0.000259904  0.000109245
 0.012239     0.0291178    0.0462907    0.0504055    0.0504055    0.0504055    0.0462907    0.0291178    0.012239
 0.152725     0.363345     0.577637     0.628984     0.628984     0.628984     0.577637     0.363345     0.152725
 0.242811     0.57767      0.918365     1.0          1.0          1.0          0.918365     0.57767      0.242811
 0.242811     0.57767      0.918365     1.0          1.0          1.0          0.918365     0.57767      0.242811
 0.242811     0.57767      0.918365     1.0          1.0          1.0          0.918365     0.57767      0.242811
 0.152725     0.363345     0.577637     0.628984     0.628984     0.628984     0.577637     0.363345     0.152725
 0.012239     0.0291178    0.0462907    0.0504055    0.0504055    0.0504055    0.0462907    0.0291178    0.012239
 0.000109245  0.000259904  0.000413188  0.000449917  0.000449917  0.000449917  0.000413188  0.000259904  0.000109245
```
---
window_gaussian(arr::AbstractArray;
                          offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

This is a wrapper for 
`window_gaussian(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_gaussian

"""
    window_radial_gaussian([T=Float64], size::NTuple; 
                          offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))  

A multidimensional radial window with a Gaussian transition between the borders (`border_out`) to one (`border_in`).
By default the standard deviation `sigma` of the Gaussian is adjusted such that the `2 sigma` level is reached at `border_out`.
However, this window is not clipped at the outer border, thus allowing the sigma to be adjusted by placing `border_out` closer to `border_in`.
Note that `border_in` and `border_out` are scalar values for `radial` type windows. The elypticity can be adjusted via the `scale` parameter.
See `?window_radial_linear` for more details on the arguments.

---
window_radial_gaussian(arr::AbstractArray;
                                  offset=CtrFT, scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N))
This is a wrapper for 
`window_radial_gaussian(eltype(arr), size(arr), scaling=scaling, offset=offset, border_in=border_in, border_out=border_out)`.
"""
window_radial_gaussian
