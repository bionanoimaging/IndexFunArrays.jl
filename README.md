# IndexFunArrays.jl
Fun with Indices (and functions on them.)
This package allows to generate complex array expressions which do not allocate memory but instead are generated once they are accessed.


| **Documentation**                       | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][CI-img]][CI-url] | [![][codecov-img]][codecov-url] |



## Installation
Not registered yet,
Type `]`in the REPL to get to the package manager:
```julia
julia> ] add https://github.com/RainerHeintzmann/GeneratorArrays.jl
```


## Quick Examples
```julia
julia> rr2((4,4), offset=CtrMid)  # GeneratorArray containing the square of the radius to the mid position
  4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
   4.5  2.5  2.5  4.5
   2.5  0.5  0.5  2.5
   2.5  0.5  0.5  2.5
   4.5  2.5  2.5  4.5

julia> rr2((3, 3), offset=(1, 1)) # square distance to the top left pixel
  3×3 GeneratorArray{Float64, 2, GeneratorArrays.var"#4#5"{Float64, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}:
   0.0  1.0  4.0
   1.0  2.0  5.0
   4.0  5.0  8.0

julia> rr((4,4), scale=ScaUnit)  # distance (not square) to the Fourier-center with unity pixel scaling
  4×4 GeneratorArray{Float64, 2, GeneratorArrays.var"#9#10"{Float64, Tuple{Float64, Float64}, Tuple{Int64, Int64}}}:
   2.82843  2.23607  2.0  2.23607
   2.23607  1.41421  1.0  1.41421
   2.0      1.0      0.0  1.0
   2.23607  1.41421  1.0  1.41421

julia> GeneratorArray(x -> sum(abs2.(x)), (3, 3))   # directly using the constructor and supplying a function to store in the array
  3×3 GeneratorArray{Int64, 2, var"#184#185"}:
    2   5  10
    5   8  13
   10  13  18
```


## Why this package?
In image processing and other applications you often encounter position-dependent functions some of which can be a bit of work to code.
It helps the thinking to picture such functions as arrays, which contain the index-dependent values. A good examples are windowing functions.
Another more complicated example is a complex-values free-space propagator.
Yet storing such arrays can be memory intensive and slow and one would ideally perform such calculations "on-the-fly", e.g. only when applying the filter
to the Fourier-transformation. Julia has a great mechanism for this: syntactic loop fusion and broadcasting (e.g. using ".*").
Using CartesianIndices() it is possible to write such index-expressions yet they do not "feel" like arrays.
GeneratorArrays allow index-based calculations to look like arrays and to take part in loop fusion. This eases the writing of more complicated expressions without loss in speed
due to Julia's syntactic loop fusion mechanism.
You can think of a GeneratorArray of being an array that stores an expression calculating with indices inside.
This also means you cannot assing to such arrays which also precludes using range indices. However views are possible and range indices can be applied to such views.
Of course such arrays can generate any datatype. See `?GeneratorArray` for more detail.


[docs-dev-img]: https://img.shields.io/badge/docs-dev-pink.svg 
[docs-dev-url]: https://bionanoimaging.github.io/GeneratorArrays.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-darkgreen.svg 
[docs-stable-url]: https://bionanoimaging.github.io/GeneratorArrays.jl/stable/

[CI-img]: https://github.com/bionanoimaging/GeneratorArrays.jl/actions/workflows/ci.yml/badge.svg
[CI-url]: https://github.com/bionanoimaging/GeneratorArrays.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/bionanoimaging/GeneratorArrays.jl/branch/master/graph/badge.svg?token=P0YYCPKXI1
[codecov-url]: https://codecov.io/gh/bionanoimaging/GeneratorArrays.jl

