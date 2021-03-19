# GeneratorArrays.jl
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


## Why this package?
In image processing and other applications you often encounter position-dependent functions some of which can be a bit of work to code.
It helps the thinking to picture such functions as arrays, which contain the index-dependent values. A good examples are windowing functions.
Another more complicated example is a complex-values free-space propagator.
Yet storing such arrays can be memory intensive and slow and one would ideally perform such calculations "on-the-fly", e.g. only when applying the filter
to the Fourier-transformation. Julia has a great mechanism for this: syntactic loop fusion and broadcasting (e.g. using ".*").
Using CartesianIndices() it is possible to write such index-expressions yet they do not "feel" like arrays.


[docs-dev-img]: https://img.shields.io/badge/docs-dev-pink.svg 
[docs-dev-url]: https://bionanoimaging.github.io/GeneratorArrays.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-darkgreen.svg 
[docs-stable-url]: https://bionanoimaging.github.io/GeneratorArrays.jl/stable/

[CI-img]: https://github.com/bionanoimaging/GeneratorArrays.jl/actions/workflows/ci.yml/badge.svg
[CI-url]: https://github.com/bionanoimaging/GeneratorArrays.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/bionanoimaging/GeneratorArrays.jl/branch/master/graph/badge.svg?token=P0YYCPKXI1
[codecov-url]: https://codecov.io/gh/bionanoimaging/GeneratorArrays.jl

