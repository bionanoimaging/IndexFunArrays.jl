# IndexFunArrays.jl

Here you can find the docstrings of all functions.
We also provide several concrete generators.


## IndexFunArray Interface

The abstract `IndexFunArray` definition
```@docs
IndexFunArray
```

## Indices with certain type
```@docs
idx
cpx
```

## Helpful Array Functions
In addition to normal `size` one can imagine a `selectsizes` which returns the sizes
of several dimensions simultaneously.

```@docs
selectsizes
IndexFunArrays.single_dim_size
```
