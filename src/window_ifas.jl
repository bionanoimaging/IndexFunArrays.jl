
function generate_window_functions_expr()
    # x_exprW = :(clamp.(- border_in .+ abs.(scale .* (x .- offset))./(border_out-border_in),0,1))
    x_exprW = :(clamp.(1 .-(abs.(scale .* (x .- offset)).-border_in)./(border_out .- border_in),0,1)) 
    x_exprRW = :(clamp.(1 .-(sqrt.(sum((scale .* (x .- offset)).^2)).-border_in)./(border_out .- border_in),0,1)) 

    functions = [  # see https://en.wikipedia.org/wiki/Window_function 
        (:(window_linear),  :(x -> T(prod(($x_exprW))))),
        (:(window_edge),  :(x -> T(prod(($x_exprW).>0.5)))),
        (:(window_hanning),  :(x -> T(prod(sinpi.(0.5 .* ($x_exprW)).^2)))),
        (:(window_half_cos),  :(x -> T(prod(sinpi.(0.5 .* ($x_exprW)))))),
        (:(window_hamming),  :(x -> T(prod(0.54 .-0.46.*cospi.(($x_exprW)))))),
        (:(window_blackman_harris),  :(x -> T(prod(0.35875 .- 0.48829.*cospi.($x_exprW).+0.14128.*cospi.(2 .*$x_exprW).-0.01168.*cospi.(3 .*$x_exprW))))),
        (:(window_radial_linear),  :(x -> T($x_exprRW))),
        (:(window_radial_edge),  :(x -> T(($x_exprRW).>0.5))),
        (:(window_radial_hanning),  :(x -> T(sinpi.(0.5 .* ($x_exprRW)).^2))),
        (:(window_radial_hamming),  :(x -> T((0.54 .-0.46.*cospi.(($x_exprRW)))))),
        (:(window_radial_blackman_harris),  :(x -> T((0.35875 .-0.48829.*cospi.($x_exprRW).+0.14128.*cospi.(2 .*$x_exprRW).-0.01168.*cospi.(3 .*$x_exprRW))))),
    ]
    return functions
end

# we automatically generate the functions for different windows like hanning 
for F in generate_window_functions_expr() 
    # default functions with certain offset and scaling behavior
    @eval function $(F[1])(::Type{T}, size::NTuple{N, Int};
                           offset=CtrFT,
                           scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N)) where{N, T} 
        scale_init = get_scale(size, scale)
        scale = ntuple(i -> i ∈ dims ? scale_init[i] : zero(scale_init[1]), N)
        offset = get_offset(size, offset)
        IndexFunArray(T, $(F[2]), size) 
    end
    
    # change order of offset and scale
    @eval function $(F[1])(size::NTuple{N, Int}; 
                           offset=CtrFT, 
                           scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N)) where{N} 
        T = $DEFAULT_T 
        $(F[1])(T, size, scale=scale, offset=offset, border_in=border_in, border_out=border_out, dims=dims) 
    end

    # convenient wrapper to provide an array as input
    @eval function $(F[1])(arr::AbstractArray{T, N}; 
                           offset=CtrFT, 
                           scale=ScaFTEdge, border_in=0.8, border_out=1.0, dims=ntuple(+, N)) where{N, T} 
        $(F[1])(default_type(T, $DEFAULT_T), size(arr), scale=scale, offset=offset, border_in=border_in, border_out=border_out, dims=dims)
    end

    @eval export $(F[1])
end


