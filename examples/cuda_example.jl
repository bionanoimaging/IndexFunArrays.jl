using CUDA
using IndexFunArrays

function f(arr, x)
    arr .= arr.^2 .+ arr.^3 .+ sqrt.(x)
    return arr
end

function f2(arr)
    arr .= arr.^2 .+ arr.^3 .+ sqrt.(arr .^2) 
    return arr
end

function f3(arr)
    arr = arr.^2 .+ arr.^3 .+ sqrt.(arr .^2) 
    return arr
end

function test(s)
    arr = randn(Float32, s)
    arr_c = CuArray(arr)
    x = rr2(Float32, s)
    x_c = CuArray(rr(Float32, s))

    @info "f on CPU"
    @time f(arr, x);
    @time f(arr, x);
    @info "f2 on CPU"
    @time f2(arr);
    @time f2(arr);
    @info "f on GPU"
    CUDA.@time f(arr_c, x);
    CUDA.@time f(arr_c, x);
    @info "f2 on GPU"
    CUDA.@time f2(arr_c);
    CUDA.@time f2(arr_c);
    @info "f3 on GPU"
    CUDA.@time f3(arr_c);
    CUDA.@time f3(arr_c);

    return 

end
