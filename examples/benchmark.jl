using GeneratorArrays

function foo()
    x = ones((10000, 10000))
    @time y = rr2(x) .+ sqrt.(1.2 .* rr2(x));
    @time y .= rr2(x) .+ sqrt.(1.2 .* rr2(x));
    @time y .= rr2(x) .+ sqrt.(1.2 .* rr2(x));
end
