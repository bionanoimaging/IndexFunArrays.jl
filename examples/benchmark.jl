using GeneratorArrays

function rr_test(s)
    x = ones(s)
    @time y = rr2(x) .+ sqrt.(1.2 .* rr2(x, scale=(1,1), offset=CtrCorner));
    @time y .= rr2(x) .+ sqrt.(1.2 .* rr2(x, scale=(1,1), offset=CtrCorner));
    @time y .= rr2(x) .+ sqrt.(1.2 .* rr2(x, scale=(1,1), offset=CtrCorner));
    @time y .= rr2(x) .+ sqrt.(1.2 .* rr2(x, scale=(1,1), offset=CtrCorner));
    return 
end


function xx_test(s)
    x = ones(s)
    @time y = xx(x) .+ sqrt.(1.2 .* abs.(xx(x)));
    @time y .= xx(x) .+ sqrt.(1.2 .* abs.(xx(x)));
    @time y .= xx(x) .+ sqrt.(1.2 .* abs.(xx(x)));
    return 
end
