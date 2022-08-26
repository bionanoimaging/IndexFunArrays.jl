using IndexFunArrays
using BenchmarkTools

function rr_test(s)
    x = ones(s)
    y = rr2(x) .+ sqrt.(1.2 .* rr2(x, scale=(1,1), offset=CtrCorner));
    @btime $y .= rr2($x) .+ sqrt.(1.2 .* rr2($x, scale=(1,1), offset=CtrCorner));
    @btime $y .= rr2($x) .+ sqrt.(1.2 .* rr2($x, scale=(1,1), offset=CtrCorner));
    @btime $y .= rr2($x) .+ sqrt.(1.2 .* rr2($x, scale=(1,1), offset=CtrCorner));
    return 
end


function xx_test(s)
    x = ones(s)
    y = xx(x) .+ sqrt.(1.2 .* abs.(xx(x)));
    @btime $y .= xx($x) .+ sqrt.(1.2 .* abs.(xx($x)));
    @btime $y .= xx($x) .+ sqrt.(1.2 .* abs.(xx($x)));
    return 
end

function compare_to_CartesianIndices()
    x = ones(1000,1000);
    sz = size(x)

    qq2(index) = sum(Tuple(index).^2);
    rng_start = .-(sz.÷2)
    rng_stop = rng_start .+ sz .-1
    rng = Tuple(sta:sto for (sta,sto) in zip(rng_start, rng_stop))

    ci = CartesianIndices(rng);
    ww2 = qq2.(CartesianIndices(rng));
    y = rr2(x) .+ sqrt.(1.2.*rr2(x).*rr2(x));

    @info "rr2 based"
    @btime $y .= rr2($x) .+ sqrt.(1.2.*rr2($x).*rr2($x));
    @btime $y .= rr2($x) .+ sqrt.(1.2.*rr2($x).*rr2($x));
    @info "CartesianIndices based"
    @btime $y .= ($qq2).(CartesianIndices($rng)) .+ sqrt.(1.2.*($qq2).(CartesianIndices($rng)).*($qq2).(CartesianIndices($rng)));
    @btime $y .= ($qq2).(CartesianIndices($rng)) .+ sqrt.(1.2.*($qq2).(CartesianIndices($rng)).*($qq2).(CartesianIndices($rng)));
    @info "precomputed CartesianIndices based"
    @btime $y .= ($qq2).(ci) .+ sqrt.(1.2.*($qq2).(ci).*($qq2).(ci));
    @btime $y .= ($qq2).(ci) .+ sqrt.(1.2.*($qq2).(ci).*($qq2).(ci));
    @info "Applying to precomputed arrays (here generated via CartesianIndices)"
    @btime $y .= $ww2 .+ sqrt.(1.2.*$ww2.*$ww2);
    @btime $y .= $ww2 .+ sqrt.(1.2.*$ww2.*$ww2);
end 

