using GeneratorArrays

function foo()
    x = ones(10000,10000);
    qq2(index)=sum(Tuple(index));
    ww2=qq2.(CartesianIndices(x));
    @time y = rr2(x) .+ sqrt.(1.2.*rr2(x).*rr2(x));
    @time y .= rr2(x) .+ sqrt.(1.2.*rr2(x).*rr2(x));
    @time y .= rr2(x) .+ sqrt.(1.2.*rr2(x).*rr2(x));
    @time y .= qq2.(CartesianIndices(x)) .+ sqrt.(1.2.*qq2.(CartesianIndices(x)).*qq2.(CartesianIndices(x)));
    @time y .= qq2.(CartesianIndices(x)) .+ sqrt.(1.2.*qq2.(CartesianIndices(x)).*qq2.(CartesianIndices(x)));
    @time y .= ww2 .+ sqrt.(1.2.*ww2.*ww2);
    @time y .= ww2 .+ sqrt.(1.2.*ww2.*ww2);
end

foo();

