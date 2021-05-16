using View5D, IndexFunArrays

sz = (50,50)
N = 1
offsets = 1 .+ (sz.-1) .* rand(2,N)
# sigmas = 2.0 .*(0.3 .+rand(2,N))
mygaussians(of) = gaussian(sz,offset = of, weight=1.0, sigma=2.0);   # , sigma=sigmas
data =  mygaussians(offsets);
@vv data  # use "e" to advance spectral channels, "n" in the bottom right (element) panel to toggle normalization
@ve  mygaussians(offsets)

loss(of) = sum(abs2.(data .- mygaussians(of)))  # a loss function
loss(offsets)  # is zero

using Zygote, IndexFunArrays
# gradient(loss,offsets)

# f(idx,of,sc) = (idx[1].*sc - of[1] + of[2])^2
function c(of,sc)
    f(idx) = (idx[1].*sc - of[1] + of[2])^2
    IndexFunArray(f,(10,))
end
c((2.2,1.1),1.0)
loss(of,sc) = sum(c(of,sc))
gradient(loss,(2.2,1.3),1.2)


using Zygote, IndexFunArrays
f(of,s) = sum(xx((10,),offset=of, scale=s))  # Float64(idx[1]+idx[2])*
f(2.2,3.3)
gradient(f, 1, 2)


ntuple(i -> i ∈ dims ? scale[i] : zero(scale[1]), N)


# by "r" "g" "b" you can choose colors for elements. "v" toggles the color in and out of overlay,
# color overlay is accessed by "C" (shift-"c")

