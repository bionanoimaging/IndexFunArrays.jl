### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ a565c4ee-d425-11eb-0a5c-d7d322fb95c6
using IndexFunArrays, Plots, Colors, ImageShow, PlutoUI

# ╔═╡ eb26d519-61f9-49c8-bf13-fb0bd5b01be6
md"#### IndexFunArrays.jl | Fun with indices (and functions on them)

Main features:
* Arrays without memory allocations
* entries are calculated once accessed

"

# ╔═╡ 44af166f-a01f-43c2-a774-9a09f5d42712
hello = IndexFunArray(x -> "hello", (3, 3))

# ╔═╡ da4c455d-f657-4fc8-b8d5-ac1497e29071
struct IndexFunArray_demo{T, N, F} <: AbstractArray{T, N} where {F}
    # stores the generator function to be applied to the indices. 
    generator::F     
    # output size of the array 
    size::NTuple{N, Int}
end

# ╔═╡ 3a2ddf05-0d2a-401e-bbad-adccf7ab21bf
hello[1, 3]

# ╔═╡ 5d5835da-4331-4f06-a24c-e60177f38cba
img = IndexFunArray(x -> x[1] + x[2], (3, 3))

# ╔═╡ df51963f-ceea-4279-aa51-6d8c8afaa9a0
distance_to_corner = IndexFunArray(x -> sqrt((x[1]-1)^2 + (x[2]-1)^2), (3, 3))

# ╔═╡ c0eb6cff-a47c-485e-8983-2bf4070d5223
md"## Defining apertures in 1 LOC"

# ╔═╡ d4e01cd2-3eff-404f-86c6-fcd4b772f478
@bind radius Slider(0:1:200, show_value=true)

# ╔═╡ a9ab1f95-42d2-4d64-b015-423ec4fcfc8a
@bind c_ind Select(["1" => "CtrFT", "2" => "CtrEnd", "3" => "CtrFFT"])

# ╔═╡ 04b7c231-0868-493e-99ef-3eb4d89b7ad4
ctr = [CtrFT, CtrEnd, CtrFFT][parse(Int, c_ind)];

# ╔═╡ 4512c50f-e094-4ef9-a177-27a9fdc58f09
Gray.(rr((300, 300), offset=ctr) .< radius)

# ╔═╡ 2d23b2d8-5dac-422a-804d-e2ccd27de12e
Gray.(0.5 .+ xx((300, 300), offset=ctr, scale=ScaNorm))

# ╔═╡ 465522c9-d209-4476-867f-f3ba1f1563ca
md"## Different kinds of Window functions"

# ╔═╡ fabeb838-fb8e-46b9-ad97-9271d4ada93b
md"
$(@bind border_out Slider(0:0.01:1, show_value=true))

$(@bind border_in Slider(0:0.01:1, show_value=true))

"

# ╔═╡ 54023a37-4843-48a2-bb94-262590f310c3
begin
	w_gauss = window_gaussian((256, 256); border_in, border_out)
	w_linear = window_linear((256, 256); border_in, border_out)
end;

# ╔═╡ fedff1f5-5595-4ac9-b53f-a11bd864554c
[Gray.(w_gauss) Gray.(w_linear) ]

# ╔═╡ Cell order:
# ╠═a565c4ee-d425-11eb-0a5c-d7d322fb95c6
# ╟─eb26d519-61f9-49c8-bf13-fb0bd5b01be6
# ╠═44af166f-a01f-43c2-a774-9a09f5d42712
# ╠═da4c455d-f657-4fc8-b8d5-ac1497e29071
# ╠═3a2ddf05-0d2a-401e-bbad-adccf7ab21bf
# ╠═5d5835da-4331-4f06-a24c-e60177f38cba
# ╠═df51963f-ceea-4279-aa51-6d8c8afaa9a0
# ╟─c0eb6cff-a47c-485e-8983-2bf4070d5223
# ╟─d4e01cd2-3eff-404f-86c6-fcd4b772f478
# ╟─a9ab1f95-42d2-4d64-b015-423ec4fcfc8a
# ╟─04b7c231-0868-493e-99ef-3eb4d89b7ad4
# ╠═4512c50f-e094-4ef9-a177-27a9fdc58f09
# ╠═2d23b2d8-5dac-422a-804d-e2ccd27de12e
# ╟─465522c9-d209-4476-867f-f3ba1f1563ca
# ╟─fabeb838-fb8e-46b9-ad97-9271d4ada93b
# ╠═54023a37-4843-48a2-bb94-262590f310c3
# ╠═fedff1f5-5595-4ac9-b53f-a11bd864554c
