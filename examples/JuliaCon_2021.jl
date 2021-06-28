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
using IndexFunArrays, Plots, Colors, ImageShow, PlutoUI, MappedArrays

# ╔═╡ eb26d519-61f9-49c8-bf13-fb0bd5b01be6
md"#### IndexFunArrays.jl | Fun with indices (and functions on them)
* _Rainer Heintzmann, Leibniz Institute of Photonic Technology in Jena, Germany_
* _Felix Wechsler, Leibniz Institute of Photonic Technology in Jena, Germany_

Main features:
* Arrays without memory allocations
* entries are calculated once accessed
* feels like an arrays but is some kind of generator

"

# ╔═╡ 840c51ce-1300-41b0-be75-4b298c6c3390
md"### Motivation

* Calculation with indices works efficiently with `CartesianIndices`
* However, as you can see it can get complicated
"

# ╔═╡ 1d880f4b-71fb-4b83-8d23-4a46a1481423
function rr_cartesian_indices(s)
	# entry wise function
	rr_fun(index) = √(sum(abs2.(Tuple(index) .- (s .÷ 2 .+ 1))));
	# array of tuple indices
	arr_indices = CartesianIndices(s)
	# apply rr_fun elementwise
	rr_fun.(arr_indices)
end

# ╔═╡ ab66c67d-832b-40ff-bbc2-f52761912b2d
rr_cartesian_indices((5,6))

# ╔═╡ f09a33fa-dd83-4dbd-948b-78303d3fb76d
md"### MappedArrays
* in principle offers similar functionality
* but doesn't provide the explicit functions (like `rr2`, `rr`)
"

# ╔═╡ e70c0bdb-aa02-4017-ae32-ba3b9fb50906
function rr_mapped_arrays(s)
	# entry wise function
	rr_fun(index) = √(sum(abs2.(Tuple(index) .- (s .÷ 2 .+ 1))));
	# array of tuple indices
	mappedarray(rr_fun, CartesianIndices(s))
end

# ╔═╡ bee02ceb-24a4-4d48-bfe0-c24fdfde0472
rr_mapped_arrays((5,6))

# ╔═╡ 01651a42-dce0-4015-a9d0-a92dc02cbd89
rr((5, 6))

# ╔═╡ 3faba71d-ffa3-4bad-96b8-76f488dfc600
md"### Data Structure Behind It"

# ╔═╡ da4c455d-f657-4fc8-b8d5-ac1497e29071
struct IndexFunArray_demo{T, N, F} <: AbstractArray{T, N} where {F}
    # stores the generator function to be applied to the indices. 
    generator::F     
    # output size of the array 
    size::NTuple{N, Int}
end

# ╔═╡ 5b7b9e0b-9f93-43ed-aba5-56cbf9f1bea8
Base.getindex(A::IndexFunArray{T,N}, I::Vararg{Int, N}) where {T,N} = 
    return A.generator(I)

# ╔═╡ 8604d201-c4ba-4314-bc0b-df73f8d3e24a
md"### How To Use It"

# ╔═╡ 44af166f-a01f-43c2-a774-9a09f5d42712
hello = IndexFunArray(x -> "hello", (3, 3))

# ╔═╡ 3a2ddf05-0d2a-401e-bbad-adccf7ab21bf
hello[2, 3]

# ╔═╡ 5d5835da-4331-4f06-a24c-e60177f38cba
img = IndexFunArray(x -> x[1] + x[2], (3, 3))

# ╔═╡ df51963f-ceea-4279-aa51-6d8c8afaa9a0
distance_to_corner = IndexFunArray(x -> sqrt((x[1]-1)^2 + (x[2]-1)^2), (3, 3))

# ╔═╡ dc58b7a5-8c64-4c3f-aa7b-8261647d2513
md"### What is already defined?

#### Window Functions
* Gaussian 
* Hamming
* Hanning
* Linear
* ...

#### Distance Functions
* radius
* radius squared
* 1D ramps
* ...
"

# ╔═╡ 455d650b-cb51-4dc3-a02f-e6c832e3e5e0


# ╔═╡ d48829c1-c04f-47d3-aa0d-af656b6a5276


# ╔═╡ fbfbc63b-4188-4c03-9bc3-b89b4c6c6298


# ╔═╡ c0eb6cff-a47c-485e-8983-2bf4070d5223
md"## Defining apertures in 1 LOC"

# ╔═╡ d4e01cd2-3eff-404f-86c6-fcd4b772f478
@bind radius Slider(0:1:200, show_value=true)

# ╔═╡ a9ab1f95-42d2-4d64-b015-423ec4fcfc8a
@bind c_ind Select(["1" => "CtrFT", "2" => "CtrEnd", "3" => "CtrFFT"])

# ╔═╡ 04b7c231-0868-493e-99ef-3eb4d89b7ad4
ctr = [CtrFT, CtrEnd, CtrFFT][parse(Int, c_ind)];

# ╔═╡ 4512c50f-e094-4ef9-a177-27a9fdc58f09
[Gray.(rr((300, 300), offset=ctr) .<= radius) Gray.(disc((300, 300), radius, offset=ctr))]

# ╔═╡ d7878c97-068d-4040-9f73-2400e05636bf


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

# ╔═╡ f0c86a16-f255-4b85-b32a-e0adba958504
md"# Gauss"

# ╔═╡ 5c691197-0516-4ab6-9804-0bb6d0df5a66
md"
$(@bind σ1 Slider(0:0.01:200, show_value=true))
$(@bind σ2 Slider(0:0.01:200, show_value=true))
"

# ╔═╡ 7f72a169-5c33-49f4-982e-e85181425374
w_gauss2 = gaussian((256, 256), sigma=(σ1, σ2));

# ╔═╡ 0fa4b07a-5fe2-4294-8753-9adb1ed7c6f7
Gray.(w_gauss2)

# ╔═╡ Cell order:
# ╠═a565c4ee-d425-11eb-0a5c-d7d322fb95c6
# ╟─eb26d519-61f9-49c8-bf13-fb0bd5b01be6
# ╟─840c51ce-1300-41b0-be75-4b298c6c3390
# ╠═1d880f4b-71fb-4b83-8d23-4a46a1481423
# ╠═ab66c67d-832b-40ff-bbc2-f52761912b2d
# ╟─f09a33fa-dd83-4dbd-948b-78303d3fb76d
# ╠═e70c0bdb-aa02-4017-ae32-ba3b9fb50906
# ╠═bee02ceb-24a4-4d48-bfe0-c24fdfde0472
# ╠═01651a42-dce0-4015-a9d0-a92dc02cbd89
# ╟─3faba71d-ffa3-4bad-96b8-76f488dfc600
# ╠═da4c455d-f657-4fc8-b8d5-ac1497e29071
# ╠═5b7b9e0b-9f93-43ed-aba5-56cbf9f1bea8
# ╟─8604d201-c4ba-4314-bc0b-df73f8d3e24a
# ╠═44af166f-a01f-43c2-a774-9a09f5d42712
# ╠═3a2ddf05-0d2a-401e-bbad-adccf7ab21bf
# ╠═5d5835da-4331-4f06-a24c-e60177f38cba
# ╠═df51963f-ceea-4279-aa51-6d8c8afaa9a0
# ╟─dc58b7a5-8c64-4c3f-aa7b-8261647d2513
# ╟─455d650b-cb51-4dc3-a02f-e6c832e3e5e0
# ╟─d48829c1-c04f-47d3-aa0d-af656b6a5276
# ╟─fbfbc63b-4188-4c03-9bc3-b89b4c6c6298
# ╟─c0eb6cff-a47c-485e-8983-2bf4070d5223
# ╟─d4e01cd2-3eff-404f-86c6-fcd4b772f478
# ╟─a9ab1f95-42d2-4d64-b015-423ec4fcfc8a
# ╟─04b7c231-0868-493e-99ef-3eb4d89b7ad4
# ╠═4512c50f-e094-4ef9-a177-27a9fdc58f09
# ╠═d7878c97-068d-4040-9f73-2400e05636bf
# ╠═2d23b2d8-5dac-422a-804d-e2ccd27de12e
# ╟─465522c9-d209-4476-867f-f3ba1f1563ca
# ╟─fabeb838-fb8e-46b9-ad97-9271d4ada93b
# ╠═54023a37-4843-48a2-bb94-262590f310c3
# ╠═fedff1f5-5595-4ac9-b53f-a11bd864554c
# ╟─f0c86a16-f255-4b85-b32a-e0adba958504
# ╟─5c691197-0516-4ab6-9804-0bb6d0df5a66
# ╠═7f72a169-5c33-49f4-982e-e85181425374
# ╠═0fa4b07a-5fe2-4294-8753-9adb1ed7c6f7
