using Test
using Random

using IndexFunArrays

include("test_center_scaling.jl")
include("constructors.jl")
include("concrete_generators.jl")
include("utils.jl")

@testset "Test defined array operation" begin
    f(ind) = ind[1]
    arr = IndexFunArray(Int, f, (10, 10))
    @test arr[2] == arr[12] == arr[2, 2] == 2
    
    arr2 = similar(arr, (11, 11))
    @test arr2[2] == arr2[13] == arr2[2, 2] == 2

    for s in [(10,), (10,1,2), (100, 100), (20,2,2,2,2)]
        @test s == size(IndexFunArray(ComplexF32, x -> zero(ComplexF32), s))
    end

end


@testset "Check errors" begin
    @test_throws ArgumentError IndexFunArray(Int, x -> zero(Float32), (10, 10))
    
    a = IndexFunArray(Int, x -> zero(Int), (10, 10))
    @test_throws ErrorException (a[1] = 1) 
end


return
