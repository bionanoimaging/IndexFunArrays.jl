
@testset "normal distribution, integral normed to 1" begin
    @test 1 / √(2π)^2 / 1 * exp.(.- rr2((5,5)) / 2 .* 1^2) ≈ normal((5,5))
    @test 1 / √(2π)^2 / 1 * exp.(.- rr2((8,5)) / 2 .* 1^2) ≈ normal((8,5))
    @test 1 / √(2π)^2 / 1 * exp.(.- rr2((8,5)) / 2 .* 1^2) ≈ normal(zeros(8,5))
    @test eltype(normal(Float32, (8,5))) === Float32
end

@testset "Gaussian function, normed to 1" begin
    arr = 1 / √(2π)^2 / 1 * exp.(.- rr2((5,5)) / 2 .* 1^2)
    arr ./= maximum(arr)
    @test arr ≈ gaussian((5,5))
    @test eltype(gaussian(Float32, (8,5))) === Float32
end

@testset "propagator" begin
    p = propagator((5,7))
    @test p[1,1] == 1.0
    @test p[3,4] == 1.0
    @test eltype(p) === ComplexF64
    @test eltype(propagator(Float32, (5,7))) === ComplexF32
end
