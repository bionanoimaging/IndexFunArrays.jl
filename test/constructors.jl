

@testset "Check constructors" begin 
    float_types_to_test = [Float16, Float32, Float64, ComplexF32, ComplexF64]
    int_types_to_test = [Int16, Int32, Int64]
    
    
    function test_arr(arr, f, T, s, offset, scale)
        @test typeof(arr[s...]) == T
        @test last(arr) ≈ sqrt(sum(abs2.(scale .* (s .- offset))))
        @test arr.size == s
        @test arr.offset === offset
        @test arr.scale === scale
        @test first(arr) ≈ sqrt(sum(abs2.(scale .* (ntuple(x -> 1, length(s)) .- offset))))
        @test typeof(arr) <: GeneratorArray{T, length(s), F} where F
    end
    
    
    @testset "Check GeneratorArray initializer with gen and size for types and boundaries" begin
    
        sizes = [(10, 10), (2,2,2), (3,3,3,3)]
        for s in sizes
            for T in float_types_to_test 
                f(ind) = T(sqrt(sum(abs2.(ind))))
                arr = GeneratorArray(f, s)
                ofs = ntuple(x -> T(0), length(s)) 
                scale = ntuple(x -> T(1), length(s)) 
                test_arr(arr, f, T, s, ofs, scale)
                arr = GeneratorArray(T, f, s)
                test_arr(arr, f, T, s, ofs, scale)
            end
        end
    end    
    
    @testset "Check GeneratorArray initializer with gen, size, offset, scale for types and boundaries" begin
    
        sizes = [(10, 10), (2,2,2) ]
        for s in sizes
            for T in float_types_to_test 
                f(ind) = T(sqrt(sum(abs2.(ind))))
                
                offset = ntuple(x -> T(rand(1:10)), length(s))
                scale = ntuple(x -> T(rand(1:10)), length(s))
                
                arr = GeneratorArray(f, offset, scale, s)
                test_arr(arr, f, T, s, offset, scale)
               
                arr = GeneratorArray(T, f, offset, scale, s)
                test_arr(arr, f, T, s, offset, scale)
            end
        end
    end    
    
    @testset "Check GeneratorArray initializer with gen, size and offset  for types and boundaries" begin
    
        sizes = [(10, 10), (2,2,2) ]
        for s in sizes
            for T in float_types_to_test 
                f(ind) = T(sqrt(sum(abs2.(ind))))
                
                offset = ntuple(x -> T(rand(1:10)), length(s))
                scale = ntuple(x -> T(1), length(s))
                
                arr = GeneratorArray(f, offset, s)
                test_arr(arr, f, T, s, offset, scale)
               
                arr = GeneratorArray(T, f, offset, s)
                test_arr(arr, f, T, s, offset, scale)
            end
        end
    end    
end
