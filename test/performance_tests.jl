@testset "Some performance testing with memory allocations" begin

    
    @testset "idx performance" begin
        function f(s, T)
            res = @benchmark idx($T, $s) samples = 2 evals = 2;
            @test res.memory < 500
            
            res = @benchmark exp.(sum.(idx($T, $s))) samples = 2 evals = 2;
            res2 = @benchmark randn($T, $s) samples=2 evals = 2

            
            # maximum 10% overhead in memory
            @test res.memory < res2.memory * 1.1
        end

        f((100, 10), Float32)
        f((256, 256), ComplexF64)
    end

    @testset "rr2 performance" begin
        function f(s, T)
            res = @benchmark rr2($T, $s) samples = 2 evals = 2;
            @test res.memory < 500
            
            res = @benchmark exp.((rr2($T, $s))) samples = 2 evals = 2;
            res2 = @benchmark randn($T, $s) samples=2 evals = 2

            
            # maximum 10% overhead in memory
            @test res.memory < res2.memory * 1.1
        end

        f((256, 256), ComplexF64)
    end

    @testset "exp_ikx performance" begin
        function f(s, T)
            res = @benchmark exp_ikx($T, $s) samples = 2 evals = 2;
            @test res.memory < 1000
            
            res = @benchmark exp.((exp_ikx($T, $s))) samples = 2 evals = 2;
            res2 = @benchmark randn($T, $s) samples=2 evals = 2

            
            # maximum 10% overhead in memory
            @test res.memory < res2.memory * 1.1
        end

        f((256, 256), ComplexF64)
    end

    @testset "cpx performance" begin
        function f(s, T)
            res = @benchmark cpx($T, $s) samples = 2 evals = 2;
            @test res.memory < 1000
            
            res = @benchmark exp.((cpx($T, $s))) samples = 2 evals = 2;
            res2 = @benchmark randn(Complex{$T}, $s) samples=2 evals = 2

            
            # maximum 10% overhead in memory
            @test res.memory < res2.memory * 1.1
        end

        f((256, 256), Float32)
    end

end
