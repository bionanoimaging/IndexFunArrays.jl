@testset "Test generated concrete generator functions" begin

    @testset "Test rr method" begin
    
    
        @test rr((5,)) == [2.0, 1.0, 0.0, 1.0, 2.0]
        @test rr(Int, (5,)) == [2, 1, 0, 1, 2]
        @test rr((4, 4)) == [2.8284271247461903 2.23606797749979 2.0 2.23606797749979; 2.23606797749979 1.4142135623730951 1.0 1.4142135623730951; 2.0 1.0 0.0 1.0; 2.23606797749979 1.4142135623730951 1.0 1.4142135623730951]
    
        @test rr((2, 3), offset=(1, 1)) == [0.0 1.0 2.0; 1.0 1.4142135623730951 2.23606797749979]
    
        @test rr((2, 3), offset=(1, 1)) == [0.0 1.0 2.0; 1.0 1.4142135623730951 2.23606797749979]
    
        @test rr(Float64, (3, 5), offset = CtrCorner, scale = ScaUnit) == [0.0 1.0 2.0 3.0 4.0; 1.0 1.4142135623730951 2.23606797749979 3.1622776601683795 4.123105625617661; 2.0 2.23606797749979 2.8284271247461903 3.605551275463989 4.47213595499958]
    
        @test rr2(Int, (10, )) == [25; 16; 9; 4; 1; 0; 1; 4; 9; 16]
        @test rr2(Int, (10, 1)) == reshape([25; 16; 9; 4; 1; 0; 1; 4; 9; 16], (10, 1))
    end


    @testset "Test rr2 method for an array as input" begin
        x = rand((1:100), (10,))
        @test rr2(x) == [25; 16; 9; 4; 1; 0; 1; 4; 9; 16]
        x = rand((1:100), (10,1))
        @test rr2(x) == reshape([25; 16; 9; 4; 1; 0; 1; 4; 9; 16], (10, 1))

    end

    @testset "Test xx and yy method" begin
        a = xx(Float64, (3, 5), offset = CtrCorner, scale = ScaUnit)
        b = yy(Float64, (5, 3), offset = CtrCorner, scale = ScaUnit)'
        @test  a == [0.0 0.0 0.0 0.0 0.0; 1.0 1.0 1.0 1.0 1.0; 2.0 2.0 2.0 2.0 2.0]
        @test a == b

        @test xx(Int, (3, 5), offset = CtrCorner, scale = ScaUnit) == [0 0 0 0 0; 1 1 1 1 1; 2 2 2 2 2] 
    
    @testset "Test different scalings" begin
        @test rr((10,), scale = ScaNorm) == [0.5555555555555556, 0.4444444444444444, 0.3333333333333333, 0.2222222222222222, 0.1111111111111111, 0.0, 0.1111111111111111, 0.2222222222222222, 0.3333333333333333, 0.4444444444444444]
        @test rr((10,), scale = ScaFT) == [0.5, 0.4, 0.30000000000000004, 0.2, 0.1, 0.0, 0.1, 0.2, 0.30000000000000004, 0.4]
    
        @test rr((11,), scale = ScaFTEdge) == [1.0, 0.8, 0.6000000000000001, 0.4, 0.2, 0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0]
        
        @test rr((10,), scale = ScaFTEdge) == [1.0, 0.8, 0.6000000000000001, 0.4, 0.2, 0.0, 0.2, 0.4, 0.6000000000000001, 0.8]


        @test rr((11,), scale = (5,)) == [25.0, 20.0, 15.0, 10.0, 5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0]

        @test rr((4,), scale = (5,)) == [10.0, 5.0, 0.0, 5.0]
    end


    end

    # simple hanning test
    @testset "Hanning window test" begin
        hanning((10, ), border_in=0, border_out=1) ≈ sinpi.(range(0, 1, length=11)[1:end-1]).^2
        hanning((10, ), border_in=0.4, border_out=0.499) ≈ [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
        
        x = hanning(Int, (10, ), border_in=0.4, border_out=0.499) 
        @test x ≈ [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
        @test typeof(x) <: GeneratorArray{Int, 1, T} where T
        
        x = hanning(ComplexF32, (10, ), border_in=0.4, border_out=0.499) 
        @test x ≈ [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
        @test typeof(x) <: GeneratorArray{ComplexF32, 1, T} where T
        

        x = hanning(Float32.(randn((10, ))), border_in=0.4, border_out=0.499) 
        @test x ≈ [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
        @test typeof(x) <: GeneratorArray{Float32, 1, T} where T
        
    end


end
