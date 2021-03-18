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

    @testset "Test xx and yy method" begin
        a = xx(Float64, (3, 5), offset = CtrCorner, scale = ScaUnit)
        b = yy(Float64, (5, 3), offset = CtrCorner, scale = ScaUnit)'
        @test  a == [0.0 0.0 0.0 0.0 0.0; 1.0 1.0 1.0 1.0 1.0; 2.0 2.0 2.0 2.0 2.0]
#        @test a == b

#        @test xx(Int, (3, 5), offset = CtrCorner, scale = ScaUnit) == [0 0 0 0 0; 1 1 1 1 1; 2 2 2 2 2] 
    



    end

end
