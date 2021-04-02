@testset "Check center and scaling methods" begin

    @testset "Check center" begin

        @test IndexFunArrays.get_offset((5,), CtrCorner) == (1.0,)
        @test IndexFunArrays.get_offset((4,), CtrCorner) == (1.0,)
        @test IndexFunArrays.get_offset((5, 2), CtrFT) == (3.0, 2.0)
        @test IndexFunArrays.get_offset((5, 1), CtrFT) == (3.0, 1.0)
        @test IndexFunArrays.get_offset((5, 2), CtrFFT) == (1.0, 1.0)
        
        @test IndexFunArrays.get_offset((5, 2), CtrMid) == (3.0, 1.5)
        @test IndexFunArrays.get_offset((5, 3), CtrMid) == (3.0, 2.0)


        @test IndexFunArrays.get_offset((5, 2), CtrRFT) == (1, 2)
        @test IndexFunArrays.get_offset((5, 3), CtrRFT) == (1, 2)

        @test IndexFunArrays.get_offset((5, 2), CtrRFFT) == (1, 1)
        @test IndexFunArrays.get_offset((5, 3), CtrRFFT) == (1, 1)
    
        @test IndexFunArrays.get_offset((1123, 3123, 2), 123.123) == (123.123, 123.123, 123.123)

        @test IndexFunArrays.get_offset((5, 4), CtrEnd) == (5.0, 4.0)
    end


    @testset "Check scaling" begin
        @test IndexFunArrays.get_scale((5, 4), ScaFT) === (0.25, 0.25)
        @test IndexFunArrays.get_scale((5, 4), ScaUnit) == (1, 1)
        @test IndexFunArrays.get_scale((6, 5), ScaNorm) == (0.2, 0.25)
        @test IndexFunArrays.get_scale((6, 5), ScaNorm) == (0.2, 0.25)
        @test IndexFunArrays.get_scale((5,3), (2,11)) == (2,11)

        @test IndexFunArrays.get_scale((12, 1, 2), 2.123f0) == (2.123f0, 2.123f0, 2.123f0)


    end
end
