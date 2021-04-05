@testset "utils" begin
    @testset "single_dim_size" begin

        @test (1,1,1,3) ==IndexFunArrays.single_dim_size(4, 3)
  
        @test (1,1,1,5) == IndexFunArrays.single_dim_size(4, 5)
  
        @test (1,5) == IndexFunArrays.single_dim_size(2, 5)
    end


    @testset "Test own size method" begin
        x = ones((2,4,6,8, 10));
        @test selectsizes(x, (2, 3, 4)) == (1, 4, 6, 8, 1)
        @test selectsizes(x, (4, 3, 2), keep_dims = false) == (8, 6, 4)
        
        x = ones((10));
        @test selectsizes(x, (1,), keep_dims=false) == (10,)
        @test selectsizes(x, (1,), keep_dims=true) == (10,)
        @test selectsizes(x, 1, keep_dims=false) == (10,)
        @test selectsizes(x, 1, keep_dims=true) == (10,)
        
    
        x = ones((10, 10));
        @test selectsizes(x, (1, 2), keep_dims=false) == (10, 10)
        @test selectsizes(x, (1, 2), keep_dims=true) == (10, 10)
       
        x = ones((1,1,3));
        @test selectsizes(x, (1,2), keep_dims=false) == (1,1)
        @test selectsizes(x, (3,3), keep_dims=false) == (3,3)
        @test selectsizes(x, (3,3), keep_dims=true) == (1, 1, 3)
    
    end



end
