@testset "Test generated concrete generator functions" begin

    @testset "Type conversion for default type" begin
        def_T = Float64

        Types_changed = [Bool, Int8, Int16, Int32, Int64, Int128, UInt8, UInt16, UInt32, UInt64, UInt128]
        Types_unchanged = [ComplexF32, ComplexF64, Complex{Int}, Float32, Float64]

        for T in Types_changed
            IndexFunArrays.default_type(T, Float64) = Float64
            for T2 in Types_changed
                @test IndexFunArrays.default_type(T, def_T) == def_T
            end
        end
        
        for T in Types_unchanged
            for T2 in Types_changed
                @test IndexFunArrays.default_type(T, T2) == T
            end
        end
    end

    @testset "Test rr method" begin
   
        offsets = [(3,3), (5,2)]
        scales = [(1.0,1.0), (2.0,0.5)]
        s = ((7,7))
        @test rr(s, offset=offsets, scale=scales) ≈ reduce(.+, rr(s, offset=o, scale=sc) for (o,sc) in zip(offsets,scales))
        
        @test rr2(s, offset=offsets, scale=scales) ≈ reduce(.+, rr2(s, offset=o, scale=sc) for (o,sc) in zip(offsets,scales))

    
        @test rr((5,)) == [2.0, 1.0, 0.0, 1.0, 2.0]
        @test rr(Int, (5,)) == [2, 1, 0, 1, 2]
        @test rr((4, 4)) == [2.8284271247461903 2.23606797749979 2.0 2.23606797749979; 2.23606797749979 1.4142135623730951 1.0 1.4142135623730951; 2.0 1.0 0.0 1.0; 2.23606797749979 1.4142135623730951 1.0 1.4142135623730951]
    
        @test rr((2, 3), offset=(1, 1)) == [0.0 1.0 2.0; 1.0 1.4142135623730951 2.23606797749979]
    
        @test rr((2, 3), offset=(1, 1)) == [0.0 1.0 2.0; 1.0 1.4142135623730951 2.23606797749979]
    
        @test rr(Float64, (3, 5), offset = CtrCorner, scale = ScaUnit) == [0.0 1.0 2.0 3.0 4.0; 1.0 1.4142135623730951 2.23606797749979 3.1622776601683795 4.123105625617661; 2.0 2.23606797749979 2.8284271247461903 3.605551275463989 4.47213595499958]
    
        @test rr2(Int, (10, )) == [25; 16; 9; 4; 1; 0; 1; 4; 9; 16]
        @test rr2(Int, (10, 1)) == reshape([25; 16; 9; 4; 1; 0; 1; 4; 9; 16], (10, 1))
        @test rr2(Int, (10, 1)) == reshape([25; 16; 9; 4; 1; 0; 1; 4; 9; 16], (10, 1))
    end


    @testset "Test rr2 method for an array as input" begin
        x = rand((1:100), (10,))
        @test rr2(x) == [25; 16; 9; 4; 1; 0; 1; 4; 9; 16]
        x = rand((1:100), (10,1))
        @test rr2(x) == reshape([25; 16; 9; 4; 1; 0; 1; 4; 9; 16], (10, 1))

    end

    @testset "Test xx, yy, zz, ee and tt methods" begin
        a = xx(Float64, (3, 5), offset = CtrCorner, scale = ScaUnit)
        b = yy(Float64, (5, 3), offset = CtrCorner, scale = ScaUnit)
        c = zz(Float64, (1, 5, 3), offset = CtrCorner, scale = ScaUnit)
        d = ee(Float64, (1, 1, 5, 3), offset = CtrCorner, scale = ScaUnit)
        e = tt(Float64, (1, 1, 1, 5, 3), offset = CtrCorner, scale = ScaUnit)
        @test a == [0.0 0.0 0.0 0.0 0.0; 1.0 1.0 1.0 1.0 1.0; 2.0 2.0 2.0 2.0 2.0]
        @test a == b'
        @test (b.+0)[:] == (c.+0)[:]
        @test (b.+0)[:] == (d.+0)[:]
        @test (b.+0)[:] == (e.+0)[:]
        @test size(c) == (1,5,3)
        @test size(d) == (1,1,5,3)
        @test size(e) == (1,1,1,5,3)

        @test xx(Int, (3, 5), offset = CtrCorner, scale = ScaUnit) == [0 0 0 0 0; 1 1 1 1 1; 2 2 2 2 2] 
    end

    @testset "Test idx_min and idx_max methods" begin
        sz = (13, 15)
        a = idx_min(Int, sz , offset = CtrCorner, scale = ScaUnit)
        b = idx_max(Int, sz, offset = CtrCorner, scale = ScaUnit)
        a == min.(xx(sz, offset = CtrCorner), yy(sz, offset = CtrCorner))
        b == max.(xx(sz, offset = CtrCorner), yy(sz, offset = CtrCorner))
    end

    @testset "Test the ramp method" begin
        a = yy(Float64, (1, 5), offset = CtrCorner, scale = ScaUnit)
        b = ramp(Int64, 2, 5, offset = CtrCorner, scale = ScaUnit)
        @test a == [0,1,2,3,4]'
        @test size(ramp(7,7)) == (1,1,1,1,1,1,7)
    end

    @testset "Test delta method" begin
        a = delta(Float64, (3, 5), offset = CtrCorner, scale = ScaUnit)
        @test  a == [1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
    end

    @testset "Test different scalings" begin
        @test rr((10,), scale = ScaNorm) == [0.5555555555555556, 0.4444444444444444, 0.3333333333333333, 0.2222222222222222, 0.1111111111111111, 0.0, 0.1111111111111111, 0.2222222222222222, 0.3333333333333333, 0.4444444444444444]
        @test rr((10,), scale = ScaFT) == [0.5, 0.4, 0.30000000000000004, 0.2, 0.1, 0.0, 0.1, 0.2, 0.30000000000000004, 0.4]
    
        @test rr((11,), scale = ScaFTEdge) == [ 0.9090909090909092, 0.7272727272727273, 0.5454545454545454, 0.36363636363636365, 0.18181818181818182, 0.0, 0.18181818181818182, 0.36363636363636365, 0.5454545454545454, 0.7272727272727273, 0.9090909090909092]
        # 1.0, 0.8, 0.6000000000000001, 0.4, 0.2, 0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0]
        
        @test rr((10,), scale = ScaFTEdge) == [1.0, 0.8, 0.6000000000000001, 0.4, 0.2, 0.0, 0.2, 0.4, 0.6000000000000001, 0.8]

        @test rr((11,), scale = (5,)) == [25.0, 20.0, 15.0, 10.0, 5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0]

        @test rr((4,), scale = (5,)) == [10.0, 5.0, 0.0, 5.0]
    end

    @testset "Test rr for different offsets" begin
        @test rr2(Int, (5, 5), offset = CtrFFT) == [0 1 4 9 16; 1 2 5 10 17; 4 5 8 13 20; 9 10 13 18 25; 16 17 20 25 32]
        @test rr2(Int, (4, 4), offset = CtrFFT) == [0 1 4 9; 1 2 5 10; 4 5 8 13; 9 10 13 18]
        @test rr2((5, 5), offset = CtrMid) == [8.0 5.0 4.0 5.0 8.0; 5.0 2.0 1.0 2.0 5.0; 4.0 1.0 0.0 1.0 4.0; 5.0 2.0 1.0 2.0 5.0; 8.0 5.0 4.0 5.0 8.0]
        @test rr2((4, 4), offset = CtrMid) == [4.5 2.5 2.5 4.5; 2.5 0.5 0.5 2.5; 2.5 0.5 0.5 2.5; 4.5 2.5 2.5 4.5]
        
        @test rr2((4, 4), offset = CtrEnd) == [18.0 13.0 10.0 9.0; 13.0 8.0 5.0 4.0; 10.0 5.0 2.0 1.0; 9.0 4.0 1.0 0.0]
        @test rr2((5, 6), offset = CtrEnd) == [41.0 32.0 25.0 20.0 17.0 16.0; 34.0 25.0 18.0 13.0 10.0 9.0; 29.0 20.0 13.0 8.0 5.0 4.0; 26.0 17.0 10.0 5.0 2.0 1.0; 25.0 16.0 9.0 4.0 1.0 0.0]

    end
end


@testset "Test idx method" begin
    @test idx(Float32, (5,)) == Tuple{Float32}[(-2.0,), (-1.0,), (0.0,), (1.0,), (2.0,)]
    @test idx(Int8, (5, 2)) == Tuple{Int8, Int8}[(-2, -1) (-2, 0); (-1, -1) (-1, 0); (0, -1) (0, 0); (1, -1) (1, 0); (2, -1) (2, 0)] 
    
    @test idx(randn((5,2))) == Tuple{Int8, Int8}[(-2, -1) (-2, 0); (-1, -1) (-1, 0); (0, -1) (0, 0); (1, -1) (1, 0); (2, -1) (2, 0)] 
end

@testset "Test cpx method" begin
    @test cpx(Float32, (5,)) == [-2.0f0 + 0.0f0im, -1.0f0 + 0.0f0im, 0.0f0 + 0.0f0im,1.0f0 + 0.0f0im,2.0f0 + 0.0f0im]
    @test cpx(Int8, (5, 2)) == Complex{Int8}[(-2-1im) (-2+0im);(-1-1im) (-1+0im);(0-1im) (0+0im);(1-1im) (1+0im);(2-1im) (2+0im)]

    @test cpx(randn((5,2))) == Complex{Int8}[(-2-1im) (-2+0im);(-1-1im) (-1+0im);(0-1im) (0+0im);(1-1im) (1+0im);(2-1im) (2+0im)]
    @test cpx((5,2)) == Complex{Int8}[(-2-1im) (-2+0im);(-1-1im) (-1+0im);(0-1im) (0+0im);(1-1im) (1+0im);(2-1im) (2+0im)]
end

@testset "Test exp_ikx method" begin
    @test exp_ikx(ComplexF32, (5,)) ≈ ComplexF32[0.3090170656104166 - 0.9510564931493436im, -0.8090170163879177 + 0.5877852219942177im, 1.0 - 0.0im, -0.8090170163879177 - 0.5877852219942177im, 0.3090170656104166 + 0.9510564931493436im] 
    @test exp_ikx(randn((2,2))) ≈ [1.0 (-1.0);-1.0 1.0]
    @test exp_ikx(((2,2))) ≈ [1.0 (-1.0);-1.0 1.0]
    @test eltype(exp_ikx(ComplexF32, ((2,2)))) == ComplexF32 
    @test exp_ikx(ComplexF32, ((2,2)))[1] isa ComplexF32 
end

@testset "Check dims argument" begin
    @test abs.(rr((5,5), dims=1)) == abs.(xx((5,5)))
    @test abs.(rr((6,6), dims=1)) == abs.(xx((6,6)))

    @test idx((3, 3), dims = 2) == [(-0.0, -1.0) (-0.0, 0.0) (-0.0, 1.0); (0.0, -1.0) (0.0, 0.0) (0.0, 1.0); (0.0, -1.0) (0.0, 0.0) (0.0, 1.0)]
    @test idx((4, 4), dims = (1, 2)) == [(-2.0, -2.0) (-2.0, -1.0) (-2.0, 0.0) (-2.0, 1.0); (-1.0, -2.0) (-1.0, -1.0) (-1.0, 0.0) (-1.0, 1.0); (0.0, -2.0) (0.0, -1.0) (0.0, 0.0) (0.0, 1.0); (1.0, -2.0) (1.0, -1.0) (1.0, 0.0) (1.0, 1.0)]
end



@testset "Test window functions" begin

    # simple window_hanning test
    @testset "window_hanning window test" begin
        window_hanning((10, ), border_in=0, border_out=1) ≈ sinpi.(range(0, 1, length=11)[1:end-1]).^2
        window_hanning((10, ), border_in=0.4, border_out=0.499) ≈ [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
        
        x = window_hanning(Int, (10, ), border_in=0.4, border_out=0.499) 
        @test x ≈ [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
        @test typeof(x) <: IndexFunArray{Int, 1, T} where T
        
        x = window_hanning(ComplexF32, (10, ), border_in=0.4, border_out=0.499) 
        @test x ≈ [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
        @test typeof(x) <: IndexFunArray{ComplexF32, 1, T} where T
        

        x = window_hanning(Float32.(randn((10, ))), border_in=0.4, border_out=0.499) 
        @test x ≈ [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
        @test typeof(x) <: IndexFunArray{Float32, 1, T} where T
        
    end

    @testset "window_linear and window_radial_linear" begin

        @test window_linear((4, 4), border_in = 0, border_out = 1) == [0.0 0.0 0.0 0.0; 0.0 0.25 0.5 0.25; 0.0 0.5 1.0 0.5; 0.0 0.25 0.5 0.25]
        @test window_linear((4, 4), border_in = 0, border_out = 1, offset = CtrCorner, scale = ScaNorm) == [1.0 0.6666666666666667 0.33333333333333337 0.0; 0.6666666666666667 0.44444444444444453 0.22222222222222227 0.0; 0.33333333333333337 0.22222222222222227 0.11111111111111113 0.0; 0.0 0.0 0.0 0.0]

        @test window_radial_linear((4, 4), border_in = 0, border_out = 1, offset = CtrCorner, scale = ScaNorm) == [1.0 0.6666666666666667 0.33333333333333337 0.0; 0.6666666666666667 0.5285954792089683 0.2546440075000701 0.0; 0.33333333333333337 0.2546440075000701 0.057190958417936644 0.0; 0.0 0.0 0.0 0.0]
    end

    @testset "window_half_cos" begin
        @test sum(abs.(window_half_cos((1, 14), border_in = 0, border_out = 1) .- [0.0  0.222521  0.433884  0.62349  0.781831  0.900969  0.974928  1.0  0.974928  0.900969  0.781831  0.62349  0.433884  0.222521])) < 1e-5
    end

    @testset "box and disc" begin
        sz = (77,202)
        boxsize = (14,99)
        sc = (1.3,2.3)
        @test box(sz,boxsize) == (abs.(xx(sz) .+ 0.25) .< (boxsize[1]./2)) .* (abs.(yy(sz) .+ 0.25) .< (boxsize[2]./2))
        @test box(zeros(sz),boxsize,scale=sc) == (abs.(xx(sz,scale=sc) .+ 0.25) .< (boxsize[1]./2)) .* (abs.(yy(sz,scale=sc) .+ 0.25) .< (boxsize[2]./2))
        @test box(sz,boxsize, offset=CtrCorner) == (abs.(xx(sz, offset=CtrCorner) .+ 0.25) .< (boxsize[1]./2)) .* (abs.(yy(sz, offset=CtrCorner) .+ 0.25) .< (boxsize[2]./2))
        @test box(zeros(sz),boxsize, offset=CtrCorner) == (abs.(xx(sz, offset=CtrCorner) .+ 0.25) .< (boxsize[1]./2)) .* (abs.(yy(sz, offset=CtrCorner) .+ 0.25) .< (boxsize[2]./2))
        @test box((10,10)) == box((10,10),(5,5))
        @test sum(box((10,12))) == 5*6
        @test sum(box1((10,12), scale=0.3)) == 49
        disc_rad = 33.0
        @test disc(sz,disc_rad) ≈ (rr(sz) .<= disc_rad)
        @test disc(zeros(sz),disc_rad,scale=sc) ≈ (rr(sz, scale=sc) .<= disc_rad)
        @test disc(sz,disc_rad, offset=CtrCorner) ≈ (rr(sz, offset=CtrCorner) .<= disc_rad)
        @test disc(zeros(sz),disc_rad, offset=CtrCorner) ≈ (rr(sz, offset=CtrCorner) .<= disc_rad)
        @test sum(disc((10,12))) == 89
        @test sum(disc1((10,11), scale=0.3)) == 37
    end
    @testset "axes1d" begin
        sz = (2,3,4)
        x,y,z = axes1d(sz, offset=(0,0,0))
        @test x == 1:sz[1]
        @test y == 1:sz[2]
        @test z == 1:sz[3]
        x,y,z = axes1d(sz, offset=CtrCorner, keepdims=true)
        @test size(x) == (sz[1],)
        @test size(y) == (1,sz[2])
        @test size(z) == (1,1,sz[3])
        @test x[:] == 0:sz[1]-1
        @test y[:] == 0:sz[2]-1
        @test z[:] == 0:sz[3]-1
    end
end



