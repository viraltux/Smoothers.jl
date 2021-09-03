using Test
using Smoothers
import Smoothers: filter

@testset "filter" begin

    # Types' promotion
    fx = filter(Float16.(b),Float16.(a),Float16.(x))
    @test eltype(fx) == Float16
    fx = filter(Float32.(b),Float32.(a),Float32.(x))
    @test eltype(fx) == Float32
    fx = filter(BigFloat.(b),BigFloat.(a),BigFloat.(x))
    @test eltype(fx) == BigFloat
    fx = filter(Float16.(b),Int32.(a),Float64.(x))
    @test eltype(fx) == Float64

    # Moving Average Filter
    w = 3; 
    b = ones(w)/w;
    a = [1];
    x = [1,2,3,4,5];

    # Standard
    fx = filter(b,a,x)
    @test eltype(fx) == Float64
    @test fx[1:2] ≈ [1/3,1]
    @test fx[3:5] ≈ sma(x,3)


    # Initial State
    si = [1/3,1]
    fx = filter(b,a,x,si)
    @test fx[1:2] ≈ [2/3,2]
    @test fx[3:5] ≈ sma(x,3)
    
end
