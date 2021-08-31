using Test
using Smoothers
import Smoothers: filter

@testset "filter" begin

    # Moving Average Filter
    w = 3; 
    b = ones(w)/w;
    a = [1];
    x = [1,2,3,4,5];

    fx = filter(b,a,x)
    @test fx[1:2] ≈ [1/3,1]
    @test fx[3:5] ≈ sma(x,3)

    si = [1/3,1]
    fx = filter(b,a,x,si)
    @test fx[1:2] ≈ [2/3,2]
    @test fx[3:5] ≈ sma(x,3)
    
end
