using Test
using Smoothers

@testset "sma" begin
    
    @test_throws AssertionError sma(1:5,10)
    
    x = sma(1:5,3)
    @test length(x) == 3
    @test x == [2.0, 3.0, 4.0]

    x = sma(1:5,3,true)
    @test length(x) == 5
    @test x == [1.0, 2.0, 3.0, 4.0, 5.0]
    

    x = sma(1:5,3,false)
    @test length(x) == 3
    @test x == [2.0, 3.0, 4.0]

end
