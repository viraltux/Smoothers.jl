using Test
using Smoothers

@testset "stl" begin

    @test_throws AssertionError stl(rand(100),10; ns=5)
    
    x = stl(rand(100),10)
    @test x isa Matrix
    @test size(x) == (100,3)
    
end
