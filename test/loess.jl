using Test
using Smoothers

@testset "loess" begin

    @test_throws AssertionError loess(rand(5),rand(5); d=3)
    
    fx = loess(sort(rand(10)),rand(10))
    x = fx(rand(20))
    @test length(x) == 20

    fx = loess(sort(rand(10)), rand(10); d=1)
    x = fx(rand(20))
    @test length(x) == 20

    xv = sort(sin.(collect(1:5.)))
    fx = loess(xv,sin.(collect(1:5.)); d=2)
    x = fx(xv)
    @test isapprox(x,[0.8414709848078966, 0.9095525828909523, 0.1420295878550277, -0.7350481509902661, -0.9589242746631387]; atol=.1)

    xv = sort(sin.(collect(1:5.)))
    fx = loess(xv,sin.(collect(1:5.)); d=1)
    x = fx(xv)
    @test isapprox(x,[0.8414709848078965, 0.9092974268256818, 0.14112000805986713, -0.7568024953079286, -0.9589242746631383]; atol=.1)

    fx = loess(1:100,rand(1:100,100))
    x = fx(Int64.(round.(100*rand(20))))
    @test length(x) == 20

    fx = loess(1:2:200, vcat(rand(80),repeat([missing],20)))
    x = fx(1.:2:40)
    @test length(x) == 20
    
end
