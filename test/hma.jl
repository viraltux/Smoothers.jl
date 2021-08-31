using Test
using Smoothers
import Smoothers: hmaSymmetricWeights

@testset "hma" begin

    # weights
    w = hmaSymmetricWeights(8)
    @test ismissing(w[end])

    w = hmaSymmetricWeights(7)
    @test sum(w) ≈ 1.0

    w = hmaSymmetricWeights(13)
    @test sum(w) ≈ 1.0
    # (https://www.mathworks.com/help/econ/seasonal-adjustment-using-snxd7m-seasonal-filters.html)
    @test length(w) == 13
    @test round.(w,digits=3) ≈ [-0.019, -0.028, 0.0, 0.065, 0.147, 0.214, 0.24,
                                 0.214, 0.147, 0.065, 0.0, -0.028, -0.019]  

    # 13-term moving average
    x = hma(sin.(1:100), 13)
    @test length(x) == 100
    @test first(x)  ≈ 0.5636212576875559
    @test last(x)   ≈ -0.6658500507547161
    @test x[50]     ≈ -0.043655947941143455

end
