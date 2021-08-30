module Smoothers

using Dierckx

export hma, loess, sma, stl

# Methods
include("hma.jl")
include("loess.jl")
include("sma.jl")
include("stl.jl")

end
