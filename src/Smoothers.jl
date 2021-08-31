module Smoothers

using Dierckx

export filter, hma, loess, sma, stl

# Methods
include("filter.jl")
include("hma.jl")
include("loess.jl")
include("sma.jl")
include("stl.jl")

end
