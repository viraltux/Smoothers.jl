module Smoothers

using Dierckx
import Base: filter

export filter, hma, loess, sma, stl

# Methods
include("filter.jl")
include("hma.jl")
include("loess.jl")
include("sma.jl")
include("stl.jl")

"""
Collection of smoothing heuristics, models and smoothing related applications. The current available smoothers and applications are:

    `filter`: Linear Time-invariant Difference Equation Filter (Matlab/Octave)
    `hma`:    Henderson Moving Average Filter
    `loess`:  Locally Estimated Scatterplot Smoothing
    `sma`:    Simple Moving Average
    `stl`:    Seasonal and Trend decomposition based on Loess

"""
Smoothers

end
