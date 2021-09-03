# MIT License

# Copyright (c) 2021 Fran Urbano
# Copyright (c) 2021 Val Lyashov
# Copyright (c) 2020 Bryan Palmer

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Package: Forecast

    hma(s, n)

Applies the Henderson moving average filter to dataset `s` with `n`-term.

"Henderson moving averages are filters which were derived by Robert Henderson in 1916 
for use in actuarial applications. They are trend filters, commonly used in time series 
analysis to smooth seasonally adjusted estimates in order to generate a trend estimate.

They are used in preference to simpler moving averages because they can reproduce 
polynomials of up to degree 3, thereby capturing trend turning points.

The ABS uses Henderson moving averages to produce trend estimates from a seasonally 
adjusted series. The trend estimates published by the ABS are typically derived using 
a 13 term Henderson filter for monthly series, and a 7 term Henderson filter for quarterly series.

Henderson filters can be either symmetric or asymmetric. Symmetric moving averages can be applied 
at points which are sufficiently far away from the ends of a time series. In this case, the 
smoothed value for a given point in the time series is calculated from an equal number of values 
on either side of the data point." - Australian Bureau of Statistics (www.abs.gov.au)

# Arguments
- `s`: Observations' support.
- `n`: Observation values. Tipically 13 or 7 for quarterly data, larger values may need a BigInt type.

# Returns
An array of Henderson filter smoothed values provided in `s`.

# Examples
```julia-repl
julia> hma(rand(1000), BigInt(303)))
1000-element Vector{BigFloat}:
[...]
```
"""
@inline function hma(s::AbstractVector{<:T}, n::Integer) where T<:Real

    @assert isodd(n) "n must be odd"
    @assert n >= 5 "n must be greater or equal 5"
    @assert length(s) >= n "length(s) must be greater or equal to n"

    # Promote to equivalent in size T Integer type andd Float type
    n = typeof(Integer(T(1.0)))(n)
    s,_ = Base.promote(collect(s),[Base.promote_op(/,T,T)(1.0)])
    
    w = hmaSymmetricWeights(n)
    
    m = (n-1) ÷ 2
    ls = length(s)

    function hmai(i)
        if i - 1 < m
            u = hmaAsymmetricWeights(m + i, w)[end:-1:1]
            sum(s[1:i + m] .* u)
        elseif i - 1 + m >= ls
            u = hmaAsymmetricWeights(m + ls - i + 1, w)
            sum(s[i-m:ls] .* u)
        else
            sum(s[i-m:i+m] .* w)
        end
    end
    
    hmav = similar(s)
    for i in 1:ls
        @inbounds hmav[i] = hmai(i)
    end
    hmav
    
end

"""
package: Smoothers

    hmaSymmetricWeights(n,T)

Caluclate the hma symmetric weights for 'n'

# Arguments
- `n`: Number of symmetric weights

# Returns

 Vector of symmetrical hma weights with type related to the bit size of `n`

# Refenrences
- "A Guide to Interpreting Time Series" ABS (2003), page 41.
"""
@inline function hmaSymmetricWeights(n::T) where T<:Integer

    m = (n-1)÷2
    m1 = (m+1)^2
    m2 = (m+2)^2
    m3 = (m+3)^2
    d = 315/(8*(m+2)*(m2-1)*(4*m2-1)*(4*m2-9)*(4*m2-25))

    P = Base.promote_op(/,T,T)
    w = Vector{P}(undef,m+2)
    m,m1,m2,m3,d = promote(m,m1,m2,m3,d,P(1.0))

    for (i,x) in enumerate(0:m+1)
        @inbounds w[i] = (m1-x^2)*(m2-x^2)*(m3-x^2)*(P(3.0)*m2-P(11.0)*x^2-P(16.0))*d
    end
    u = vcat(w[end-1:-1:1],w[2:end-1])

    return mod(n, 2) != 0 ? u : vcat(u, missing)

end

"""
package: Smoothers

    hmaAsymmetricWeights(m,w)

Calculate the hma asymmetric end-weights for 'm' given 'w'

# Arguments
- `m`: Number of asymmetric weights (m < length(w))
- `w`: Vector of symmetrical hma weights

# Returns

Vector of asymmetrical hma weights

# References

- Mike Doherty (2001), "The Surrogate Henderson Filters in X-11",Aust, NZ J of Stat. 43(4), 2001, 385–392
"""
@inline function hmaAsymmetricWeights(m::Integer, w::AbstractVector{<:T}) where T<:Real

    n = length(w)

    @assert m <= n "The m argument must be less than w"
    @assert m >= (n-1)÷2 "The m argument must be greater than (n-1)/2"

    sumResidual = sum(w[range(m + 1, n, step = 1)]) / T(m)
    
    sumEnd = T(0.0)
    for i in m+1:n
        @inbounds sumEnd += (i-(m+T(1.0))/T(2.0))*w[i]
    end
    
    ic = n < 13 ? T(1.0) : (13 <= n < 15 ? T(3.5) : T(4.5))
    
    b2s2 = T(4.0)/T(pi)/ic^2
    denominator = T(1.0) + ((m*(m-T(1.0))*(m+T(1.0)) / T(12.0) ) * b2s2)

    aw = Vector{T}(undef,m)
    for r in 1:m
        numerator = (T(r) - (m+T(1.0)) / T(2.0)) * b2s2
        @inbounds aw[r] = w[r] + sumResidual +  numerator / denominator * sumEnd
    end
    aw
end
