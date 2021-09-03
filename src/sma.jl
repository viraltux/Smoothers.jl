"""
Package: Smoothers

    sma(x, n)
    sma(x, n, center)

Smooth a vector of data using a simple moving average

# Arguments
- `x`: Vector of data.
- `n`: Size of the moving average.
- `center`: returns a vector of the same size of x with missing values on the tails.

# Returns
Vector of moving average values

# Examples
```julia-repl
julia> sma(1:5,3)
3-element Vector{Float64}:
 2.0
 3.0
 4.0

julia> sma(1:5,3,true)
5-element Vector{Union{Missing, Float64}}:
  missing
 2.0
 3.0
 4.0
  missing
```
"""
@inline function sma(x::AbstractVector{T}, n::Integer, center::Bool=false) where T<:Real

    n == 1 && return x
    N = length(x)
    @assert 1 <= n <= N
    
    V = Base.promote_op(/, T, T)
    res = Vector{V}(undef, N-n+1)
    
    # initial moving average value
    res[1] = ma = sum(x[1:n])/n
    for i in 1:N-n
        @inbounds res[1+i] = ma += (x[n+i] - x[i]) / T(n)
    end

    center ? vcat(repeat([missing], n÷2),res,repeat([missing], n-n÷2-1)) : res

end
