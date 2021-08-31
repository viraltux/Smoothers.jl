"""
Package: Smoothers

    sma(x, n)
    sma(x, n, extend)

Smooth a vector of data using a simple moving average

# Arguments
- `x`: Vector of data.
- `n`: Size of the moving average.
- `extend`: if true extends the tails of the moving average using a linear spline.

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
5-element Vector{Float64}:
 1.0
 2.0
 3.0
 4.0
 5.0
```
"""
@inline function sma(x::AbstractVector{T}, n::Integer) where T<:Real

    n == 1 && return x
    N = length(x)
    @assert 1 <= n <= N
    
    V = Base.promote_op(/, T, typeof(n))
    res = Vector{V}(undef, N-n+1)
    
    # initial moving average value
    res[1] = ma = sum(x[1:n])/n
    for i in 1:N-n
        @inbounds res[1+i] = ma += (x[n+i] - x[i]) / n
    end

    res

end

function sma(x::AbstractVector{T}, n::Integer, extend::Bool) where T<:Real

    ma = sma(x,n)
    !extend && return ma

    n == 1 && return x    
    N = length(x)
    n = n-1
    tails = Spline1D(n÷2+1:N-(n-n÷2), ma; k=1, bc="extrapolate")
    vcat(n÷2>0 ? tails(1:n÷2) : [],ma,tails(N-(n-n÷2)+1:N))

end
