"""
Package: Smoothers

    filter(b,a,x,si)

Apply a digital filter to x using the following linear, time-invariant difference equation

```math
y(n) = \\sum_{k=0}^M d_{k+1} \\cdot x_{n-k}-\\sum_{k=1}^N c_{k+1} \\cdot y_{n-k} 
\\\\ \\forall n | 1\\le n \\le \\left\\| x \\right\\| | c=a/a_1, d=b/a_1
```

The implementation follows the description of [Octave filter function](https://octave.sourceforge.io/octave/function/filter.html)

# Arguments
- `a`: Vector of numerator coefficients for the filter rational transfer function.
- `b`: Vector of denominator coefficients for the filter rational transfer function.
- `x`: Vector of data.
- `si`: Vector of initial states.

# Returns

Vector of filtered data

# Examples
```julia-repl
using Plots

t = Array(LinRange(-pi,pi,100));
x = sin.(t) .+ 0.25*rand(length(t));

# Moving Average Filter
w = 5; 
b = ones(w)/w;
a = [1];

plot(t,x,label="sin(x)",legend=:bottomright)
y1 = filter(b,a,x)
si = x[1:4] .+ .1;
y2 = filter(b,a,x,si)
plot!(t,y1,label="MA")
plot!(t,y2,label="MA with si")
```
"""
@inline function filter(b::AbstractVector{A},
                        a::AbstractVector{B},
                        x::AbstractVector{C},
                        si::AbstractVector{D}=zeros(C,max(length(a),length(b))-1)
                        ) where {A<:Real,B<:Real,C<:Real,D<:Real}

    @assert a[1] != 0 "a[1] must not be zero"

    T = Base.promote_op(/,B,A)
    a,b,x,si,_ = Base.promote(a,b,x,si,[T(1.0)])
    
    Na,Nb,Nx = length(a),length(b),length(x)
    Nsi = max(Na,Nb)-1
    @assert Nsi == length(si) "length(si) must be max(length(a),length(b))-1)"
    
    N,M = Na-1,Nb-1
    c,d = a/a[1],b/a[1]
    
    y = zeros(T,Nx)
    y[1:Nsi] = si

    for n in 1:Nx
        for k in 0:min(n-1,M)
            @inbounds y[n] += d[k+1]*x[n-k]
        end
        for k in 1:min(n-1,N)
            @inbounds y[n] -= c[k+1]*y[n-k]
        end
    end
    
    return y
end
