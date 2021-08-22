"""
Package: Forecast

    loess(xv, yv;
          d = 0,
          q = 3*sum((!ismissing).(yv))÷4,,
          rho = fill(1.0,length(yv)),
          exact = [])

Smooth a vector of observations using locally weighted regressions.

Although loess can be used to smooth observations for any given number of independent variables, this implementation is univariate. The speed of loess can be greatly increased by using fast aproximations for the linear fitting calculations, however this implementation calculates allows as well for exact results.

The loess functionality and nomenclature follows the descriptions in:

"STL: A Seasonal, Trend Decomposition Procedure Based on Loess"
Robert B. Cleveland, William S. Cleveland, Jean E. McRae, and Irma Terpenning.
Journal of Official Statistics Vol. 6. No. 1, 1990, pp. 3-73 (c) Statistics Sweden.

# Arguments
- `xv`: Observations' support.
- `yv`: Observation values.
- `d`: Degree of the linear fit, it accepts values 0, 1 or 2, if 0 an estimation of `d` is calculated.
- `q`: As q increases loess becomes smoother, when q tends to infinity loess tends to an ordinary least square poynomial fit of degree `d`. It defaults to the rounding of 3/4 of xv's non-missing values length.
- `rho`: Weights expressing the reliability of the observations (e.g. if yi had variances σ^2*ki where ki where known, the rhoi could be 1/ki). It defaults to 1.
- `exact`: Vector containing support values where exact loess will be calculated, the least values the faster will be loess. It defaults to an emtpy array. For 10 values or less exact loess is guaranteed.

# Returns
A function returning the exact loess for the values contained in `exact` and either exact loess or fast approximations for the rest.


# Examples
```julia-repl
n = 1000
x = sort(rand(n)*2*pi);
y = sin.(x) + randn(n);
f = Smoothers.loess(x,y)
#f(π) ≈ 0.0
[...]
```
"""
loess

function loess(yv::AbstractVector{<:Union{Missing,S}};
               d::Integer=0, # 0 -> auto
               q::Integer=3*sum((!ismissing).(yv))÷4,
               rho::AbstractVector{<:Union{Missing,T}}=fill(1.0,length(yv)),
               exact::AbstractVector{S}=S[]) where {S<:Real, T<:Real}

    loess(S(1):length(yv),yv;d,q,rho,exact)
end

function loess(xv::AbstractVector{R},
               yv::AbstractVector{<:Union{Missing,S}};
               d::Integer=0, # 0 -> auto
               q::Integer=3*sum((!ismissing).(xv))÷4,
               rho::AbstractVector{<:Union{Missing,T}}=fill(1.0,length(xv)),
               exact::AbstractVector{R}=R[]) where {R<:Real, S<:Real, T<:Real}
        
    @assert d in 0:2 "Linear Regression must be of degree 1 or 2 or 0 (auto)"
    @assert length(xv) == length(yv)

    # Removing missing values
    myi = findall(x -> !ismissing(x),yv)
    xv = xv[myi]
    yv = yv[myi]
    rho = rho[myi]

    @assert 1<=q<=length(xv) "q should be smaller of equal to the number of non missing values in xv"
    
    # Promote to same type
    P = promote_type(R,S,T)
    xv = R == P ? xv :  P.(xv)
    yv =  P.(yv)
    rho = P.(rho)
    exact = P.(exact)
    
    ## Ax = b prediction
    A = hcat(xv,fill(P(1),length(xv)))
    b = yv

    ## Order Estimation
    k = khat(yv)
    d = (d == 0) ? min(k,2) : 2  #loess accepts 1 or 2
    d == 2 && (A = hcat(xv .^ 2, A))
    
    # Select collection
    n = max(sr(length(xv)),10)
    m,M = extrema(xv)
    gap = (M-m)/n
    exact = vcat(exact, xv[Int64.(round.(LinRange(1,length(xv),10)))])
    predict = sort(unique(vcat(collect(m:gap:M),M,exact)))
    predictx = vcat(m-P(10)*gap,m-gap,m-P(0.1)*gap,
                    predict,
                    M+P(0.1)*gap,M+gap,M+P(10)*gap)

    # Predict collection
    res = similar(predictx)
    for (i,xi) in enumerate(predictx)
        res[i] = ghat(xi;A,b,d,q,rho)
    end
    Spline1D(predictx, res; k=d, bc="extrapolate", s=0)
end

# loess estimation
function ghat(x::T;
              A::AbstractMatrix{T},
              b::AbstractVector{T},
              d::Integer=2,
              q::Integer,
              rho::AbstractVector{T}) where T<:Real
              
    xv = A[:,d]
    yv = b

    ## λ_q
    n = length(xv)
    xvx = @. abs(xv-x)
    qidx = sortperm(xvx)[1:q]
    qdist = abs(xv[last(qidx)]-x)*max(T(1),q/n)

    ## upsilon
    w = zeros(n)
    for wi in qidx
        aq = abs(xv[wi]-x)/qdist
        w[wi] = max((1-aq^3)^3,T(0))
    end

    A = @. A*(w*rho)
    b = @. b*(w*rho)

    lsq_x = A\b

    d == 1 ? [x,T(1)]'*lsq_x : [x^2,x,T(1)]'*lsq_x

end

#Sturges' Rule
sr(n) = Int(round(1+3.322*log(n)))

# absolute standard deviation
function astd(x::AbstractVector{<:Real})
    n = length(x)
    sum(abs.(x.-sum(x)/n))/n
end

# k order estimation for Splines
function khat(y)
    dy = diff(y); ady = astd(dy)
    ady < 1e-10 && return 1
    mi = 1
    for i in 1:3
        dy = diff(dy)
        nady = astd(dy); nady < 1e-10 && break
        mi += (ady > nady) ? 1 : break
        ady = nady
    end
    mi+1
end
