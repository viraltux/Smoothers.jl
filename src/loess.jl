"""
Package: Forecast

    loess(xv, yv;
          d = 0,
          q = 3*sum((!ismissing).(yv))÷4,
          rho = fill(1.0,length(yv)),
          exact = [], extra = [])

    loess(yv;
          d = 0,
          q = 3*sum((!ismissing).(yv))÷4,
          rho = fill(1.0,length(yv)),
          exact = [], extra = [])

Return a funcion to smooth a vector of observations using locally weighted regressions.

Although loess can be used to smooth observations for any given number of independent variables, this implementation is univariate. The speed of loess can be greatly increased by using fast aproximations for the linear fitting calculations, however this implementation calculates allows as well for exact results.

The loess functionality and nomenclature follows the descriptions in:

"STL: A Seasonal, Trend Decomposition Procedure Based on Loess"
Robert B. Cleveland, William S. Cleveland, Jean E. McRae, and Irma Terpenning.
Journal of Official Statistics Vol. 6. No. 1, 1990, pp. 3-73 (c) Statistics Sweden.

# Arguments
- `xv`: Observations' support, if not provided 1:length(yv) is used instead.
- `yv`: Observation values.
- `d`: Degree of the linear fit, it accepts values 0, 1 or 2, if 0 an estimation of `d` is calculated.
- `q`: As q increases loess becomes smoother, when q tends to infinity loess tends to an ordinary least square poynomial fit of degree `d`. It defaults to the rounding of 3/4 of xv's non-missing values length.
- `rho`: Weights expressing the reliability of the observations (e.g. if yi had variances σ^2*ki where ki where known, the rhoi could be 1/ki). It defaults to 1.
- `exact`: Vector containing support values where exact loess will be calculated, the least values the faster will be loess. It defaults to an emtpy array. For series of 10 values or less exact loess is guaranteed.
- `extra`: Vector containing support values where, besideds those values chosen by the heuristic, exact loess will be calculated.

# Returns
A function returning the exact loess for the values contained in `exact` and either exact loess or fast approximations for the rest.


# Examples
```julia-repl
n = 1000
x = sort(rand(n)*2*pi);
y = sin.(x) + randn(n);
f = Smoothers.loess(x,y)
isapprox(f(pi),0;atol=0.1) #true
[...]
```
"""
loess

function loess(yv::AbstractVector{<:Union{Missing,T}};
               d::Integer=0, # 0 -> auto
               q::Integer=3*length(yv)÷4,
               rho::AbstractVector{<:Union{Missing,T}}=fill(T(1),length(yv)),
               exact::AbstractVector{T}=T[],
               extra::AbstractVector{T}=T[]) where {T<:Real,S<:Real}

    loess(T(1):length(yv),yv;d,q,rho,exact,extra)
end

function loess(xv::AbstractVector{R},
               yv::AbstractVector{<:Union{Missing,T}};
               d::Integer=0,
               q::Integer=3*length(yv)÷4,
               rho::AbstractVector{<:Union{Missing,T}}=fill(T(1.0),length(xv)),
               exact::AbstractVector{R}=R[],
               extra::AbstractVector{R}=R[]) where {R<:Real,T<:Real}

    @assert issorted(xv) "xv values must be sorted"
    @assert d in 0:2 "Linear Regression must be of degree 1 or 2 or 0 (auto)"
    @assert length(xv) == length(yv) "xv and yv must have the same length"

    # Removing missing values
    myi = findall(x -> !ismissing(x),yv)
    @assert length(myi) > 0 "at least one yv value must be non missing"
    q = Int(floor(q*length(myi)/length(yv)))
    xv = xv[myi]
    yv = Vector{T}(yv[myi])
    rho = rho[myi]

    length(xv) == 1 && return x -> repeat([yv[1]],length(x))
    length(xv) == 2 && begin d=(yv[2]-yv[1])/(xv[2]-xv[1]); return x -> @. d*x+(yv[1]-d*xv[1]) end
    
    # Promote to same Float type
    xv,yv,rho,exact,extra,_ = Base.promote(xv,yv,rho,collect(exact),collect(extra),
                                           [Base.promote_op(/,T,T)(1.0)])
    
    ## Ax = b prediction
    P = eltype(xv)
    A = hcat(xv,fill(P(1.0),length(xv)))
    b = yv

    ## Order Estimation
    d = (d == 0) ? min(khat(yv),2) : d  # loess accepts 1 or 2
    d == 2 && (A = hcat(xv .^ 2, A))
    
    # Select collection
    predictx = exact
    if isempty(exact)
        N = length(xv)
        n = 2*sr(N)*N÷q
        m,M = extrema(xv)
        gap = P((M-m)/n)
        #exact = vcat(exact, xv[Int64.(round.(LinRange(1,length(xv),20)))])
        predict = sort(unique(vcat(collect(m:gap:M),M,extra)))
        predictx = vcat(m-P(10.0)*gap,m-gap,
                        predict,
                        M+gap,M+P(10.0)*gap)
    else
        @assert length(exact) > d "length(exact) should be greater than d"
    end
    
    # Predict collection
    res = similar(predictx)
    for (i,xi) in enumerate(predictx)
        res[i] = ghat(xi;A,b,d,q,rho)
    end

    Spline1D(predictx, res; k=d, bc="extrapolate",s=0)

end

# loess estimation
function ghat(x::T;
              A::AbstractMatrix{T},
              b::AbstractVector{T},
              d::Integer,
              q::Integer,
              rho::AbstractVector{T}) where T<:Real
              
    xv = A[:,d]
    yv = b

    ## λ_q
    n = length(xv)
    xvx = @. abs(xv-x)
    qidx = sortperm(xvx)[1:min(q,n)]
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

    d == 1 ? sum([x,T(1)].*lsq_x) : sum([x^2,x,T(1)].*lsq_x)

end

#Sturges' Rule
sr(n) = Int(round(1+3.322*log(n)))

# absolute standard deviation
function astd(x::AbstractVector{<:Union{Missing,Real}})
    astd(collect(skipmissing(x)))
end

function astd(x::AbstractVector{<:Real})
    n = length(x)
    sum(abs.(x.-sum(x)/n))/n
end


# k order estimation for Splines
function khat(y::AbstractVector{<:Union{Missing,Real}})
    khat(collect(skipmissing(y)))
end

function khat(y::AbstractVector{<:Real})
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
