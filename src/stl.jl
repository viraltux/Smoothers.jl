"""
Package: Smoothers

    stl(Yv, np; robust=false, 
                nl=nextodd(np), 
                ns=10*length(Yv)+1,
                nt=nextodd(1.5*np/(1-1.5/ns)), 
                ni=robust ? 1 : 2,
                no=0,
                spm=false,
                qsmp=max(div(np,7),2),
                cth = 0.01,
                verbose=false)

Decompose a time series into trend, seasonal, and remainder components.

"STL has a simple design that consists of a sequence of applications of the loess smoother; the simplicity allows analysis of the properties of the procedure and allows fast computation, even for very long time series and large amounts of trend and seasonal smoothing. Other features of STL  are specification of amounts of seasonal and trend smoothing that range, in a nearly continous way, from very small amount of smoothing to a very large amount; robust estimates of the trend and seasonal components that are not distorted by aberrant behavior in the data; specification of the period of the seasonal component to any intenger multiple of the time sampling interval greater than one; and the ability to decompose time series with missing values."*

All default values are chosen following the recommendations of the original paper when those were recommended. `ns` is recommended to be chosen of the basis of knowledge of the time series and on the basis of diagnostic methods; it must nonethelessbe  always odd and at least 7.

for `no` the authors advise 5 ("safe value") or 10 ("near certainty of convergence") cycles  or a convergence criterion when robustness is required, in this case when `robust` is true computations stop when convergence is achieved in trend and seasonality.

for `qsmp` the authors do not advise a default but they use a value close to div(`np`,7).

# Arguments
- `Yv`: Time series.
- `np`: Seasonality.
- `robust`: Robust stl.
- `nl`: Smoothing parameter of the low-pass filter.
- `ns`: Smoothing parameter for the seasonal component.
- `nt`: Smoothing parameter for the trend decomposition.
- `ni`: Number of inner loop cycles.
- `no`: Number of outer loop cycles.
- `spm`: Seasonal post-smoothing.
- `qsmp`: Loess q window for Seasonal post-smoothing.
- `cth`: Corvengence threshold for Seasonal and Trend.
- `verbose`: If true shows updates for the Seasonal and Trend convergence.

# Returns
A three columns matrix with the Seasonal, Trend and Remainder values for Yv.

* STL: A Seasonal, Trend Decomposition Procedure Based on Loess" Robert B. Cleveland, William S. Cleveland, Jean E. McRae, and Irma Terpenning. Journal of Official Statistics Vol. 6. No. 1, 1990, pp. 3-73 (c) Statistics Sweden.

            
# Examples
```julia-repl
x = stl(rand(100),10)
100×3 Matrix{Float64}:
  0.0659717   0.512964  -0.3029
 -0.0822641   0.502129   0.0710054
  0.217383    0.491153   0.145449
[...]
```
"""
function stl(Yv::AbstractVector{<:Union{Missing,T}},
             np::Integer;
             robust=true,
             nl=nextodd(np),
             ns=7,
             nt=nextodd(1.5*np/(1-1.5/ns)),
             ni=robust ? 1 : 2,
             no=0,
             spm=true,
             qsmp=max(div(np,7),2),
             cth = 0.01,
             verbose=false) where {T<:Real}

    @assert mod(ns,2)==1 & (ns>=7) "`ns` is chosen on the knowledge of the time series and diagnostic methods; must always be odd and at least 7"

    function B(u)
        (u < 1) & !ismissing(u) ? (1.0-u^2)^2 : 0.0
    end

    N = length(Yv)
    # initial robustness weights
    rhov = ones(T,N)
    # intial trend
    Tv = Tv0 = zeros(T,N)
    Sv = Sv0 = zeros(Union{Missing,T},N)
    Rv = Vector{T}(undef,N)
    Cv = Vector{T}(undef,N+2*np)
    scnv = false # seasonal convergence flag
    tcnv = false # trend convergence flag
    #for o in 0:no
    o = 0
    while robust | (o <= no)
        for k in 1:ni
            # Updating sesonal and trend components
            ## 1. Detrending (Yv = Tv + Sv)
            Sv = Yv - Tv

            ### Seasonal convergence criterion
            Md = maximum(abs.(skipmissing(Sv-Sv0)))
            M0 = maximum(skipmissing(Sv0))
            m0 = minimum(skipmissing(Sv0))
            scnv = (Md/(M0-m0) < cth)
            if verbose
                println("Outer loop: " * string(o) * " - " * "Inner loop: " * string(k))
                println("Seasonal Convergence: " * string(Md/(M0-m0)))
            end
            Sv0 = Sv

            # Seasonal Smoothing 2,3,4
            ## 2. Cycle-subseries Smoothing
            for csi in 1:np
                srange = csi:np:N
                floess = loess(srange,Sv[srange];q=ns,d=1,rho=rhov[srange])
                Cv[csi:np:N+2*np] = floess(csi-np:np:N+np)

            end
            ## 3. Low-Pass Filtering of Smoothed Cycle-Subseries
            ### centered support instead 1:N to balance out machine error
            floess = loess((1:N).-N÷2,sma(sma(sma(Cv,np),np),3),d=1,q=nl,rho=rhov)
            Lv = floess((1:N).-N÷2)                  

            ## 4. Detreending of Smoothed Cycle-Subseries
            ### Lv is substracted to prevent low-frenquency power
            ### from entering the seasonal component.
            Sv = Cv[np+1:end-np] - Lv

            ## 5. Deseasonalizing
            Dv = Yv - Sv

            # Trend Smoothing
            ## 6. Trend Smoothing
            ### centered support instead 1:N to balance out machine error
            ### (floor isntead ceil like in Lv to balance out even lengths)
            floess = loess((1:N).-N÷2,Dv,q=nt,d=1,rho=rhov)
            Tv = floess((1:N).-N÷2)

            ### Trend convergence criterion
            Md = maximum(abs.(skipmissing(Tv-Tv0)))
            M0 = maximum(skipmissing(Tv0))
            m0 = minimum(skipmissing(Tv0))
            tcnv = (Md/(M0-m0) < cth)
            if verbose
                println("Trend    Convergence: " * string(Md/(M0-m0)) * "\n")
            end
            Tv0 = Tv

            scnv & tcnv ? break : nothing
        end
        # Computation of robustness weights
        ## These new weights will reduce the influence of transient behaviour
        ## on the trend and seasonal components.

        # Tv and Sv defined everywhere but Rv is not defined
        # where Yv has missing values

        Rv = Yv - Tv - Sv

        if scnv & tcnv
            @info "Corvengence achieved (< " * string(cth) * "); Stopping computation..."    
            break 
        end
        
        if 0 < o <= no
            smRv = skipmissing(Rv)
            h = 6*median(abs.(smRv))
            rhov = B.(abs.(Rv)/h)
        end
        o += 1
    end

    if spm
        # Using div(np,7) as default to approximate
        # the q=51 for a np=365 chosen in the original paper
        floess = loess((1:N).-N÷2,Sv; d=2, q=qsmp, exact = (1:N).-N÷2)
        Sv = floess((1:N).-N÷2)
        Rv = Yv - Tv - Sv
    end
    
    if !scnv
        @warn "Seasonal convergence not achieved (>= " * string(cth) * ")
         Consider a robust estimation"
    end

    if !tcnv
        @warn "Trend convergence not achieved (>= " * string(cth) * ")
         Consider a robust estimation"
    end

    [Sv Tv Rv]
end

"""
Package: Smoothers

    nextodd(x)

Return the smallest odd integer greater than or equal to `x`.        
"""
function nextodd(x::Real)::Integer
    cx = Integer(ceil(x))
    mod(cx,2) == 0 ? cx+1 : cx
end
