# Provides dynamic risk measures for SDEs and tree processes

# In the future also more general lattices will be dealt with


function dynamicRM(process::SDEProblem,eval = x-> sum(x), RM=0.0, u0::Flux.Chain=Flux.Chain(), σᵀ∇u::Flux.Chain=Flux.Chain())
    # get information on forward process from SDEProblem

    # get information on RM from class RM which riskMeasures introduces


    f(X,u,σᵀ∇u,p,t) = abs.(σᵀ∇u)[1]                 # nonlinear riskMeasure part
    prob = TerminalPDEProblem(eval, f, prob2.f, prob2.g, prob2.u0, prob2.tspan);
    
    d = length(prob2.u0)
    hls  = 10 + d #hide layer size
    opt = Flux.ADAM(0.001)

    if (u0 == Flux.Chain())
        u0 = Flux.Chain(Dense(d,hls,relu),
                Dense(hls,hls,relu),
                Dense(hls,hls,relu),
                Dense(hls,1))
    end
    if (σᵀ∇u == Flux.Chain())
        σᵀ∇u = Flux.Chain(Dense(d+1,hls,relu),
                  Dense(hls,hls,relu),
                  Dense(hls,hls,relu),
                  Dense(hls,hls,relu),
                  Dense(hls,d))
    end


    pdealg = NNPDENS(u0, σᵀ∇u, opt=opt)
    ans = solve(prob, pdealg, verbose=true, maxiters=50, trajectories=100,alg=EM(), dt=0.2, pabstol = 1f-6)
    return(ans)
end