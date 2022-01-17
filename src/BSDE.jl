# Provides dynamic risk measures for SDEs and tree processes

# In the future also more general lattices will be dealt with


function dynamicRM(prob2::SDEProblem,eval = x-> sum(x), RM=0.0, u0::Flux.Chain=Flux.Chain(), σᵀ∇u::Flux.Chain=Flux.Chain())
    # get information on forward process from SDEProblem

    # get information on RM from class RM which riskMeasures introduces

    dt = 0.01   # time step
    time_steps = div(prob2.tspan[2]-prob2.tspan[1],dt, RoundNearest)


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
        σᵀ∇u = [Flux.Chain(Dense(d,hls,relu),
            Dense(hls,hls,relu),
            Dense(hls,d)) for i in 1:time_steps]
    end



    pdealg = NNPDEHan(u0,σᵀ∇u;opt=Flux.ADAM(0.1))
    ans = solve(prob, pdealg, verbose=true, maxiters=100, trajectories=200, dt=dt)
    return(ans)
end