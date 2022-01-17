# Provides dynamic risk measures for SDEs


function dynamicRM(prob2::SDEProblem,eval = x-> sum(x), RM=1.0, dt = 0.01, MCs = 100, iter = 100
                    u0::Flux.Chain=Flux.Chain() )


    time_steps = div(prob2.tspan[2]-prob2.tspan[1],dt, RoundNearest)


    f(X,u,σᵀ∇u,p,t) = RM*abs.(σᵀ∇u)[1]                 # nonlinear riskMeasure part



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

    σᵀ∇u = [Flux.Chain(Dense(d,hls,relu),
        Dense(hls,hls,relu),
        Dense(hls,d)) for i in 1:time_steps]




    pdealg = NNPDEHan(u0,σᵀ∇u;opt=Flux.ADAM(0.1))
    ans = solve(prob, pdealg, verbose=true, maxiters=iter, trajectories=MCs, dt=dt)

    return(ans)
end


# Provide a jump diffusion version of dynamicRM