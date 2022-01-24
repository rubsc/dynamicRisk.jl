# Provides dynamic risk measures for SDEs

"""
    dynamicRM(problem::SDEProblem,
                eval = x-> sum(x),
                RM=1.0,
                dt = 0.01,
                MCs = 200,
                iter = 100,
                u0::Flux.Chain=Flux.Chain() )
implements the unique dynamic risk measure for the SDE ``problem``associated to the nonlinear term
```math
s_\\rho(t,x) \\lVert \\sigma(t,x) \\nabla u(t,x) \\rVert. 
```
For more details on dynamic risk measures for SDEs, its relationship to BSDEs and parabolic PDEs as well as discrete time dynamic risk measures
see the Section [BSDE]().
"""
function dynamicRM(problem::SDEProblem,eval = x-> sum(x), RM=1.0, dt = 0.01, MCs = 200, iter = 100,
                    u0::Flux.Chain=Flux.Chain() )


    time_steps = div(problem.tspan[2]-problem.tspan[1],dt, RoundNearest)


    f(X,u,σᵀ∇u,p,t) = RM*abs.(σᵀ∇u)[1]                 # nonlinear riskMeasure part



    BSDE = TerminalPDEProblem(eval, f, problem.f, problem.g, problem.u0, problem.tspan);
    
    d = length(problem.u0)
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
    ans = solve(BSDE, pdealg, verbose=true, maxiters=iter, trajectories=MCs, dt=dt)

    return(ans)
end


# Provide a jump diffusion version of dynamicRM