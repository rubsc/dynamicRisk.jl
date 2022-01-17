# Provides dynamic risk measures for SDEs and tree processes

# In the future also more general lattices will be dealt with


################
# First a test




X0 = [0.0f0]            # initial value of stochastic state
tspan = (0.0f0,1.0f0)


g(X) = sum(X)                                         # terminal condition
μ_f(X,p,t) = zero(X) #Vector d x 1                    # drift of forward 
σ_f(X,p,t) = 1.0f0 #Matrix d x d    # diffusion part of forward --> 0.5*σ^2 is the full factor


process = SDEProblem(μ_f,σ_f,X0,tspan)

#Example:
# dynamicRisk(prob,eval=g, RM= ???)

#function dynamicRM(process::SDEProblem,eval, RM, u0::Flux.Chain=Flux.Chain(), σᵀ∇u::Flux.Chain=Flux.Chain())
    # get information on forward process from SDEProblem


    # get information on RM from class RM which riskMeasures introduces
    #eval(X) = sum(X)       # same as g above


    f(X,u,σᵀ∇u,p,t) = abs.(σᵀ∇u)[1]                 # nonlinear riskMeasure part
    prob = TerminalPDEProblem(eval, f, process.f, process.g, process.u0, process.tspan);
    
    d = length(process.u0)
    hls  = 10 + d #hide layer size
    opt = Flux.ADAM(0.001)

    u0 = Flux.Chain(Dense(d,hls,relu),
                Dense(hls,hls,relu),
                Dense(hls,hls,relu),
                Dense(hls,1))

    σᵀ∇u = Flux.Chain(Dense(d+1,hls,relu),
                  Dense(hls,hls,relu),
                  Dense(hls,hls,relu),
                  Dense(hls,hls,relu),
                  Dense(hls,d))


    pdealg = NNPDENS(u0, σᵀ∇u, opt=opt)
    ans = solve(prob, pdealg, verbose=true, maxiters=150, trajectories=100,alg=EM(), dt=0.2, pabstol = 1f-6)
    #return(nothing)
#end