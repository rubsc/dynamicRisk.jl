# Provides dynamic risk measures for SDEs and tree processes

# In the future also more general lattices will be dealt with


################
# First a test
using NeuralPDE
using Flux
using DifferentialEquations
using LinearAlgebra


d = 1                   # number of dimensions
X0 = [0.0f0]            # initial value of stochastic state
tspan = (0.0f0,1.0f0)
r = 0.05f0
sigma = 1.0f0
f(X,u,σᵀ∇u,p,t) = abs.(σᵀ∇u)[1] #0.0f0 *r * (u - sum(X.*σᵀ∇u))       # nonlinear part
g(X) = sum(X)                                         # terminal condition
μ_f(X,p,t) = zero(X) #Vector d x 1                    # drift of forward 
σ_f(X,p,t) = sigma #Matrix d x d    # diffusion part of forward --> 0.5*σ^2 is the full factor
prob = TerminalPDEProblem(g, f, μ_f, σ_f, X0, tspan);



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


ans = solve(prob, pdealg, verbose=true, maxiters=550, trajectories=1000,
                            alg=EM(), dt=0.2, pabstol = 1f-6)



function dynamicRisk(process, RM)



    return(nothing)
end