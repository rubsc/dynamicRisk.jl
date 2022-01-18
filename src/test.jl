using Flux, LinearAlgebra

using NeuralPDE

using Random
Random.seed!(100)

# one-dimensional heat equation
x0 = [11.0f0]  # initial points
tspan = (0.0f0,5.0f0)
dt = 0.01   # time step
time_steps = div(tspan[2]-tspan[1],dt, RoundNearest)
d = 1      # number of dimensions
m = 10     # number of trajectories (batch size)

g(X) = sum(X.^2)   # terminal condition
f(X,u,σᵀ∇u,p,t) = 0.0  # function from solved equation
μ_f(X,p,t) = 0.0
σ_f(X,p,t) = 1.0
prob = TerminalPDEProblem(g, f, μ_f, σ_f, x0, tspan)

hls = 10 + d #hidden layer size
opt = Flux.ADAM(0.005)  #optimizer
#sub-neural network approximating solutions at the desired point
u0 = Flux.Chain(Dense(d,hls,relu),
                Dense(hls,hls,relu),
                Dense(hls,1))
# sub-neural network approximating the spatial gradients at time point
σᵀ∇u = [Flux.Chain(Dense(d,hls,relu),
                  Dense(hls,hls,relu),
                  Dense(hls,d)) for i in 1:time_steps]

alg = NNPDEHan(u0, σᵀ∇u, opt = opt)

ans = solve(prob, alg, verbose = true, abstol=1e-8, maxiters = 200, dt=dt, trajectories=m)