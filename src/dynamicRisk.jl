module dynamicRisk

using Reexport
@reexport using DifferentialEquations
@reexport using ScenTreesMakie

import riskMeasures: VaR, AVaR, EVaR, EVaR2, CTE

using LinearAlgebra, NeuralPDE, Flux
using Test


include("BSDE.jl")
include("tree.jl")
include("lattice.jl")


export dynamicRM, Var, CTE, AVaR, EVaR, EVaR2


end # module