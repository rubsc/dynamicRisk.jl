module dynamicRisk

using Reexport
@reexport using DifferentialEquations
@reexport using ScenTreesMakie

import RiskMeasures

using LinearAlgebra, NeuralPDE, Flux
using Test


include("BSDE.jl")
include("tree.jl")
include("lattice.jl")


export dynamicRM, Var, CTE, AVaR, EVaR, EVaR2, Expectation, entropic,
        mSD, meanVariance, meanDeviation, meanSemiVariance, meanSemiDevi


end # module