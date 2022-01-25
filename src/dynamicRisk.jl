module dynamicRisk

#using Pkg; Pkg.add(url="https://github.com/rubsc/riskMeasures.jl")
using Reexport
@reexport using DifferentialEquations
@reexport using ScenTreesMakie

import riskMeasures

using LinearAlgebra, NeuralPDE, Flux
using Test


include("BSDE.jl")
include("tree.jl")
include("lattice.jl")


export dynamicRM, Var, CTE, AVaR, EVaR, EVaR2, Expectation, entropic,
        mSD, meanVariance, meanDeviation, meanSemiVariance, meanSemiDevi


end # module
