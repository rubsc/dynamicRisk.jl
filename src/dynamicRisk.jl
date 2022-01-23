module dynamicRisk

import Pkg; Pkg.add("https://github.com/rubsc/ScenTreesMakie.jl"); Pkg.add("https://github.com/rubsc/riskMeasures.jl")

using Reexport
@reexport using DifferentialEquations
@reexport using ScenTreesMakie
@reexport using riskMeasures

using LinearAlgebra, NeuralPDE, Flux
using Test


include("BSDE.jl")


export dynamicRM


end # module
