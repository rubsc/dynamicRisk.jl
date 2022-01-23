module dynamicRisk


using Reexport
@reexport using DifferentialEquations
@reexport using ScenTreesMakie
@reexport using riskMeasures

using LinearAlgebra, NeuralPDE, Flux
using Test


include("BSDE.jl")


export dynamicRM


end # module
