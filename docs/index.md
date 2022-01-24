```@meta
CurrentModule = riskMeasures
```

# riskMeasures.jl

We present `riskMeasures.jl` package for evaluating various risk measures on samples as well as distribution objects from `Distribution.jl`. The main focus of this package is on static coherent risk measures but some other statistics are provided for completeness. 
The companion package `dynamicRisk.jl` builds on this package to implement (discrete time) dynamic risk measures for general trees and lattices as well as dynamic risk measures in continuous time for SDEs and certain Lévy processes.

## Main features of the package

1. We provide several implementations of well-known risk measures for discrete random variables in a static setting. 

2. We also provide a general template for coherent and convex risk measures based on the dual formulation from convex analysis. The user must only specify the set of feasible densities and in the convex case a convex function ``\\rho^*`` with the specified domain. ``\\rho^*`` is then the convex conjugate of specified convex risk measures.

To understand the general template it is recommended to have some deeper understanding of convex analysis.

## Installation

The package `riskMeasures.jl` can be installed in Julia REPL as follows:

```julia
julia> using Pkg
julia> Pkg.add("https://github.com/rubsc/riskMeasures.jl")
julia> using riskMeasures
```
