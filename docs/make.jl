using Documenter
using dynamicRisk

push!(LOAD_PATH,"../src/")
makedocs(sitename="dynamicRisk.jl Documentation",
         pages = [
            "Home" => "index.md",
            "Mathematical Background" =>Any[
                "back/DiscreteTime.md",
                "back/Limit.md",
                "back/Continuous.md"
            ],
            "Tutorials" => Any["tut/basic.md",
				    "tut/Tree.md",
				    "tut/Lattice.md",
                    "tut/BSDE.md"
            ],
            "Function Library" => "library.md",
         ],
         format = Documenter.HTML(prettyurls = false)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/rubsc/dynamicRisk.jl.git",
    devbranch = "main"
)