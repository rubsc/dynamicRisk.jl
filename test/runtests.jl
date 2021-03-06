using dynamicRisk
using Test


@testset "BSDE" begin
    # Write your tests here.
    X0 = [0.0f0];          # initial value of stochastic state
    tspan = (0.0f0,1.0f0);
    g(X) = sum(X);         # needs to be number not [x]                                # terminal condition
    μ_f(X,p,t) = 0.0f0; #Vector d x 1                    # drift of forward 
    σ_f(X,p,t) = 1.0f0; #Matrix d x d    # diffusion part of forward --> 0.5*σ^2 is the full factor

    prob2 = SDEProblem(μ_f,σ_f,X0,tspan);
    
    @test  dynamicRM(prob2) ≈ 1.0 atol=0.1
end


@testset "Tree" begin
    # Write your tests here.
    trr = Tree(402)
    @test AVaR(trr::Tree,alpha::Float64)[1] <= maximum(trr.state)
    

end
