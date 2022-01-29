

""" 
    VaR(lat::Lattice,alpha::Float32)
implements the nested version of the Value-at-Risk at level ``\\alpha`` for the lattice ``lat``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version. 
"""
function VaR(lat::Lattice,alpha::Float32)
    # look for smallest x such that  P(-states <= x) >alpha , i.e.
    T0 = deepcopy(lat);

    for i=(length(T0.state)-1):-1:1
		for j=1:length(T0.state[i])
			T0.state[i][j] = RiskMeasures.VaR(T0.state[i+1],T0.probability[i+1][j,:,1],alpha)
		end
	end
    return(T0.state[1],T0)
end

""" 
    CTE(lat::Lattice,alpha::Float32)
implements the nested version of the Conditional Tail Expectation level ``\\alpha`` for the lattice ``lat``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version. 
"""
function CTE(lat::Lattice,alpha::Float32)
    T0 = deepcopy(lat);

	for i=(length(T0.state)-1):-1:1
		for j=1:length(T0.state[i])
			T0.state[i][j] = RiskMeasures.CTE(T0.state[i+1],T0.probability[i+1][j,:,1],alpha)
		end
	end
    
    return(T0.state[1],T0)
end



"""
    EVaR2(lat::Lattice,beta::Float32)
implements the nested version of the Entropic Value-at-Risk at level ``\\beta`` for the lattice ``lat``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version. Here the optimization is done with a simple 
goldenSearch algorithm.
"""
function EVaR2(lat::Lattice,beta::Float32)
	T0 = deepcopy(lat);

	for i=(length(T0.state)-1):-1:1
		for j=1:length(T0.state[i])
			T0.state[i][j] = RiskMeasures.EVaR2(T0.state[i+1],T0.probability[i+1][j,:,1],beta)[1]
		end
	end

    return(T0.state[1],T0)
end


"""
    EVaR(lat::Lattice,beta::Float32)
implements the nested version of the Entropic Value-at-Risk at level ``\\beta`` for the lattice ``lat``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version. The optimization routine uses JuMP and Ipopt.
"""
function EVaR(lat::Lattice,beta::Float32)
	T0 = deepcopy(lat);

    for i=(length(T0.state)-1):-1:1
		for j=1:length(T0.state[i])
			T0.state[i][j] = RiskMeasures.EVaR(T0.state[i+1],T0.probability[i+1][j,:,1],beta)[1]
		end
	end
    return(T0.state[1],T0)
end



"""
    AVaR(lat::Lattice,alpha::Float32)
implements the nested version of the Average Value-at-Risk at level ``\\alpha`` for the lattice ``lat``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version. 
"""
function AVaR(lat::Lattice, alpha::Float32)
    T0 = deepcopy(lat);

    for i=(length(T0.state)-1):-1:1
		for j=1:length(T0.state[i])
			T0.state[i][j] = RiskMeasures.AVaR(T0.state[i+1],T0.probability[i+1][j,:,1],alpha)[1]
		end
	end
    return(T0.state[1],T0)
end