"""
    dynamicRM(trr::tree, RM::Array{Any,1} ) 

Implements the nested version of RM for the tree specified by trr.
The user has to specify a vector of risk measures of lenght N, where N is the number of stages of trr. 
Furthermore the risk measures RM[i] should be of the form RM[i](states, prob). 

"""
function dynamicRM(trr::Tree; RM::Array{Any,1})

    if !(length(RM) == maximum(stage(trr)) )
        println("Error, number of risk measures should be equal to the number of stages")
        return(nothing)
    end
# recursively go through tree and apply RM
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]
        T0.state[i] = RM[i](states,prob)
    end


    return(T0[1],T0)
end


""" 
    VaR(trr::Tree,alpha::Float32)
implements the nested version of the Value-at-Risk at level ``\\alpha`` for the Tree ``trr``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version.
"""
function VaR(trr::Tree,alpha::Float32)
    # look for smallest x such that  P(-states <= x) >alpha , i.e.
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = riskMeasures.VaR(states,prob,alpha)
    end
    return(T0.state[1],T0)
end

""" 
    CTE(trr::Tree,alpha::Float32)
implements the nested version of the Conditional Tail Expectation at level ``\\alpha`` for the Tree ``trr``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version.
"""
function CTE(trr::Tree,alpha::Float32)
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]
        println(states)
        println(prob)
        T0.state[i] = riskMeasures.CTE(states,prob,alpha)
        println("risk value:")
        println(T0.state[i])
    end
    return(T0.state[1],T0)
end



""" 
    EVaR2(trr::Tree,beta::Float32)
implements the nested version of the Entropic Value-at-Risk at level ``\\beta`` for the Tree ``trr``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version.
"""
function EVaR2(trr::Tree,beta::Float32)
	T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = riskMeasures.EVaR2(states,prob,beta)[1]
    end
    return(T0.state[1],T0)
end


""" 
    EVaR(trr::Tree,beta::Float32)
implements the nested version of the Entropic Value-at-Risk at level ``\\beta`` for the Tree ``trr``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version. This version uses JuMP and Ipopt.
"""
function EVaR(trr::Tree,beta::Float32)
	T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = riskMeasures.EVaR(states,prob,beta)[1]
    end
    return(T0.state[1],T0)
end



""" 
    AVaR(trr::Tree,beta::Float32)
implements the nested version of the Average Value-at-Risk at level ``\\alpha`` for the Tree ``trr``. 
See [here](https://rubsc.github.io/riskMeasures.jl/dev/) for the univariate version.
"""
function AVaR(trr::Tree, alpha::Float32)
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]
        println(states)
        println(prob)

        T0.state[i] = riskMeasures.AVaR(states,prob,alpha)[1]
    end
    return(T0.state[1],T0)
end

