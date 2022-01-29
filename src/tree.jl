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

        T0.state[i] = RiskMeasures.VaR(states,prob,alpha)
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
        T0.state[i] = RiskMeasures.CTE(states,prob,alpha)
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

        T0.state[i] = RiskMeasures.EVaR2(states,prob,beta)[1]
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

        T0.state[i] = RiskMeasures.EVaR(states,prob,beta)[1]
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

        T0.state[i] = RiskMeasures.AVaR(states,prob,alpha)[1]
    end
    return(T0.state[1],T0)
end

#########################################
# From basicRM.jl
"""
    Expectation(states,prob)

Computes the expectation for the random variable ``Y`` with values given by `states` under
the probability measure ``Q`` given by `prob`:

```math
\\mathbb{E} Y = states \\cdot prob
```
"""
function Expectation(trr::Tree, alpha::Float32)
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = RiskMeasures.Expectation(states,prob)
    end
    return(T0.state[1],T0)
end



"""
    entropic(states, prob,theta::Float32)

implements the entropic risk measure defined as
```math
\\rho_\\theta(Y) = \\theta \\cdot \\log \\mathbb{E} e^{\\frac{Y}{\\theta}},
```
where ``\\theta`` is greater 0 and ``Y`` is the random variable defined by `states` and `prob`.
"""
function entropic(trr::Tree, theta::Float32)
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = RiskMeasures.entropic(states,prob,theta)
    end
    return(T0.state[1],T0)
end




#########################################################
#Mean Semi-deviation
#########################################################
""" 
    mSD(states,prob,beta::Float32,p::Float32)

implements the mean semi-deviation of order `p` which is a coherent risk measure defined by
```math
mSD_\\beta^p (Y) = \\mathbb{E} Y  + \\beta \\lvert \\left( Y - \\mathbb{E}Y \\right)_+ \\rvert_p,
```
for the random variable ``Y`` defined by `states` and `prob`.
"""
function mSD(trr::Tree, beta::Float32,p::Float32)
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = RiskMeasures. mSD(states,prob,beta::Float32,p::Float32)
    end
    return(T0.state[1],T0)
end


"""
    meanVariance(states, prob,c::Float32)

implements the mean Variance risk measure defined by
```math
\\rho_c(Y) := \\mathbb{E}Y + c \\cdot \\mathbb{E} \\left( Y- \\mathbb{E}Y)^2 \\right),
```
where ``c >0`` and ``Y`` is the random variable defined by `states` and `prob`.
"""
function meanVariance(trr::Tree, c::Float32)
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = RiskMeasures.meanVariance(states, prob,c::Float32)
    end
    return(T0.state[1],T0)
end


"""
    meanDeviation(states, prob,c::Float32,p::Float32)

implements the mean deviation risk measure of order ``p\\geq 1`` defined by
```math
\\rho_c^p(Y) := \\mathbb{E}Y + c \\cdot \\lVert \\left( Y- \\mathbb{E}Y \\right)^2 \\rVert_p,
```
where ``c >0`` and ``Y`` is the random variable defined by `states` and `prob`.
"""
function meanDeviation(trr::Tree, c::Float32,p::Float32)
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = RiskMeasures.meanDeviation(states, prob,c::Float32,p::Float32)
    end
    return(T0.state[1],T0)
end



"""
    meanSemiVariance(states, prob,c::Float32,t::Float32)

implements the mean upper-semi variance risk measure from a target ``t`` defined by
```math
\\rho_{c,t}(Y) = \\mathbb{E}Y + c \\cdot \\mathbb{E} \\left( Y - t \\right)^2_+ ,
```
where ``c >0`` and ``Y`` is the random variable defined by `states` and `prob`.
"""
function meanSemiVariance(trr::Tree, c::Float32,t::Float32)
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = RiskMeasures.meanSemiVariance(states, prob,c::Float32,t::Float32)
    end
    return(T0.state[1],T0)
end



"""
    meanSemiDevi(states, prob,c::Float32,target::Float32,p::Float32)

implements the mean upper-semi variance risk measure of order ``p \\geq 1`` from a target ``t`` defined by
```math
\\rho_{c,t}^p(Y) = \\mathbb{E}Y + c \\cdot \\lVert \\left( Y - t \\right)_+ \\rVert_p ,
```
where ``c >0`` and ``Y`` is the random variable defined by `states` and `prob`.
"""
function meanSemiDevi(trr::Tree, c::Float32,target::Float32,p::Float32)
    T0 = deepcopy(trr);

    # For every parent node going backwards
    tmp = reverse(unique(T0.parent)); pop!(tmp)
    for i in tmp
        #Note that the probs are the conditional ones
        #Risk measure w/ adjustment of values
        states = T0.state[T0.children[i+1]]
        prob = T0.probability[T0.children[i+1]]

        T0.state[i] = RiskMeasures.meanSemiDevi(states, prob,c::Float32,target::Float32,p::Float32)
    end
    return(T0.state[1],T0)
end