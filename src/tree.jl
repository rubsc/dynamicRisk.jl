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



#function dynamicRM(tree0::lattice,eval = x-> sum(x), RM=1.0 )

    # recursively go through lattice and apply RM
    
    
#    return(nothing)
#end



""" 
    VaR(states,prob,alpha)
implements the Value-at-Risk at level ``\\alpha`` defined by
```math
VaR_\\alpha (Y) = \\arg \\min_x \\left( x\\in \\mathbb{R} : F_Y(x) \\geq \\alpha \\right),
```
for the random variable ``Y`` defined by `states` and `prob`.
"""
function VaR(trr::Tree,alpha::Float64)
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
    CTE(states,prob,alpha)
implements the Conditional Value-at-Risk at level ``\\alpha`` defined by
```math
CTE_\\alpha (Y) = VaR_\\alpha(Y) + \\frac{1}{1-\\alpha} \\mathbb{E} \\left( Y- VaR_\\alpha (Y) \\right)_+ ,
```
for the random variable ``Y`` defined by `states` and `prob`.
"""
function CTE(trr::Tree,alpha::Float64)
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
    EVaR2(states,prob,beta)
Solves the optimization problem associated with the primal formulation of the Entropic Value-at-Risk:
```math
EVaR_\\alpha(Y) = \\min_{x >0} \\frac{1}{x} \\left( \\beta +  \\log\\mathbb{E} e^{xY} \\right),
```
where ``Y`` is the discrete random variable defined by `states` and `prob`.
Here the optimization is done via the goldenSearch optimization routine implemented as part of this package. 
"""
function EVaR2(trr::Tree,beta::Float64)
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
    EVaR(states,prob,beta)
Solves the optimization problem associated with the primal formulation of the Entropic Value-at-Risk:
```math
EVaR_\\alpha(Y) = \\min_{x >0} \\frac{1}{x} \\left( \\beta +  \\log\\mathbb{E} e^{xY} \\right),
```
where ``Y`` is the discrete random variable defined by `states` and `prob`. Here, the optimization is done using JuMP and Ipopt.  
"""
function EVaR(trr::Tree,beta::Float64)
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
    AVaR(states,prob,alpha)
Solves the optimization problem associated with the primal formulation of the Average Value-at-Risk:
```math
AVaR_\\alpha(Y) = \\min_{x\\in \\mathbb{R}} x + \\frac{1}{1-\\alpha} \\mathbb{E} \\left( Y - x \\right)_+,
```
where ``Y`` is the discrete random variable defined by `states` and `prob`.
"""
function AVaR(trr::Tree, alpha::Float64)
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

