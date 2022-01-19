"""
    dynamicRM(tree0::tree,eval = x-> sum(x), RM=1.0 ) 

Implements the nested version of RM for the tree specified by tree0.

"""
function dynamicRM(trr,eval = x-> sum(x); RM=1.0 )

    println(typeof(RM))
#    if typeof(RM) == func

# recursively go through tree and apply RM
#T0 = deepcopy(trr);

# For every parent node going backwards
#tmp = reverse(unique(T0.parent)); pop!(tmp)
#for i in tmp
    #Note that the probs are the conditional ones
    #Risk measure w/ adjustment of values
#    states2 = T0.state[T0.children[i+1]]
#    prob2 = T0.probability[T0.children[i+1]]



#    T0.state[i] = mSD(states2,prob2,0.95,5)
#end


return(nothing)
end



#function dynamicRM(tree0::lattice,eval = x-> sum(x), RM=1.0 )

    # recursively go through lattice and apply RM
    
    
#    return(nothing)
#end