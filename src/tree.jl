"""
    dynamicRM(tree0::tree,eval = x-> sum(x), RM=1.0 ) 

Implements the nested version of RM for the tree specified by tree0.

"""
function dynamicRM(tree0::tree,eval = x-> sum(x), RM=1.0 )

# recursively go through tree and apply RM
Tree1 = deepcopy(Tree0);

# For every parent node going backwards
tmp = reverse(unique(Tree2.parent)); pop!(tmp)
for i in tmp
    #Note that the probs are the conditional ones
    #Risk measure w/ adjustment of values
    states2 = Tree2.state[Tree2.children[i+1]]
    prob2 = Tree2.probability[Tree2.children[i+1]]

    Tree2.state[i] = mSD(0.95,5,states2,prob2)
end


return(nothing)
end



#function dynamicRM(tree0::lattice,eval = x-> sum(x), RM=1.0 )

    # recursively go through lattice and apply RM
    
    
#    return(nothing)
#end