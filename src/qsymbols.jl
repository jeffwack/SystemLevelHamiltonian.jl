function onesparams(H)
    params = get_numsymbols(H)
    rules = []
    for sym in params
        push!(rules,Pair(sym,1))
    end
    return Dict(rules)
end

function op_creator(op,space,idx,Ncutoff)
    aon = op.aon
    if space isa FockSpace
        if aon == idx
            if op isa Create
                return QuantumToolbox.create(Ncutoff)
            elseif op isa Destroy
                return QuantumToolbox.destroy(Ncutoff)
            else
                error("Unknown operator type $op")
            end
        else
            return QuantumToolbox.qeye(Ncutoff)
        end
    elseif space isa NLevelSpace && length(space.levels) == 2
        if aon == idx
            if op isa Transition
                if op.j == space.GS && op.i != space.GS
                    return QuantumToolbox.sigmap()
                elseif op.i == space.GS && op.j != space.GS
                    return QuantumToolbox.sigmam()
                elseif op.i != space.GS && op.j != space.GS
                    return QuantumToolbox.sigmaz()
                else 
                    error("Don't know how to classify this transition $op")
                end
            else
                error("Unknown operator type $op")
            end
        else
            return QuantumToolbox.qeye(2)
        end
    else
        error("Unknown Hilbert space $space")
    end
end

function standard_initial_state(H)
    states = []
    for space in H.dimensions.to
        if space.size == 2 
            push!(states,spin_state(1//2,-1//2)) 
        else
            push!(states, fock(space.size,0))
        end
    end
    return QuantumToolbox.tensor(states...)
end


function convert_to_QT(exprlist,Ncutoff,paramrules)
    oldops = collect(get_qsymbols(sum(exprlist)))#TODO: replace get qsymbols with function from QuantumCumulants

    if check_hilberts(sum(oldops)) #TODO: replace with compatible hilberts from QuantumCumulants
        hilb = first(oldops).hilbert
    else
        error("Expressions do not have homogenous Hilbert space")
    end

    newops = []

    if hasfield(typeof(hilb),:spaces)
        subspaces = hilb.spaces

        
        for oldop in oldops
            newsubops = [op_creator(oldop,subspaces[idx],idx,Ncutoff) for idx in eachindex(subspaces)]
            newop = QuantumToolbox.tensor(newsubops...)
            push!(newops,newop)
        end
    else
        for oldop in oldops
            newop = op_creator(oldop,hilb,1,Ncutoff) 
            push!(newops,newop)
        end
    end

    #now we need to substitute these operators into the symbolic expression of the Hamiltonian
    oprules = Dict([old => new for (old,new) in zip(oldops,newops)])
    #=
    println("Oprules")
    println(oprules)
    println("exprlist")
    println(exprlist)
    println("paramrules")
    println(paramrules)
    println("done")
    =#
    num_expr = [substitute(X,paramrules) for X in exprlist]
    new_expr = [substitute(X,oprules) for X in num_expr]

    return new_expr
end

#Wrote this with chatgpt 
#TODO: see if this can be replaced with QuantumCumulants' fundamental_operators()
#does not seem to be the case, fundamental_operators() works on a hilbert space
function get_qsymbols(expr)
    if Symbolics.istree(expr)
        return union(get_qsymbols.(Symbolics.arguments(expr))...)
    else
        return (expr isa Number) || (expr isa SymbolicUtils.Symbolic) ? Set() : Set([expr])
    end
end 

function get_numsymbols(expr)
    if Symbolics.istree(expr)
        return union(get_numsymbols.(Symbolics.arguments(expr))...)
    else
        return (expr isa SymbolicUtils.Symbolic) ? Set([expr]) : Set() 
    end
end

function check_hilberts(expr)
    operators = get_qsymbols(expr)
    hilberts = [op.hilbert for op in operators]
    return all(y->y==hilberts[1],hilberts)
end

function get_additive_terms(expr)
    """
    Extract additive terms from a quantum operator expression.
    
    Takes an expression containing quantum operators and returns a list of terms
    that contain no addition, only multiplication. Summing all returned terms
    results in the original expression.
    
    Args:
        expr: A symbolic expression containing quantum operators
        
    Returns:
        Vector: List of terms without addition operators
    """
    if Symbolics.istree(expr)
        op = Symbolics.operation(expr)
        args = Symbolics.arguments(expr)
        
        if op === (+)
            # If this is an addition, recursively get terms from each argument
            terms = []
            for arg in args
                append!(terms, get_additive_terms(arg))
            end
            return terms
        else
            # If this is not an addition (multiplication, function call, etc.), 
            # return the whole expression as a single term
            return [expr]
        end
    else
        # If not a tree (atom), return as single term
        return [expr]
    end
end