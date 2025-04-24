module SystemLevelHamiltonian
#=
using GLMakie
include("QuantumMakie.jl")

#TODO update these in light of the nice access to solution by sol[variable]
export plotall, plotlist
=#
using QuantumCumulants
using Symbolics

using QuantumOptics

export SLH, get_qsymbols, promote, concatenation, feedbackreduce, QOHamiltonian, check_hilberts, get_numsymbols,onesparams

function onesparams(H)
    params = get_numsymbols(H)
    rules = []
    for sym in params
        push!(rules,Pair(sym,1))
    end
    return Dict(rules)
end

function op_selector(op,idx)
    aon = op.aon
    if idx == aon
        if op isa Create
            return create
        elseif op isa Destroy
            return destroy
        else
            error("Unknown operator type")
        end
    else
        return one
    end
end

function QOHamiltonian(H,paramrules,Ncutoff)
    oldops = collect(get_qsymbols(H))
    if check_hilberts(H)
        hilb = first(oldops).hilbert
    else
        error("Expression does not have homogenous Hilbert space")
    end

    #now we want to contruct a dictionary mapping the 'atomic' hilbert spaces of QC to the 'atomic' hilbert spaces of QOHamiltonian
    #first just get it working for only Fockstates
    subspaces = []
    for space in hilb.spaces
        if space isa FockSpace{Symbol}
            push!(subspaces,FockBasis(Ncutoff)) #TODO - how can we embed Ncutoff in the QC atomic hilbert space?
        end
    end

    #next, create a QO operator from each QC operator
    newops = []
    for oldop in oldops
        newsubopfuncs = [op_selector(oldop,idx) for idx in eachindex(subspaces)]
        newsubops = [f(space) for (f,space) in zip(newsubopfuncs,subspaces)]
        newop = tensor(newsubops...)
        push!(newops,newop)
    end

    #now we need to substitute these operators into the symbolic expression of the Hamiltonian
    oprules = Dict([old => new for (old,new) in zip(oldops,newops)])
    rules = merge(oprules,paramrules)
    new_H = substitute(H,rules)

    return new_H
end

struct SLH
    inputs #must have unique elements
    outputs #must have unique elements
    S #size nxn
    L #size n
    H #has operators which act on hilbert
    hilbert 
    operators # contains a list of all operators appearing above
end

SLH(inout, S, L, H, hilbert, operators) = SLH(inout, inout, S,L,H,hilbert,operators)

#Wrote this with chatgpt
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


#This function is for QuantumCumulants operators
function promote(operator,product_space)
    #we identify the hilbert space which our operator acts on, which should be a subspace of the product space
    
    if hasproperty(operator.hilbert,:spaces) #then it is a product space
        subspace = operator.hilbert.spaces[operator.aon]
    else #it is not a product space
        subspace = operator.hilbert
    end

    #this identifies the operator's hilbert space as a subspace of the product space
    subspaceindex = findfirst(isequal(subspace),product_space.spaces)

    if isnothing(subspaceindex)
        error("$operator does not act on a subspace of $product_space")
    end

    #these next two lines grab the necessary data to construct a new version of the operator on the product space
    #in the case of a Fock space this is just the name of the operator
    #for an NLevelSpace this is the operator name and the names of the two levels it transitions between
    middlefieldnames = fieldnames(typeof(operator))[2:end-2]
    middlefields = [getfield(operator,name) for name in middlefieldnames]

    #this calls the operator constructor with the old 'middle data' but on the larger hilbert space
    return typeof(operator).name.wrapper(product_space, middlefields...,subspaceindex)
end


#Here we write functions to implement the composition rules for SLH systems described in Combes (2017)
function concatenation(A,B)
    #To create a combined system, we 'stack' the inputs and outputs of A on top of those of B
    inputs = cat(A.inputs,B.inputs,dims = 1)
    outputs = cat(A.outputs,B.outputs,dims = 1)

    # since the matrices may be of any type, we cannot use the built-in function cat(A[1], B[1]; dims=(1,2)) since it tries to create a matrix of zeros from type Any 
    Z12 = zeros(size(A.S)[1],size(B.S)[1])
    Z21 = zeros(size(B.S)[1],size(A.S)[1])

    S1 = cat(A.S, Z12; dims=2)
    S2 = cat(Z21, B.S; dims=2)

    S = cat(S1,S2;dims=1)

    if !isnothing(A.hilbert) 
        if !isnothing(B.hilbert)
            hilb_product = tensor(A.hilbert,B.hilbert)

            new_Aops = [promote(op,hilb_product) for op in A.operators]
            new_Bops = [promote(op,hilb_product) for op in B.operators]

            Arule = Dict([old => new for (old,new) in zip(A.operators,new_Aops)])
            Adagrule = Dict([adjoint(old) => adjoint(new) for (old,new) in zip(A.operators,new_Aops)])
            Arules = merge(Arule,Adagrule)
            new_HA = substitute(A.H,Arules)

            Brule = Dict([old => new for (old,new) in zip(B.operators,new_Bops)])
            Bdagrule = Dict([adjoint(old) => adjoint(new) for (old,new) in zip(B.operators,new_Bops)])
            Brules = merge(Brule,Bdagrule)
            new_HB = substitute(B.H,Brules)

            ops = cat(new_Aops,new_Bops;dims=1)
            H = new_HB + new_HA

            new_LA = [substitute(collapse,Arules) for collapse in A.L]
            new_LB = [substitute(collapse,Brules) for collapse in B.L]

            L = cat(new_LA, new_LB, dims=1)
        else
            L = cat(A.L,B.L,dims = 1)
            H = A.H
            hilb_product = A.hilbert
            ops = A.operators
        end
    else
        if !isnothing(B.hilbert)
            L = cat(A.L,B.L,dims = 1)
            H = B.H
            hilb_product = B.hilbert
            ops = B.operators
        else 
            L = cat(A.L,B.L,dims = 1) 
            H = 0
            hilb_product = nothing
            ops = nothing
        end
    end

    return SLH(inputs, outputs,S,L,H,hilb_product,ops)
end

#x is the output port, y is the input port
function feedbackreduce(A,output, input)

    x = findfirst(isequal(output),A.outputs)
    y = findfirst(isequal(input), A.inputs)

    newoutputs = A.outputs[eachindex(A.outputs) .!= x]
    newinputs = A.inputs[eachindex(A.inputs) .!= y]

    Sxbarybar = A.S[1:end .!= x, 1:end .!= y]
    Sxbary = A.S[1:end .!= x, y]
    Sxybar = permutedims(A.S[x, 1:end .!= y])
    Sxy = A.S[x,y]

    Sy = A.S[:,y]

    Lxbar = A.L[1:end .!= x]
    Lx = A.L[x]

    S = Sxbarybar + Sxbary*(1-Sxy)^(-1)*Sxybar
    #Have to use fill here because the usual broadcast syntax on operator does not work
    L = Lxbar +  Sxbary .* fill(((1-Sxy)^(-1)*Lx),size(Sxbary))
    #println((dagger(Sy)*A.L)[1])

    term1 = adjoint(A.L)*Sy
    term2 = (1-Sxy)^(-1)
    term3 = Lx
    termA = adjoint(Lx)
    termB = ((1-adjoint(Sxy))^(-1))
    termC = (adjoint(Sy)*A.L)

    Hprime = 1/(2im)*(term1*term2*term3-termA*termB*termC)

    H = A.H + Hprime

    return SLH(newinputs,newoutputs,S,L,H,A.hilbert,A.operators)
end


end