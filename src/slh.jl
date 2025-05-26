"""
SLH(inputs, outputs, S, L, H)

An SLH triple describes an open quantum system. See Combes, arXiv.1611.00375
"""
struct SLH
    inputs #must have unique elements
    outputs #must have unique elements
    S #size nxn
    L #size n
    H #has operators which act on hilbert
end

SLH(inout, S, L, H) = SLH(inout, inout, S,L,H)

function hilbert(sys::SLH)
    ops = operators(sys) 
    if check_hilberts(sum(ops)) #Have to sum because check_hilberts() is expecting an expression
        return first(ops).hilbert
    else
        return error("non-unique hilbert space")
    end
end

function operators(sys)
    return get_qsymbols(sys.H)
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


"""
concatenation(A::SLH,B::SLH)

creates a composite system of A and B with no interconnections. Combes eq. 59
"""
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

    if !isnothing(hilbert(A)) 
        if !isnothing(hilbert(B))
            hilb_product = QuantumCumulants.tensor(hilbert(A),hilbert(B))

            new_Aops = [promote(op,hilb_product) for op in operators(A)]
            new_Bops = [promote(op,hilb_product) for op in operators(B)]

            Arule = Dict([old => new for (old,new) in zip(operators(A),new_Aops)])
            Adagrule = Dict([adjoint(old) => adjoint(new) for (old,new) in zip(operators(A),new_Aops)])
            Arules = merge(Arule,Adagrule)
            new_HA = substitute(A.H,Arules)

            Brule = Dict([old => new for (old,new) in zip(operators(B),new_Bops)])
            Bdagrule = Dict([adjoint(old) => adjoint(new) for (old,new) in zip(operators(B),new_Bops)])
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
            hilb_product = hilbert(A)
            ops = operators(A)
        end
    else
        if !isnothing(hilbert(B))
            L = cat(A.L,B.L,dims = 1)
            H = B.H
            hilb_product = hilbert(B)
            ops = operators(B)
        else 
            L = cat(A.L,B.L,dims = 1) 
            H = 0
            hilb_product = nothing
            ops = nothing
        end
    end

    return SLH(inputs, outputs,S,L,H)
end

"""
feedbackreduce(A::SLH,output,input)

Connects the output port to the input port, reducing the number of outputs and inputs by one each. Combes eq 61.
"""
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

    return SLH(newinputs,newoutputs,S,L,H)
end


