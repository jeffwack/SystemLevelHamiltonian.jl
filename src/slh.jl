"""
SLH(name, inputs, outputs, S, L, H)

An SLH triple describes an open quantum system. See Combes, arXiv.1611.00375

The name of the system should be unique. When multiple systems are combined, the names of their inputs and outputs will 
have the system name appended to them. The inputs and outputs describe 'ports' where signals leave and enter the system.
Quantum systems must have the same number of inputs and outputs, which we denote by n.

size(S) = (n, n) <- S is an nxn matrix
size(L) = (n,)
size(H) = ()

The two ways of combining SLH systems are concatenate() and feedbackreduce()
"""
struct SLH
    name #should be a symbol
    inputs #must have unique elements
    outputs #must have unique elements
    S #size nxn
    L #size n
    H #has operators which act on hilbert
end
#TODO: you should not be able to construct this if dimensions are wrong.

"""
hilbert(sys:SLH)

returns the hilbert space of the Hamiltonian of the given system.
"""
function hilbert(sys::SLH)
    ops = operators(sys) 
    if check_hilberts(sum(ops)) #Have to sum because check_hilberts() is expecting an expression
        return first(ops).hilbert
    else
        return error("non-unique hilbert space")
    end
end

"""
operators(sys)

returns all the quantum operators contained in the system's Hamiltonian.
"""
function operators(sys)
    return get_qsymbols(sys.H)
end

"""
parameters(sys)

returns all the symbolic numbers contained in the system's Hamiltonian and coupling vector L.
"""
function parameters(sys::SLH)
    return union(get_numsymbols(sys.H),get_numsymbols(sys.L...))
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
concatenate(name, syslist::Vector{SLH})

creates a composite system with no interconnections. Combes eq. 59
"""
function concatenate(name,syslist)
    #To create a combined system, we 'stack' the inputs and outputs of A on top of those of B

    #first, we promote the names of inputs and outputs, to prevent naming collisions
    newinputs = [[Symbol(input,:_,sys.name) for input in sys.inputs] for sys in syslist]
    inputs = cat(newinputs...,dims = 1)

    newoutputs = [[Symbol(output,:_,sys.name) for output in sys.outputs] for sys in syslist]
    outputs = cat(newoutputs...,dims = 1)

    #next, we concate all the scattering matrices block diagonally
    Slist = [sys.S for sys in syslist]
    S = cat(Slist...;dims=(1,2))

    #= #this used to be here, and I have a feeling will be useful to look at when that call to cat throws an error eventually
    #since the matrices may be of any type, we cannot use the built-in function cat(A[1], B[1]; dims=(1,2)) since it tries to create a matrix of zeros from type Any 
    Z12 = zeros(size(A.S)[1],size(B.S)[1])
    Z21 = zeros(size(B.S)[1],size(A.S)[1])

    S1 = cat(A.S, Z12; dims=2)
    S2 = cat(Z21, B.S; dims=2)

    S = cat(S1,S2;dims=1)
    =#
    hilb_product = QuantumCumulants.tensor([hilbert(sys) for sys in syslist]...)

    oldops = [collect(operators(sys)) for sys in syslist]
    newops = [[promote(op,hilb_product) for op in oplist] for oplist in oldops]
    

    rules = merge([Dict([old => new for (old,new) in zip(oldoplist,newoplist)]) for (oldoplist,newoplist) in zip(oldops,newops)]...)
    newHs = [substitute(sys.H,rules) for sys in syslist]

    H = sum(newHs)

    newLs = [[substitute(collapse,rules) for collapse in sys.L] for sys in syslist]

    L = cat(newLs, dims=1)

    return SLH(name,inputs, outputs,S,L,H)
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

    return SLH(A.name,newinputs,newoutputs,S,L,H)
end
"""
convert_to_QT(sys::SLH,paramrules,Ncuttof)

return (standard_initial_state(QTH), SLH(sys.inputs, sys.outputs, S, QTL,QTH)) 
"""
function convert_to_QT(sys::SLH,Ncutoff,paramrules)
    QTL = convert_to_QT([sys.L...,sys.H], Ncutoff,paramrules)
    #QTS = TODO
    QTH = pop!(QTL)
    name = Symbol(sys.name,:_,:numeric)

    return (SLH(name,sys.inputs, sys.outputs, sys.S, QTL,QTH),standard_initial_state(QTH))
end

"""
extract_creation_annihilation_operators(sys::SLH)

Extract all creation and annihilation operators from the system's Hamiltonian and coupling operators.
Returns a vector of operators ordered as [a₁, a₁†, a₂, a₂†, ...] for each bosonic mode.
"""
function extract_creation_annihilation_operators(sys::SLH)
    # Get all operators from H and L
    H_ops = get_qsymbols(sys.H)
    L_ops = union([get_qsymbols(L_term) for L_term in sys.L]...)
    all_ops = collect(union(H_ops, L_ops))
    
    # Group operators by their base mode (assuming they have a common structure)
    creation_ops = filter(op -> op isa Create, all_ops)
    annihilation_ops = filter(op -> op isa Destroy, all_ops)
    
    # Sort by mode index or name to ensure consistent ordering
    sort!(creation_ops, by = op -> op.name)
    sort!(annihilation_ops, by = op -> op.name)
    
    # Build the operator vector [a₁, a₁†, a₂, a₂†, ...]
    operators = []
    for (create_op, destroy_op) in zip(creation_ops, annihilation_ops)
        push!(operators, destroy_op)  # annihilation operator
        push!(operators, create_op)   # creation operator
    end
    
    return operators
end

"""
check_linearity(expr)

Check if an expression contains only linear terms (up to second order in bosonic operators).
Returns true if linear, false if nonlinear (third+ order products or fermionic operators).
"""
function check_linearity(expr)
    # TODO: Implement check for third+ order products of bosonic operators
    # TODO: Implement check for fermionic operators
    # For now, assume all systems are linear
    return true
end

"""
build_drift_matrix(H, operators)

Build the drift matrix A from the Hamiltonian H for the given operator ordering.
The A matrix describes how the operators evolve: d/dt [operators] = A * [operators] + ...
Throws an error if the system is nonlinear.
"""
function build_drift_matrix(H, operators)
    # Check if the system is linear
    if !check_linearity(H)
        throw(ArgumentError("This system is not linear and thus cannot be represented in state space form."))
    end
    
    n = length(operators)
    A = zeros(Symbolics.Num, n, n)
    
    # TODO: For each operator, compute its time derivative from the Hamiltonian
    # TODO: Use commutation relation: d/dt a = (i/ℏ)[H, a] with ℏ = 1
    # TODO: Express result in terms of basis operators to fill A matrix
    
    return A
end

"""
build_coupling_matrices(L, operators)

Build the coupling matrices B and C from the dissipation operators L.
B describes input-to-state coupling, C describes state-to-output coupling.
"""
function build_coupling_matrices(L, operators)
    n = length(operators)
    m = length(L)  # number of input/output channels
    
    B = zeros(Symbolics.Num, n, m)
    C = zeros(Symbolics.Num, m, n)
    
    # TODO: For each L operator, determine its coupling to the state vector
    # TODO: B matrix: how inputs couple into the system
    # TODO: C matrix: how system operators couple to outputs
    
    return B, C
end

"""
SLH2ABCD(sys::SLH)

Convert an SLH system to ABCD state-space representation.
Returns (A, B, C, D) matrices where:
- A: drift matrix (state evolution)
- B: input coupling matrix
- C: output coupling matrix  
- D: direct feedthrough matrix (from S matrix)
"""
function SLH2ABCD(sys::SLH)
    # Extract creation/annihilation operators
    operators = extract_creation_annihilation_operators(sys)
    
    # Build matrices
    A = build_drift_matrix(sys.H, operators)
    B, C = build_coupling_matrices(sys.L, operators)
    D = sys.S  # Direct feedthrough from scattering matrix
    
    return A, B, C, D
end

