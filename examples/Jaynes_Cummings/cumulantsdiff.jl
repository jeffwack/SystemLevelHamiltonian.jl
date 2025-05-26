using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using ForwardDiff

##############################################################
### SYSTEM DEFINITION

# Define parameters and their numerical values
ps = @cnumbers Δ g κ h 
p0 = ComplexF64[0.1, 5, 0.3,0.1]

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
hilb = hf ⊗ ha

# Define the operators
a = Destroy(hilb,:a)
sm = Transition(hilb,:σ,:g,:e)
sp = sm'
sz = Transition(hilb,:σ,:e,:e)

# Hamiltonian
H = Δ*a'*a + g*(a'*sm + a*sp) + h*(a + a')

# Coupling operators
L = [κ*a]

small_ops = [a,a'*a,sz]

#####################################
#Simulate with QuantumCumulants

eq_n = meanfield(small_ops,H,L;rates=ones(size(L)),order=3)
eqs = complete(eq_n)

#operators we want to calculate expected values of
ops = eqs.operators

@named sys = ODESystem(eqs)

f_expr, other_expr = generate_function(sys, dvs = unknowns(sys), ps = parameters(sys),expression = Val{false})

f! = eval(f_expr)

#NOTE - you must set the type of these numbers to be ComplexF64 to avoid error
qcu0 = zeros(ComplexF64,length(eqs))#TODO: Are some initial conditions unphysical? For example not all combos of expected values of spin operators should be allowed?

function flatten_complex(x::Vector{ComplexF64})
    vcat(real.(x), imag.(x))
end

function fqc(params)
    _prob = ODEProblem(f!, qcu0, (0.0, 100.0), params)
    sol = solve(_prob, RK4(); saveat=[100.0])
    flatten_complex(Array(sol)[:, end])
end


dx = ForwardDiff.jacobian(fqc, collect(p0))