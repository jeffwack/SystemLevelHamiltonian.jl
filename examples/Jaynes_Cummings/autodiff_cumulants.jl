"""
This script does not work! It is an attempt to perform autodifferentiation on a system of equations generated
by QuantumCumulants. The core issue is that the package I tried to use, ForwardDiff, does not suppost systems
of complex-valued equations. It seems like one should be able to 'simply' rephrase it as a system on 2n real
variables...

ERROR: ArgumentError: Cannot create a dual over scalar type ComplexF64. If the type behaves as a scalar, define ForwardDiff.can_dual(::Type{ComplexF64}) = true.
Stacktrace:
  [1] throw_cannot_dual(V::Type)
    @ ForwardDiff ~/.julia/packages/ForwardDiff/UBbGT/src/dual.jl:41
  [2] ForwardDiff.Dual{ForwardDiff.Tag{…}, ComplexF64, 4}(value::ComplexF64, partials::ForwardDiff.Partials{4, ComplexF64})
    @ ForwardDiff ~/.julia/packages/ForwardDiff/UBbGT/src/dual.jl:18
  [3] _broadcast_getindex_evalf
    @ ./broadcast.jl:678 [inlined]
  [4] _broadcast_getindex
    @ ./broadcast.jl:651 [inlined]
  [5] getindex
    @ ./broadcast.jl:610 [inlined]
  [6] macro expansion
    @ ./broadcast.jl:973 [inlined]
  [7] macro expansion
    @ ./simdloop.jl:77 [inlined]
  [8] copyto!
    @ ./broadcast.jl:972 [inlined]
  [9] copyto!
    @ ./broadcast.jl:925 [inlined]
 [10] materialize!
    @ ./broadcast.jl:883 [inlined]
 [11] materialize!
    @ ./broadcast.jl:880 [inlined]
 [12] seed!(duals::Vector{ForwardDiff.Dual{…}}, x::Vector{ComplexF64}, seeds::NTuple{4, ForwardDiff.Partials{…}})
    @ ForwardDiff ~/.julia/packages/ForwardDiff/UBbGT/src/apiutils.jl:52
 [13] vector_mode_dual_eval!(f::typeof(fqc), cfg::ForwardDiff.JacobianConfig{…}, x::Vector{…})
    @ ForwardDiff ~/.julia/packages/ForwardDiff/UBbGT/src/apiutils.jl:23
 [14] vector_mode_jacobian(f::typeof(fqc), x::Vector{…}, cfg::ForwardDiff.JacobianConfig{…})
    @ ForwardDiff ~/.julia/packages/ForwardDiff/UBbGT/src/jacobian.jl:129
 [15] jacobian
    @ ~/.julia/packages/ForwardDiff/UBbGT/src/jacobian.jl:22 [inlined]
 [16] jacobian(f::typeof(fqc), x::Vector{…}, cfg::ForwardDiff.JacobianConfig{…})
    @ ForwardDiff ~/.julia/packages/ForwardDiff/UBbGT/src/jacobian.jl:19
 [17] jacobian(f::typeof(fqc), x::Vector{ComplexF64})
    @ ForwardDiff ~/.julia/packages/ForwardDiff/UBbGT/src/jacobian.jl:19
 [18] top-level scope
    @ ~/.julia/dev/SystemLevelHamiltonian/examples/Jaynes_Cummings/cumulantsdiff.jl:60

"""

using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using SciMLSensitivity
using Zygote

##############################################################
### SYSTEM DEFINITION

# Define parameters and their numerical values
ps = @cnumbers Δ g κ h 
p0 = ComplexF64[0.1, 5, 0.3,0.1]

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
hilb = QuantumCumulants.tensor(hf, ha)

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

function fqc(params)
    _prob = ODEProblem(f!, qcu0, (0.0, 10.0), params)
    sol = solve(_prob, RK4())
    return Array(sol)[:, end]
end

function jacobi(f, x)
  y, back = Zygote.pullback(f, x)
  back(1)[1], back(im)[1]
end

function wirtinger(f, x)
  du, dv = jacobi(f, x)
  (du' + im*dv')/2, (du + im*dv)/2
end

dx = wirtinger(fqc, collect(p0))

"""
What is this output?? It has the same dimensionality as the parameters...
"""