using QuantumCumulants
using SystemLevelHamiltonian
using QuantumOptics
using ModelingToolkit, OrdinaryDiffEq
using GLMakie

##############################################################
### SYSTEM DEFINITION
# Define parameters
@cnumbers Δ g κ h ω

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
hilb = hf ⊗ ha

# Define the fundamental operators
a = Destroy(hilb,:a)
sm = Transition(hilb,:σ,:g,:e)
sp = sm'
#sz = Transition(hilb,:σ,:e,:e)

# Hamiltonian
H = Δ*a'*a + g*(a'*sm + a*sp) + h*(a + a')

# Collapse operators
J = [a]
rates = [κ]
#####################################
#Simulate with QuantumCumulants
ops = [a,a'*a,sz]
eq_n = meanfield(ops,H,J;rates=rates,order=3)
eqs = complete(eq_n)

#All the parameters we will pass to the ODE
ps = (Δ, g, κ, h)

@named sys = ODESystem(eqs)

#NOTE - you must set the type of these numbers to be ComplexF64 to avoid error
u0 = zeros(ComplexF64,length(eqs))#TODO: how can we set physical initial conditions? First random stab led to casting complex to real error

p0 = (0.1, 5, 0.3,0.1)
prob = ODEProblem(sys,u0,(0.0,100),ps.=>p0)
sol = solve(prob,RK4())

plotops(sol,ops)