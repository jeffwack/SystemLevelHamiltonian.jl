using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq 

#define the Hilbert space and operators
hilb = FockSpace(:DARM) ⊗ FockSpace(:SRC) ⊗ NLevelSpace(:atom,(:g,:e))

a = Destroy(hilb,:a,1)
b = Destroy(hilb,:b,2)
σz = Transition(hilb,:s,:e,:e,3)
σplus = Transition(hilb,:s,:e,:g)
σminus = Transition(hilb,:s,:g,:e)

#some parameters
@cnumbers Δ
@cnumbers Ω 

@syms t::Real
@register_symbolic h(t)
@register_symbolic g(t)

#h(t) = exp(-(t-10)^2/4)
#g(t) = h(t)

#now for the Hamiltonian. This is in the rotating frame at the frequency of the two cavities, which are at the same frequency
H = Δ*σz + g(t)*(b'*σminus + b*σplus) + Ω*(b'*a + b*a') + (a + a')*h(t)
#=
J = [a,b]
@cnumbers κa κb
rates = [κa, κb]

eq_n = meanfield(σz,H,J;rates=rates,order=3)

eqs = complete(eq_n)

# Numerical solution

#All the parameters we will pass to the ODE
ps = (Δ, Ω, κa, κb)

@named sys = ODESystem(eqs)

#Initial state
u0 = zeros(ComplexF64, length(eqs))

p0 = (0.3, 0.2, 0.7,0.3)
prob = ODEProblem(sys,u0,(0.0,100.0),ps.=>p0)
sol = solve(prob,RK4())
=#