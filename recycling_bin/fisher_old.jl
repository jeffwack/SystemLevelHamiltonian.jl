"""
This script uses the finite difference method to calculate a quantum fisher information using two neighboring master equation problems
"""
using QuantumCumulants
using SystemLevelHamiltonian
using QuantumToolbox
using ModelingToolkit, OrdinaryDiffEq

##############################################################
### SYSTEM DEFINITION

# Define parameters and their numerical values
ps = @cnumbers Δ g κ h 

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
hilb = QuantumCumulants.tensor(hf,ha)

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
###########################################################
Ncutoff = 10
N_steps = 1000
T = range(0,100,N_steps)

p = (0.1, 5, 0.3,0.1)
dh = 1e-6
#######################################################
#Simulate with QuantumToolbox #TODO: get rid of SystemLevelHamiltonian dep for conversion

p1 = (0.1, 5, 0.3, 0.1+dh/2)
p2 = (0.1, 5, 0.3, 0.1-dh/2)

QTH1 = convert_to_QT([H],Dict(ps.=>p1),Ncutoff)[1] #Making a list with one element and extracting the one element is because convert_to_QO() is expecting a list
QTL1 = convert_to_QT(L,Dict(ps.=>p1),Ncutoff)

QTH2 = convert_to_QT([H],Dict(ps.=>p2),Ncutoff)[1] #Making a list with one element and extracting the one element is because convert_to_QO() is expecting a list
QTL2 = convert_to_QT(L,Dict(ps.=>p2),Ncutoff)

Ψ₀ = standard_initial_state(QTH1)

sol_me1 = mesolve(QTH1,Ψ₀,T, QTL1) #Note that you pass a pure initial state, not a density matrix
sol_me2 = mesolve(QTH2,Ψ₀,T, QTL2) #Note that you pass a pure initial state, not a density matrix

ρt_master1 = sol_me1.states
ρt_master2 = sol_me2.states
##################################################################################

rho1 = ρt_master1[end]
rho2 = ρt_master2[end]

rho_dot = (rho1 - rho2)/dh

rho = (rho1+rho2)/2

ρ = hermitian_data(rho)
ρ_dot = hermitian_data(rho_dot)

tol = 1e-4

L = sld_operator(ρ,ρ_dot,eps = tol)

qfiviaL = tr(ρ*L*L)

qfi2 = compute_qfi(ρ,ρ_dot,eps = tol)

qfi3 = compute_qfi_alt(ρ,ρ_dot,eps = tol)

println(qfiviaL)
println(qfi2)
println(qfi3)


