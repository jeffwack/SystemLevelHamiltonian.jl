using QuantumCumulants
using SystemLevelHamiltonian
using QuantumToolbox
using ModelingToolkit, OrdinaryDiffEq
using DifferentiationInterface

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

sys = SLH(:jc,[:In],[:Out],[1],L,H)
###########################################################
Ncutoff = 10
N_steps = 1000
T = range(0,100,N_steps)

paramrules = Dict([Δ=>0.1,
                    g=>5.0,
                    h=>0.1,
                    κ=>0.3])


qfi_alt = compute_qfi(sys, Ncutoff,T, paramrules, h,AutoFiniteDiff())



#TODO: plot qfi as a function of time. What is the scaling with T?
#TODO: isolate qfi of subsystems with reduced density matrix. Verify 'reverse triangle inequality' of subsystem QFI
#TODO: using the SLD, perform the optimal measurement and calculate the classical Figure

#TODO: if the states are pure and evolution is unitary, QFI is the variance of the generator of the unitary.
#TODO: Calculate time series for linear systems. Here the 'local time' measurements should not lose on any info. 
#Compare to the analytic Laplace transform. 

