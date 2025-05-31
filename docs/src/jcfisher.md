# Calculating the QFI of a Jaynes-Cummings model
```@example jcfisher

using SystemLevelHamiltonian
using QuantumCumulants
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


qfi = compute_qfi(sys, Ncutoff,T, paramrules, h,AutoFiniteDiff())
```