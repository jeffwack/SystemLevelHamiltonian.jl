using QuantumCumulants
using SystemLevelHamiltonian
using QuantumOptics

##############################################################
### SYSTEM DEFINITION
# Define parameters
@cnumbers Δ g κ

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf ⊗ ha

# Define the fundamental operators
a = Destroy(h,:a)
s = Transition(h,:σ,:g,:e)

# Hamiltonian
H = Δ*a'*a + g*(a'*s + a*s')

# Collapse operators
J = [a]
rates = [κ]
