using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using GLMakie

##############################################################
### SYSTEM DEFINITION
# Define parameters


# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:s,:e))
h = hf ⊗ ha

# Define the fundamental operators
ℰ = Destroy(h,:a)
σ(i,j) = Transition(h, :σ, i, j)

# Hamiltonian
@cnumbers Δ Ω g
H = Δ*σ(:e,:e) - (Ω*σ(:e,:s)+Ω'*σ(:s,:e) +g*ℰ*σ(:e,:g)+ g*ℰ'*σ(:g,:e))

# Collapse operators
#J = [a,ge,ge']
#rates = [κ,γ,ν]