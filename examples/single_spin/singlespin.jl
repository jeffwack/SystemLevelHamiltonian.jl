using QuantumToolbox
using SystemLevelHamiltonian
using QuantumCumulants
using GLMakie



# Define parameter
@cnumbers h

# Hilbert space for qubit
hilb = NLevelSpace(:spin, (:e, :g))

# Pauli operator
σz = Transition(hilb,:σ,:e,:e)
σx = Transition(hilb,:σ,:g,:e)

# Hamiltonian: H = λ * σz
H = h * σx'

# No dissipation
L = [0.0001*σx]

sys = SLH(:spin,nothing,nothing,[1],L,H)

Ncutoff = 10
N_steps = 1000
T = range(0,100,N_steps)

paramrules = Dict([h=>2.0])

#numsys = convert_to_QT(sys,Ncutoff,paramrules)

qfiviaL = compute_qfi_fdm(sys, Ncutoff,T, paramrules, h)

t_end = range(3,1000,100)

qfi = []
for t in t_end
    T = range(0,t,N_steps)
    push!(qfi,compute_qfi_fdm(sys, Ncutoff,T, paramrules, h))
end

plot(t_end,real.(qfi))