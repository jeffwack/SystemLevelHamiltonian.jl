using QuantumCumulants
using SystemLevelHamiltonian
using QuantumToolbox
using ModelingToolkit, OrdinaryDiffEq
using GLMakie

##############################################################
### SYSTEM DEFINITION

# Define parameters and their numerical values
ps = @cnumbers Δ g κ h 
dh = 0.001
p1 = (0.1, 5, 0.3,0.1)
p2 = (0.1, 5, 0.3,0.1+dh)

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
#######################################################
#Simulate with QuantumToolbox #TODO: get rid of SystemLevelHamiltonian dep for conversion
Ncutoff = 10

QTH1 = convert_to_QT([H],Dict(ps.=>p1),Ncutoff)[1] #Making a list with one element and extracting the one element is because convert_to_QO() is expecting a list
QTL1 = convert_to_QT(L,Dict(ps.=>p1),Ncutoff)

QTH2 = convert_to_QT([H],Dict(ps.=>p2),Ncutoff)[1] #Making a list with one element and extracting the one element is because convert_to_QO() is expecting a list
QTL2 = convert_to_QT(L,Dict(ps.=>p2),Ncutoff)

Ψ₀ = standard_initial_state(QTH)

N_steps = 1000
T = range(0,100,N_steps)
sol_me1 = mesolve(QTH1,Ψ₀,T, QTL1) #Note that you pass a pure initial state, not a density matrix
sol_me2 = mesolve(QTH2,Ψ₀,T, QTL2) #Note that you pass a pure initial state, not a density matrix

ρt_master1 = sol_me1.states
ρt_master2 = sol_me2.states
##################################################################################
T_f = 600

rho1 = ρt_master1[T_f]
rho2 = ρt_master2[T_f]

rho_dot = (rho2 - rho1)/dh

rho = (rho1+rho2)/2


using LinearAlgebra

#Function from James, translated to Julia by chatGPT
function compute_qfi(qρ, qρ_dot)
    ρ = qρ.data
    ρ_dot = qρ_dot.data
    evals, evecs = eigen(Hermitian(ρ))  # eigendecomposition

    FQ_element = zeros(eltype(ρ), size(ρ))

    for j in eachindex(evals)
        λ_j = evals[j]
        ψ_j = evecs[:, j]
        for k in eachindex(evals)
            λ_k = evals[k]
            ψ_k = evecs[:, k]

            numerator = 4 * λ_j * abs(dot(ψ_k', ρ_dot * ψ_j))^2
            denominator = (λ_k + λ_j)^2
            FQ_element[j, k] = numerator / denominator
        end
    end

    FQ = sum(skipmissing(FQ_element))  # equivalent to np.nansum
    return real(FQ)
end

qfi = compute_qfi(rho,rho_dot)

#TODO: plot qfi as a function of time. What is the scaling with T?
#TODO: isolate qfi of subsystems with reduced density matrix. Verify 'reverse triangle inequality' of subsystem QFI
#TODO: using the SLD, perform the optimal measurement and calculate the classical Figure
#TODO: compare 1) 'direct' QFI 2) QFI from SLD 3) CFI of optimal measurment. If all 3 agree we feel good.
#TODO: if the states are pure and evolution is unitary, QFI is the variance of the generator of the unitary.
#TODO: Calculate time series for linear systems. Here the 'local time' measurements should not lose on any info. 
#Compare to the analytic Laplace transform. 