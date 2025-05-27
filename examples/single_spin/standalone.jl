using QuantumToolbox
using LinearAlgebra
using FiniteDiff
using SystemLevelHamiltonian
using GLMakie

# --- Define constants ---
σz = sigmaz()
σx = sigmax()
I2 = Matrix{ComplexF64}(I, 2, 2)

# --- Initial state: |+⟩ along x ---
Ψ0 = 1/sqrt(2)*(basis(2,0) + basis(2,1))


t_list = range(0.1,2,100)
qfi = []
for t in t_list
    # --- Time at which to evaluate ---
    T = range(0,t,1000)

    # --- Time evolution of density matrix using sesolve ---
    function evolved_state(B::Float64)
        H = -0.5 * B * σz
        sol = sesolve(H, Ψ0, T)
        return ket2dm(sol.states[end])
    end

    # --- Compute derivative matrix with FiniteDiff ---
    B0 = 0.4
    rho_dot = FiniteDiff.finite_difference_derivative(evolved_state, B0)

    # --- Compute SLD and Fisher information ---
    rho = evolved_state(B0)
    L = sld_operator(hermitian_data(rho), hermitian_data(rho_dot))
    fisher_info = real(tr(hermitian_data(rho) * L * L))
    push!(qfi,fisher_info)
end

plot(t_list,real.(qfi))