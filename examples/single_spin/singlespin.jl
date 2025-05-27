using QuantumToolbox
using SystemLevelHamiltonian
using QuantumCumulants
using OrdinaryDiffEq
using LinearAlgebra
using ModelingToolkit

println("Quantum Fisher Information for spin-½ in constant magnetic field H = λ σz")
println("---------------------------------------------------------------")

# Define parameter
@cnumbers λ

# Hilbert space for qubit
hs = NLevelSpace(:spin, (:up, :down))

# Pauli operators
σz = Transition(hilb,:σ,:e,:e)

# Hamiltonian: H = λ * σz
H = λ * σz

# No dissipation
L = []

# Simulation settings
λ₀ = 1.0
dh = 1e-6
tol = 1e-8
Ncutoff = 2  # Not used but required by the interface

tvals = [0.5, 1.0, 2.0]


for t_final in tvals
    println("\n--- Time t = $t_final ---")

    tspan = (0.0, t_final)
    times = range(tspan[1], tspan[2], length=100)

    # QFI via finite-difference method
    qfi_fd = compute_qfi_fdm(H, L, [λ], Ncutoff, ψ₀, times, t_final, [λ₀], 1, dh, eps=tol)

    # Manually evolve ψ₀ at λ ± dh
    function evolve(λval)
        Hλ = substitute(H, Dict(λ => λval))
        Hexpr = flatten(Hλ)
        model = quantum_schrodinger_equation(Hexpr, ψ₀)
        prob = ODEProblem(model, ψ₀, tspan)
        sol = solve(prob, Tsit5(), saveat=t_final)
        return sol[end]
    end

    ψ₊ = evolve(λ₀ + dh)
    ψ₋ = evolve(λ₀ - dh)

    ρ₊ = ψ₊ * dagger(ψ₊)
    ρ₋ = ψ₋ * dagger(ψ₋)

    ρ̄ = (ρ₊ + ρ₋)/2
    ρ̇ = (ρ₊ - ρ₋)/(2dh)

    # Compute QFI using different methods
    L_op = sld_operator(ρ̄, ρ̇, eps=tol)
    qfi_sld = real(tr(ρ̄ * L_op * L_op))
    qfi_alt = compute_qfi(ρ̄, ρ̇, eps=tol)

    # Analytical QFI
    qfi_exact = 4 * t_final^2

    # Print results
    println("Exact QFI     = $(round(qfi_exact, digits=6))")
    println("QFI (fdm)     = $(round(qfi_fd, digits=6))")
    println("QFI (SLD)     = $(round(qfi_sld, digits=6))")
    println("QFI (alt)     = $(round(qfi_alt, digits=6))")

    # Optional: print purity of ρ̄ to verify unitary evolution
    purity = real(tr(ρ̄ * ρ̄))
    println("Purity        = $(round(purity, digits=6))")
end