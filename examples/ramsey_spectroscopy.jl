using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkit
using CairoMakie

@cnumbers Δ Ω Γ ν
@syms t::Real
@register_symbolic f(t)

# Hilbert space
h = NLevelSpace(:atom,2)

# operator
σ(i,j) = Transition(h, :σ, i, j, 2)

# Hamiltonian
H = -Δ*σ(2,2) + f(t)*Ω/2*(σ(1,2) + σ(2,1))

# Jump operators & rates
J = [σ(1,2), σ(2,2)]
rates = [Γ, ν]

eqs = meanfield([σ(2,2), σ(1,2)],H,J;rates=rates)

@named sys = ODESystem(eqs)

# Parameter
Γ_ = 1.0
Ω_ = 500Γ_
Δ_ = 0Γ_
ν_ = 0.2Γ_

tp = π/2Ω_ # π/2-pulse
tf = 1/20Γ_ # free evolution without drive

function f(t)
    if t<tp || (t>tp+tf && t<2tp+tf)
        return 1
    else
        0
    end
end

ps = [Γ; Ω; Δ; ν]
p0 = [Γ_; Ω_; Δ_; ν_]

# Initial state
u0 = zeros(ComplexF64, length(eqs))

prob = ODEProblem(sys,u0,(0.0, 2tp+tf), ps.=>p0)
sol = solve(prob,Tsit5(),maxiters=1e7)

# Plot time evolution
t = sol.t
s22 = real.(sol[σ(2,2)])
plot(t, s22, xlabel="tΓ", ylabel="⟨σ22⟩", legend=false, size=(600,300))