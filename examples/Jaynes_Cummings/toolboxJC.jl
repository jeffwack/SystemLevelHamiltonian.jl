using SystemLevelHamiltonian

using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq

import QuantumToolbox: sigmaz, sigmam, destroy, qeye, basis, fock, 
    ⊗, n_thermal, sesolve, mesolve

using GLMakie

#Example from QT docs

N = 2 # Fock space truncated dimension

ωa = 1
ωc = 1 * ωa    # considering cavity and atom are in resonance
σz = sigmaz() ⊗ qeye(N) # order of tensor product should be consistent throughout
a  = qeye(2)  ⊗ destroy(N)  
Ω  = 0.05
σ  = sigmam() ⊗ qeye(N)

Ha = ωa / 2 * σz
Hc = ωc * a' * a # the symbol `'` after a `QuantumObject` act as adjoint
Hint = Ω * (σ * a' + σ' * a)

Htot  = Ha + Hc + Hint

e_ket = basis(2,0) 
ψ0 = e_ket ⊗ fock(N, 0)

tlist = 0:2.5:1000 # a list of time points of interest

# define a list of operators whose expectation value dynamics exhibit Rabi oscillation
eop_ls = [
    a' * a,                      # number operator of cavity
    (e_ket * e_ket') ⊗ qeye(N), # excited state population in atom
]


γ = 4e-3
κ = 7e-3

L = [sqrt(κ)*a,sqrt(γ)* σ]

# `mesolve()` only has one additional keyword argument `c_ops` from `sesolve()`
sol_me  = mesolve(Htot,  ψ0, tlist,L , e_ops = eop_ls)

n_me = real.(sol_me.expect[1, :])
e_me = real.(sol_me.expect[2, :])

fig_me = Figure(size = (600, 350))
ax_me = Axis(
    fig_me[1, 1],
    xlabel = L"time $[1/\omega_a]$", 
    ylabel = "expectation value", 
    xlabelsize = 15, 
    ylabelsize = 15,
    width = 400,
    height = 220
)
lines!(ax_me, tlist, n_me, label = L"\langle a^\dagger a \rangle")
lines!(ax_me, tlist, e_me, label = L"$P_e$")
axislegend(ax_me; position = :rt, labelsize = 15)
display(fig_me);