using QuantumToolbox
using SystemLevelHamiltonian


function final_state(A)
    #Example from QT docs

    N = 4 # Fock space truncated dimension

    ωa = 1
    ωc = 1 * ωa    # considering cavity and atom are in resonance
    σz = sigmaz() ⊗ qeye(N) # order of tensor product should be consistent throughout
    a  = qeye(2)  ⊗ destroy(N)  
    Ω  = 0.05
    σ  = sigmam() ⊗ qeye(N)

    Ha = ωa / 2 * σz
    Hc = ωc * a' * a # the symbol `'` after a `QuantumObject` act as adjoint
    Hint = Ω * (σ * a' + σ' * a)
    Hdrive = (a' + a)

    omega = 2

    coef(p, t) = A*sin(omega*t)

    H_t = QobjEvo(Hdrive, coef)

    Htot  = Ha + Hc + Hint +H_t
    e_ket = basis(2,0) 
    ψ0 = e_ket ⊗ fock(N, 0)

    tlist = 0:2.5:1000 # a list of time points of interest

    κ = 0.3
    L = sqrt(κ) * a

    sol_me  = mesolve(Htot,  ψ0, tlist,[L])
    return sol_me.states[end]
end

function derivative(A,dA)
    rho1 = hermitian_data(final_state(A))
    rho2 = hermitian_data(final_state(A+dA))
    
    return (rho1, rho2, (rho2 - rho1)/dA)
end

(rho1, rho2, rhodot) = derivative(3,1)

L = sld_operator(rho1,rhodot)

qfi = tr(L*L*rho1)