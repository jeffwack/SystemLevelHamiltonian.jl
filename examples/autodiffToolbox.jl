#This script will try to implement the strategy of treating mesolve as a black box to differentiate through

using QuantumToolbox
using DifferentiationInterface
using SciMLSensitivity
using Enzyme

function final_state(p)
    f1(p, t) = p[1] * cos(p[2] * t)
    f2(p, t) = p[3] * sin(p[4] * t)
    γ(p, t)  = sqrt(p[5] * exp(-p[6] * t))

    H_t = sigmaz() + QobjEvo(sigmax(), f1) + QobjEvo(sigmay(), f2)

    c_ops = [
        QobjEvo(destroy(2), γ)
    ]

    ψ0 = basis(2, 0)
    tlist = 0:2.5:100

    sol = mesolve(H_t, ψ0, tlist, c_ops; params = p)
    return sol.states[end].data
end
Enzyme.API.strictAliasing!(false)
Enzyme.Compiler.VERBOSE_ERRORS[] = true
jac = DifferentiationInterface.jacobian(final_state,AutoEnzyme(),[0.1, 0.2, 0.3, 0.4, 0.5, 0.6])