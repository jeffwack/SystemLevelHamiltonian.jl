
using QuantumCumulants
using ModelingToolkit
using OrdinaryDiffEqLowOrderRK

using DifferentiationInterface
using SciMLSensitivity

using Enzyme
using FiniteDiff
using Zygote
using Mooncake
using ForwardDiff

ps = @cnumbers Δ g κ h

hilb = FockSpace(:cavity)⊗NLevelSpace(:spin,(:g, :e))

a = Destroy(hilb,:a)
sp = Transition(hilb, :σ, :e, :g)
sm = sp'

H = Δ*a'*a + g*(a'*sm + a*sp) + h*(a+a')

L = [κ*a]

small_ops = [a, a'*a, Transition(hilb, :σ, :e, :e)]

eq_n = meanfield(small_ops, H, L; rates = ones(size(L)), order=3)
eqs = complete(eq_n)

@named sys = ODESystem(eqs)

u0 = zeros(ComplexF64,length(eqs))

function finalstate(p0)    
    prob = ODEProblem(sys,u0,(0.0,100),ps .=> p0)
    sol = solve(prob,RK4())
    return sol.u[end]
end

jac = DifferentiationInterface.jacobian(finalstate,AutoForwardDiff(),[0.3,0.2,0.4,0.1])
#jac = DifferentiationInterface.jacobian(finalstate,AutoMooncake(;config = nothing),[0.3,0.2,0.4,0.1])
