using QuantumCumulants
using SystemLevelHamiltonian
using ModelingToolkit

################################################
#system definition
h_arm = FockSpace(:fabry)
h_src = FockSpace(:fabry)
h_atom = NLevelSpace(:atom, (:g,:e))

hilb = h_arm⊗h_src⊗h_atom

@qnumbers σ::Transition(hilb,3) a::Destroy(hilb,1) b::Destroy(hilb,2)

ps = @cnumbers Δ g h Ω κa κb Γ
p0 = (4, 1, 0.1,0.01,0.1,0.3,0.3)


H = Δ*σ(:e,:e) + g*(b'*σ(:g,:e) + b*σ(:e,:g)) + Ω*(b'*a + b*a') + (a + a')*h

L = [sqrt(κa)*a, sqrt(κb)*b,Γ*σ(:g,:e)]

ops = [σ(:e,:e),a,b]
T_final = 100

ifo = SLH(:sys, [:in_a,:in_b,:in_atom],[:out_a,:out_b,:out_atom],[1],L,H)
###################################################################################


#Define the filter cavity 
hilb = FockSpace(:cavA)
a = Destroy(hilb,Symbol(:a,:_,:cavA))

@syms t::Real
@register_symbolic κ_A(t)
κ_A(t) = exp(-(t-10)^2/4)

Δ_A = rnumbers(:Δ_A)

filter_A = SLH(:cavA,
                [:Input],
                [:Output],
                [1],
                [κ_A(t)*a],
                Δ_A*adjoint(a)*a)
