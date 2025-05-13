using SystemLevelHamiltonian
using QuantumCumulants
using GLMakie
using ModelingToolkit, OrdinaryDiffEq 

#############################################################
#Define the filter cavity 
h_fabry = FockSpace(:fabry)
@qnumbers a::Destroy(h_fabry)

@rnumbers κ Δω

FPcavity = SLH(["LEFT","RIGHT"],
                [1 0;
                 0 1],
                [sqrt(κ)*a;
                 sqrt(κ)*a],
                Δω*adjoint(a)*a)

#############################################################
#define the coherent drive (laser)
h_laser = FockSpace(:laser)
@qnumbers b::Destroy(h_laser)
@syms t::Real
@register_symbolic α(t)
α(t) = exp(-(t-5)^2/4)

Laser = SLH(["BRIGHT"],
            [1],
            [α(t)],
            1e-50*b'*b)

#############################################################

#############################################################
h_collector = FockSpace(:collector)
@qnumbers c::Destroy(h_collector)
Ω = 1
@syms t::Real
@register_symbolic β(t)
β(t) = sin(Ω*t)

collector = SLH(["INPUT"],
            [1],
            [β(t)],
            1e-50*c'*c)

#connect the subsystems

system = concatenation(Laser,FPcavity)
system = concatenation(system, collector)

system = feedbackreduce(system,"BRIGHT","LEFT")
system = feedbackreduce(system,"RIGHT","INPUT")

allops = collect(operators(system))
A = allops[4] # <- bad
C = allops[3]

ops = [A,A',A'*A,C,C'*C]

eq_n = meanfield(ops,system.H,system.L; rates=ones(length(system.L)))
eqs = complete(eq_n)

#Numerical solution

#All the parameters we will pass to the ODE
ps = (κ, Δω, α)


@named sys = ODESystem(eqs)

paramwidget(sys,ops,30)