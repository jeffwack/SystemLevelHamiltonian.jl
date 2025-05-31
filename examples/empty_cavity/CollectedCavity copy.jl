#In this script I drive an empty cavity with a Gaussian wavelet pulse.

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
h_laser = FockSpace(:laser)
@qnumbers c::Destroy(h_laser)
@syms t::Real
@register_symbolic β(t)
Ω = 3
σ = 2
t_0 = 5
β(t) = exp(-(t-t_0)^2/σ^2)*cos(Ω*(t-t_0))

laser = SLH(["BRIGHT"],
            [1],
            [β(t)],
            1e-50*c'*c)

#connect the subsystems

system = concatenate(:sys,[FPcavity, laser])

system = feedbackreduce(system,"BRIGHT","LEFT")

allops = collect(operators(system))
A = allops[1] # <- bad
C = allops[4]

ops = [A,A',A'*A,C,C'*C]

eq_n = meanfield(ops,system.H,system.L; rates=ones(length(system.L)))
eqs = complete(eq_n)

#Numerical solution

#All the parameters we will pass to the ODE
ps = (κ, Δω)


@named sys = ODESystem(eqs)

paramwidget(sys,ops,30)