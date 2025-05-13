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
@cnumbers α

Laser = SLH(["BRIGHT"],
            [1],
            [α],
            1e-50*b'*b)

#############################################################
#connect the subsystems

sys = concatenation(Laser,FPcavity)

sys = feedbackreduce(sys,"BRIGHT","LEFT")

A = collect(operators(sys))[2] # <- bad

X1 = 1/sqrt(2)*(A'+A)
X2 = 1.0im/sqrt(2)*(A'-A)


ops = [A,A',A'*A]

eq_n = meanfield(ops,sys.H,sys.L; rates=ones(length(sys.L)))
eqs = complete(eq_n)

#Numerical solution

#All the parameters we will pass to the ODE
ps = (κ, Δω, α)


@named sys = ODESystem(eqs)

paramwidget(sys,ops,10)