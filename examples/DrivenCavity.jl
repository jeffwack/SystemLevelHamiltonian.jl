include("SLH.jl")

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
                Δω*adjoint(a)*a,
                h_fabry,
                [a])

#############################################################
#define the coherent drive (laser)

@cnumbers α

Laser = SLH(["BRIGHT"],
            [1],
            [α],
            0,
            nothing,
            nothing)

#############################################################
#connect the subsystems

sys = concatenation(Laser,FPcavity)

sys = feedbackreduce(sys,"BRIGHT","LEFT")

A = sys.operators[1]

X1 = 1/sqrt(2)*(A'+A)
X2 = 1.0im/sqrt(2)*(A'-A)


ops = [X1,X2,a'*a]

eq_n = meanfield(ops,sys.H,sys.L; rates=ones(length(sys.L)))
eqs = complete(eq_n)

#Numerical solution

#All the parameters we will pass to the ODE
ps = (κ, Δω, α)


using ModelingToolkit, OrdinaryDiffEq 
@named sys = ODESystem(eqs)

#Initial state
u0 = zeros(ComplexF64, length(eqs))

p0 = (0.2, 0.1, 0.7+0.1im)
prob = ODEProblem(sys,u0,(0.0,100.0),ps.=>p0)
sol = solve(prob,RK4())

include("../QuantumMakie.jl")

plotall(sol)
