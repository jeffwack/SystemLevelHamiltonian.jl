using QuantumCumulants
using SystemLevelHamiltonian
using QuantumOptics
using ModelingToolkit, OrdinaryDiffEq
using GLMakie

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
###################################################################################

##################################################################################
#Simulate with QuantumCumulants

eq_n = meanfield(ops,H,L;rates=ones(size(L)),order=3)
eqs = complete(eq_n)

@named sys = ODESystem(eqs)

#NOTE - you must set the type of these numbers to be ComplexF64 to avoid error
u0 = zeros(ComplexF64,length(eqs))#TODO: Are some initial conditions unphysical? For example not all combos of expected values of spin operators should be allowed?

prob = ODEProblem(sys,u0,(0.0,T_final),ps.=>p0)
sol = solve(prob,RK4())
###################################################################################

####################################################################################
#Simulate with QuantumOptics
Ncutoff = 15

QOH = convert_to_QO([H],Dict(ps.=>p0),Ncutoff)[1] #Making a list with one element and extracting the one element is because convert_to_QO() is expecting a list
QOL = convert_to_QO(L,Dict(ps.=>p0),Ncutoff)
QOops = convert_to_QO(ops,Dict(ps.=>p0),Ncutoff)

Ψ₀ = standard_initial_state(basis(QOH))
ρ₀ = Ψ₀ ⊗ dagger(Ψ₀)

N_steps = 500
T = range(0,T_final,N_steps)
tout, ρt_master = timeevolution.master(T, ρ₀, QOH, QOL)

#Now calculate the expected values
traj = []
for op in QOops
    expectop(ρ) = expect(op,ρ)
    push!(traj,expectop.(ρt_master))
end
#############################################################################

#############################################################################
#Plotting
fig = Figure()
ax = Axis(fig[1,1])

curves = [Point2f.(sol.t,real.(sol[op])) for op in ops]
series!(ax,curves,labels = "QC".*string.(ops),color = :tab10)

curves = [Point2f.(tout,real.(traj[idx])) for idx in eachindex(traj)]
series!(ax,curves,labels = "QO".*string.(ops),color = :Dark2_8)

axislegend(ax)
fig