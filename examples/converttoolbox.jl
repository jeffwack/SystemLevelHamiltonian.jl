using QuantumCumulants
using SystemLevelHamiltonian
using QuantumToolbox
using ModelingToolkit, OrdinaryDiffEq
using GLMakie

##############################################################
### SYSTEM DEFINITION

# Define parameters and their numerical values
ps = @cnumbers Δ g κ h 
p0 = (0.1, 5, 0.3,0.1)

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
hilb = QuantumCumulants.tensor(hf,ha)

# Define the operators
a = Destroy(hilb,:a)
sm = Transition(hilb,:σ,:g,:e)
sp = sm'
sz = Transition(hilb,:σ,:e,:e)

# Hamiltonian
H = Δ*a'*a + g*(a'*sm + a*sp) + h*(a + a')

# Coupling operators
L = [κ*a]

small_ops = [a,a'*a,sz]
###########################################################


#####################################
#Simulate with QuantumCumulants

eq_n = meanfield(small_ops,H,L;rates=ones(size(L)),order=3)
eqs = complete(eq_n)

#operators we want to calculate expected values of
ops = eqs.operators

@named sys = ODESystem(eqs)

#NOTE - you must set the type of these numbers to be ComplexF64 to avoid error
u0 = zeros(ComplexF64,length(eqs))#TODO: Are some initial conditions unphysical? For example not all combos of expected values of spin operators should be allowed?

prob = ODEProblem(sys,u0,(0.0,100),ps.=>p0)
sol = solve(prob,RK4())
#######################################################

#######################################################
#Simulate with QuantumToolbox #TODO: get rid of SystemLevelHamiltonian dep for conversion
Ncutoff = 10

QTH = convert_to_QT([H],Dict(ps.=>p0),Ncutoff)[1] #Making a list with one element and extracting the one element is because convert_to_QO() is expecting a list
QTL = convert_to_QT(L,Dict(ps.=>p0),Ncutoff)
QTops = convert_to_QT(small_ops,Dict(ps.=>p0),Ncutoff)

Ψ₀ = standard_initial_state(QTH)

N_steps = 1000
T = range(0,100,N_steps)
sol_me = mesolve(QTH,Ψ₀,T, QTL) #Note that you pass a pure initial state, not a density matrix

ρt_master = sol_me.states
#Now calculate the expected values
traj = []
for op in QTops
    expectop(ρ) = expect(op,ρ)
    push!(traj,expectop.(ρt_master))
end
##################################################################################

##################################################################################
#Plotting
fig = Figure()
ax = Axis(fig[1,1])

curves = [Point2f.(sol.t,real.(sol[op])) for op in small_ops]
series!(ax,curves,labels = "QC".*string.(ops),color = :batlowS)

tout = sol_me.times
curves = [Point2f.(tout,real.(traj[idx])) for idx in eachindex(traj)]
series!(ax,curves,labels = "QT".*string.(ops),color = :batlowW)

axislegend(ax)
fig

