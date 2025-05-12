using QuantumCumulants
using SystemLevelHamiltonian
using QuantumOptics
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
hilb = hf ⊗ ha

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

#########################################################
#Simulate with QuantumOptics #TODO: get rid of SystemLevelHamiltonian dep for conversion
Ncutoff = 10

QOH = convert_to_QO([H],Dict(ps.=>p0),Ncutoff)[1] #Making a list with one element and extracting the one element is because convert_to_QO() is expecting a list
QOL = convert_to_QO(L,Dict(ps.=>p0),Ncutoff)
QOops = convert_to_QO(ops,Dict(ps.=>p0),Ncutoff)

Ψ₀ = standard_initial_state(basis(QOH))
ρ₀ = Ψ₀ ⊗ dagger(Ψ₀)

N_steps = 1000
T = range(0,100,N_steps)
tout, ρt_master = timeevolution.master(T, ρ₀, QOH, QOL)

#Now calculate the expected values
traj = []
for op in QOops
    expectop(ρ) = expect(op,ρ)
    push!(traj,expectop.(ρt_master))
end
##################################################################################


##################################################################################
#=
#do Gram-Schmidt
using LinearAlgebra

opdata = [op.data for op in QOops]

#reshape 
opsmat = hcat(vec.(opdata)...)

#Orthonormalize 
Q = Matrix(qr(opsmat).Q)

#reshape vectors into operators
newopsdata = [reshape(v,22,22) for v in eachcol(Q)]

newops = [Operator(basis(QOH),data) for data in newopsdata]
#congrats, but you need to do the Gram-Schmidt in QC opjects to be very useful.
#We could use the R matrix to tell us how to combine expected values!
=#
#################################################################################
#Plotting
fig = Figure()
ax = Axis(fig[1,1])

curves = [Point2f.(sol.t,real.(sol[op])) for op in ops]
series!(ax,curves,labels = "QC".*string.(ops),color = :batlow)

curves = [Point2f.(tout,real.(traj[idx])) for idx in eachindex(traj)]
series!(ax,curves,labels = "QO".*string.(ops),color = :batlow)

axislegend(ax)
fig

