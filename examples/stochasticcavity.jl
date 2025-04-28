using ModelingToolkit, StochasticDiffEq
using ModelingToolkit: t_nounits as t,  D_nounits as D

@parameters Δ=3.0 γ=0.5
@variables rea(t) = 50 ima(t) = 0.0
@brownian B_inre B_inim

eqs = [D(rea) ~ -Δ*ima - γ/2*rea - sqrt(γ)*B_inre;
       D(ima) ~ Δ*rea - γ/2*ima - sqrt(γ)*B_inim]

@mtkbuild de = System(eqs, t)
prob = SDEProblem(de, [], (0.0, 100.0), [])
sol = solve(prob, SRIW1())

plot(sol)