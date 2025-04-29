using ModelingToolkit, StochasticDiffEq
using ModelingToolkit: t_nounits as t,  D_nounits as D
using GLMakie

@variables x(t) = 0.0
@brownian B

eqs = [D(x) ~ B]

@mtkbuild de = System(eqs, t)
prob = SDEProblem(de, [], (0.0, 1000.0), [])
sols = []
for ii in 1:100
    push!(sols,solve(prob, SRIW1()))
end

fig = Figure()
ax = Axis(fig[1,1])
for sol in sols
    series!(ax,[Point2f.(sol.t,[x[1] for x in sol.u])]) #note how annoying it was to use series.
end
fig
