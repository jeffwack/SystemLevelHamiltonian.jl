using SciMLSensitivity
using OrdinaryDiffEq

function f(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + u[1]*u[2]
end

p = [1.5,1.0,3.0]
prob = ODEForwardSensitivityProblem(f,[1.0;1.0],(0.0,10.0),p)

sol = solve(prob,DP8())

x,dp = extract_local_sensitivities(sol)

#How do we put 