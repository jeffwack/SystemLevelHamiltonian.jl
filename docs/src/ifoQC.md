# Linear Interferometer Readout - QuantumCumulants

```@example ifo
using QuantumCumulants
```

```@example ifo
hilb = FockSpace(:cavity)⊗FockSpace(:mirror)

a = Destroy(hilb,:a,1)
b = Destroy(hilb,:b,2)

@cnumbers Δ g κ 

H = Δ*b'*b + g*(a*b'+a'*b)

L = [sqrt(κ)*a]
nothing #hide
```

```@example ifo
ops = [a,a',b,b']
heisenberg_eqs = meanfield(ops,H,L;rates=ones(size(L)),order=3)

print(string(heisenberg_eqs))
```

```@example ifo
c = CorrelationFunction(a', a, heisenberg_eqs)

eqs = complete(heisenberg_eqs)
```

```@example ifo
using ModelingToolkit
using OrdinaryDiffEq
# Numerical solution
ps = (Δ, g, κ)
@named sys = ODESystem(eqs)
u0 = zeros(ComplexF64, length(eqs))
p0 = (1.0, 1.5, 0.25)
prob = ODEProblem(sys,u0,(0.0,10.0),ps.=>p0)
sol = solve(prob,RK4())
```
