# Linear Quantum System - QuantumToolbox

```@example ifo
using QuantumToolbox
```

```@example ifo
N = 4
a = destroy(N) ⊗ qeye(N)
b = qeye(N) ⊗ destroy(N)

H1 = QobjEvo(b'*b, (p,t)->p[1])
H2 = QobjEvo(a'*b+a*b', (p,t)->p[2])

H = H1+H2

L = [QobjEvo(a, (p,t)->p[3])]
nothing #hide
```

```@example ifo
tlist = LinRange(0, 100, 1000)
p0 = [0.01,500,10]
ψ0= basis(N,0) ⊗ basis(N,0)

eops = [a, a'*a,b,b'*b]

sol = mesolve(H, ψ0, tlist, L;params = p0,e_ops=eops)
nothing #hide
```

```@example ifo
using CairoMakie

fig = Figure()
ax = Axis(fig[1,1])

lines!(ax,tlist,real(sol.expect[2,:]),label = L"$\langle a' a \rangle$" )
lines!(ax,tlist,real(sol.expect[4,:]),label = L"$\langle b' b \rangle$" )
axislegend()

fig
```

```@example ifo
corr = correlation_2op_1t(H,ψ0,tlist,L,a,a';params = p0)
```

```@example ifo
(omega, S) = spectrum_correlation_fft(tlist,corr)

fig = Figure()
ax = Axis(fig[1,1],xscale = log10, yscale = log10)
nz = Int(ceil(length(omega)/2+1))

lines!(ax, omega[nz:end],S[nz:end])

fig
```

