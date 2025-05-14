using Symbolics
using LinearAlgebra
using GLMakie

@variables ω Ω g κ

A = [-κ ω/2 0 0;
     -ω/2 -κ g 0;
     0 0 0 Ω/2;
     g 0 -Ω/2 0]

B = [-sqrt(2κ) 0 0;
     0 -sqrt(2κ) 0;
     0 0 0;
     0 0 1]

C = [sqrt(2κ) 0 0 0;
     0 sqrt(2κ) 0 0]

D = [1 0 0;
     0 1 0]

@variables s

G =  inv(s*I - A)
T = C*G*B + D

#now we want to extract the transfer function from the GW input to the two outputs
Q = Symbolics.build_function(T[2,3],s, ω, Ω, g, κ,expression = Val{false})
#P = Symbolics.build_function(T[2,3],s, ω, Ω, g, κ,expression = Val{false})
K = Symbolics.build_function(T[2,2],s, ω, Ω, g, κ,expression = Val{false})
R = Symbolics.build_function(T[2,1],s, ω, Ω, g, κ,expression = Val{false})

ω_val = 100*2*pi
Ω_val = 0.1*2*pi
g_val = 3
κ_val = 100

freq = logrange(1, 10000, length=1000)
s = 1.0im*freq

Qtf = Q.(s, ω_val, Ω_val, g_val, κ_val)
Ktf = K.(s,ω_val, Ω_val, g_val, κ_val)
Rtf = R.(s,ω_val, Ω_val, g_val, κ_val)

#what is going on here? Can't bring this home.
G = (Qtf .+ Ktf .+Rtf)./(Qtf) 

fig = Figure()
ax = Axis(fig[1,1],yscale = log,xscale = log)

#lines!(ax,freq, Qtf)
#lines!(ax, freq, Ktf)
#lines!(ax, freq, Rtf)
lines!(ax, freq, abs.(G))

fig