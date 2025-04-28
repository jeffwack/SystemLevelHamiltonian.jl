using Symbolics
using LinearAlgebra
using GLMakie

@variables γ Δ g ω

A = [-γ  0  Δ  0; #q_a
      0  0  0  ω; #q_b
      -Δ 2g -γ 0; #p_a
      2g -ω 0  0] #p_b

#    q_in        p_in      GW
B = [-sqrt(2*γ) 0          0;
     0          0          0;
     0          -sqrt(2*γ) 0;
     0          0          1]

C = [sqrt(2*γ) 0 0         0;
     0         0 sqrt(2*γ) 0]

D = [1 0 0; #q_out
     0 1 0] #p_out

@variables s

G =  inv(s*I - A)
T = C*G*B + D

#now we want to extract the transfer function from the GW input to the two outputs
Q = Symbolics.build_function(T[1,3],s, γ, Δ, g, ω,expression = Val{false})
P = Symbolics.build_function(T[2,3],s, γ, Δ, g, ω,expression = Val{false})

γ_val = 100*2*pi
Δ_val = 10
g_val = 200
ω_val = 100*2*pi

freq = logrange(1, 1000, length=1000)
s = 1.0im*freq

Qtf = abs.(Q.(s,γ_val,Δ_val,g_val,ω_val))
Ptf = abs.(P.(s,γ_val,Δ_val,g_val,ω_val))

fig = Figure()
ax = Axis(fig[1,1],yscale = log,xscale = log)

lines!(ax,freq, Qtf)
lines!(ax, freq, Ptf)

fig