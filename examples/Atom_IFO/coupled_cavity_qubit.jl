using QuantumOptics
using GLMakie

Ncutoff = 20 # Maximum photon number

fock = FockBasis(Ncutoff)
spin = SpinBasis(1//2)

a = destroy(fock)⊗identityoperator(fock)⊗identityoperator(spin)
at = create(fock)⊗identityoperator(fock)⊗identityoperator(spin)

b = identityoperator(fock)⊗destroy(fock)⊗identityoperator(spin)
bt = identityoperator(fock)⊗create(fock)⊗identityoperator(spin)

sm = identityoperator(fock)⊗identityoperator(fock)⊗sigmam(spin)
sp = identityoperator(fock)⊗identityoperator(fock)⊗sigmap(spin)
sz = identityoperator(fock)⊗identityoperator(fock)⊗sigmaz(spin)

Δ = 0.01
g = 3
h = 1
Ω = 1

H = Δ*sz + g*(b'*sm + b*sp) + Ω*(b'*a + b*a') + (a + a')*h

κ = 0.1
J = [sqrt(κ)*a, sqrt(κ)*b]

Ψ₀ = fockstate(fock,0) ⊗ fockstate(fock,0) ⊗ spindown(spin)
ρ₀ = Ψ₀ ⊗ dagger(Ψ₀)

N_steps = 1000
T = range(0,10,N_steps)
tout, ρt_master = timeevolution.master(T, ρ₀, H, J)

xvec = range(-5,5,length = 100)
yvec = range(-5,5,length = 100)

ρa = [ptrace(state,[2,3]) for state in ρt_master]
ρb = [ptrace(state,[1,3]) for state in ρt_master]
ρs = [ptrace(state,[1,2]) for state in ρt_master]
u = [real(rho[1]) for rho in ρs]
d = [real(rho[4]) for rho in ρs]

fig = Figure()

axa = Axis(fig[1, 1],aspect = 1, title = "signal injection cavity")
axb = Axis(fig[1, 2],aspect = 1,title = "spin interaction cavity")
axsu = Axis(fig[2, 1],aspect = 1,title = "ρ[1,1]",limits = (nothing, nothing, 0, 1))
axsd = Axis(fig[2, 2],aspect = 1,title = "ρ[2,2]",limits = (nothing, nothing, 0, 1))

supertitle = Label(fig[0, :], "Δ = $Δ, g = $g, Ω = $Ω, h = $h, κ = $κ", fontsize = 15)
# animation settings
framerate = 30

record(fig, "output.mp4", 1:N_steps;
        framerate = framerate) do ii
            Wa = wigner(ρa[ii], xvec, yvec)
            empty!(axa)
            heatmap!(axa,xvec,yvec,Wa)
            Wb = wigner(ρb[ii], xvec, yvec)
            empty!(axb)
            heatmap!(axb,xvec,yvec,Wb)

            empty!(axsu)
            lines!(axsu,tout[1:ii],u[1:ii],color = "blue")
            empty!(axsd)
            lines!(axsd,tout[1:ii],d[1:ii],color = "blue")
end
