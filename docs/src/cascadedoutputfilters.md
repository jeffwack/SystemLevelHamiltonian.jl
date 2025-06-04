# Interferometer readout 

We want to model a general interferometer readout protocol. We will model different temporal modes of the output signal by using a cascaded series of filters which collect light according to $\kappa_i(t)$, the `filter function.'

First, we will look at building the Hamiltonian of this cascaded system using the SLH formalism. For now we will model our nonlinear interferometer with the Jaynes-Cummings system. The package SystemLevelHamiltonian.jl implements the SLH composition rules, using the quantum algebra provided by the package QuantumCumulants.jl.

```@example readout
using SystemLevelHamiltonian 

ifo = qed_cavity(:ifo)

```

Now we create our chain of filters
```@example readout
filters = [cavity(symb) for symb in [:A,:B]]
```

We use the SLH composition rules to assemble our system.

```@example readout
push!(filters,ifo)
sys = concatenate(:sys,filters)
sys = feedbackreduce(sys,:Out_ifo,:In_A)
sys = feedbackreduce(sys,:Out_A,:In_B)
```

Even with two output filters, this Hamiltonian is a mess. It doesn't help that for now the symbolic system does not simplify or group terms nicely. We want to represent this Hamiltonian as a matrix to facilitate a density matrix calculation. Automating this procedure is WIP, so let's do something a bit manual to organize this Hamiltonian.

$H_0 = Δ_A a_A^\dagger a_A+Δ_B a_B^\dagger a_B +g_{ifo}(a_{ifo}^\dagger σ_{ifo12} + a_{ifo} σ_{ifo21})+Δ_{ifo} a_{ifo}^\dagger a_{ifo}$

$H_1 = h_{ifo} (a_{ifo} + a_{ifo}^\dagger)$

$H_2 = \frac{i}{2}\sqrt{κ_A κ_{ifo}} (a_A a_{ifo}^\dagger -a_A^\dagger a_{ifo})$

$H_3 = \frac{i}{2}\sqrt{κ_B κ_A}(a_A^\dagger a_B - a_A a_B^\dagger)$

$H_4 = \frac{i}{2}\sqrt{κ_B κ_{ifo}}(a_B a_{ifo}^\dagger - a_B^\dagger a_{ifo})$

Now, we will build this system and simulate its time evolution using the package QuantumToolbox. We will send $\Delta_A$ and $\Delta_B$ to zero, as well as make the coupling rates $\kappa_A$ and $\kappa_B$ functions of time.

```@example readout
using QuantumToolbox
```

We first need to establish an order on the subspaces of our Hilbert space, as well as a cutoff for our Fock states.

$\text{spin} \otimes \text{ifo} \otimes \text{A} \otimes \text{B}$

```@example readout
N = 4
σ  = sigmam() ⊗ qeye(N) ⊗ qeye(N) ⊗ qeye(N)
a_ifo = qeye(2) ⊗ destroy(N) ⊗ qeye(N) ⊗ qeye(N) 
a_A = qeye(2) ⊗ qeye(N) ⊗ destroy(N) ⊗ qeye(N) 
a_B = qeye(2) ⊗ qeye(N) ⊗ qeye(N) ⊗ destroy(N)
nothing # hide
```

We can now define our constants and build the first part of the Hamiltonian
```@example readout
g = 2
Δ = 0.5

H_0 = Δ*a_ifo'*a_ifo + g*(a_ifo'*σ + a_ifo*σ')
```

We now define the first time dependent part
```@example readout



function signal(p,t)
    return exp(-(t-p.d)^2/p.s^2)*p.A*sin(p.w*t)
end

H_1 = QobjEvo(a_ifo+a_ifo', signal)
nothing # hide
```

```@example readout
kappaifo = 0.3

function kappaA(p,t)
    return Complex(exp(-(t-p.dA)^2/p.sA^2)*p.AA*sin(p.wA*t))
end

function kappaB(p,t)
    return Complex(exp(-(t-p.dB)^2/p.sB^2)*p.AB*sin(p.wB*t))
end

function f2(p,t)
    return 1.0im/2*sqrt(kappaA(p,t)*kappaifo)
end
H_2 = QobjEvo(a_A*a_ifo' - a_A'*a_ifo,f2)

function f3(p,t)
    return 1.0im/2*sqrt(kappaB(p,t)*kappaA(p,t))
end
H_3 = QobjEvo(a_A'*a_B - a_A*a_B',f3)

function f4(p,t)
    return 1.0im/2*sqrt(kappaB(p,t)*kappaifo)
end
H_4 = QobjEvo(a_B*a_ifo' - a_B'*a_ifo,f4)

H_total = H_0 + H_1 + H_2 + H_3 + H_4
```

We now have the total Hamiltonian, which we can simulate by providing an initial state as well as parameters for all the time dependent functions.

```@example readout
ψ0 = basis(2,0) ⊗ fock(N, 0) ⊗ fock(N, 0) ⊗ fock(N, 0)
tlist = 0:0.1:100 

p = (
    d = 5, #center of signal pulse
    s = 2, #width of signal pulse
    A = 4, #amplitude of signal pulse
    w = 3*2*pi, #frequency of signal pulse 
    dA = 3, #center 
    sA = 2, #width 
    AA = 2, #amplitude 
    wA = 2*2*pi, #frequency 
    dB = 8, #center 
    sB = 2, #width 
    AB = 2, #amplitude 
    wB = 4*2*pi #frequency 
)

sol_me  = mesolve(H_total,  ψ0, tlist, params = p)
```
Now, we want to calculate the quantum Fisher information of the final state left in the filters with respect to the frequency of the signal (for example). To do this, we are going to wrap this solver call in a function which takes a frequency, time evolves, extracts the final state, and traces out the interferometer. This will then allow us to use the finite difference method to calculate the derivative with respect to a parameter and calculate the QFI.

```@example readout
function final_state(omega)
    p = (
        d = 5, #center of signal pulse
        s = 2, #width of signal pulse
        A = 4, #amplitude of signal pulse
        w = omega, #frequency of signal pulse 
        dA = 3, #center 
        sA = 2, #width 
        AA = 2, #amplitude 
        wA = 2*2*pi, #frequency 
        dB = 8, #center 
        sB = 2, #width 
        AB = 2, #amplitude 
        wB = 4*2*pi #frequency 
    )

    sol_me  = mesolve(H_total,  ψ0, tlist, params = p)
    rho = ptrace(sol_me.states[end],(3,4))
    return hermitian_data(rho)
end

final_state(3*2*pi)
```


```@example readout

function derivative(A,dA)
    rho1 = final_state(A)
    rho2 = final_state(A+dA)
    
    return (rho1, rho2, (rho2 - rho1)/dA)
end

(rho1, rho2, rhodot) = derivative(3*2*pi,0.001)

L = sld_operator(rho1,rhodot)

qfi = tr(L*L*rho1)

```

We can wrap all of this into a single function which will output the QFI of the system as a whole as well as each subsystem

```@example readout
function all_qfi(p,eps)

    sol_me = mesolve(H_total,  ψ0, tlist, params = p)
    rho_final = ket2dm(sol_me.states[end])

    pp= merge(p,(w=p.w+eps,))
    sol_me_prime = mesolve(H_total,  ψ0, tlist, params = pp)
    rho_final_prime = ket2dm(sol_me_prime.states[end])

    rho_dot_full = (hermitian_data(rho_final_prime) - hermitian_data(rho_final))/eps

    L_full = sld_operator(hermitian_data(rho_final),rho_dot_full)
    qfi_full = real(tr(L_full*L_full*hermitian_data(rho_final)))
    println("qfi_full = $qfi_full")

    rho_ifo = ptrace(rho_final,(1,2))
    rho_ifo_prime = ptrace(rho_final_prime,(1,2))
    rho_dot_ifo = (hermitian_data(rho_ifo_prime) - hermitian_data(rho_ifo))/eps
    L_ifo = sld_operator(hermitian_data(rho_ifo),rho_dot_ifo)
    qfi_ifo = real(tr(L_ifo*L_ifo*hermitian_data(rho_ifo)))
    println("qfi_ifo = $qfi_ifo")

    rho_filters = ptrace(rho_final,(3,4))
    rho_filters_prime = ptrace(rho_final_prime,(3,4))
    rho_dot_filters = (hermitian_data(rho_filters_prime) - hermitian_data(rho_filters))/eps
    L_filters = sld_operator(hermitian_data(rho_filters),rho_dot_filters)
    qfi_filters = real(tr(L_filters*L_filters*hermitian_data(rho_filters)))
    println("qfi_filters = $qfi_filters")

    rho_A = ptrace(rho_final,3)
    rho_A_prime = ptrace(rho_final_prime,3)
    rho_dot_A = (hermitian_data(rho_A_prime) - hermitian_data(rho_A))/eps
    L_A = sld_operator(hermitian_data(rho_A),rho_dot_A)
    qfi_A = real(tr(L_A*L_A*hermitian_data(rho_A)))
    println("qfi_A = $qfi_A")

    rho_B = ptrace(rho_final,4)
    rho_B_prime = ptrace(rho_final_prime,4)
    rho_dot_B = (hermitian_data(rho_B_prime) - hermitian_data(rho_B))/eps
    L_B = sld_operator(hermitian_data(rho_B),rho_dot_B)
    qfi_B = real(tr(L_B*L_B*hermitian_data(rho_B)))
    println("qfi_B = $qfi_B")


    return
end

eps = 0.0005
nothing #hide
```

## Nominal parameters
```@example readout
p_nominal = (
    d = 5, #center of signal pulse
    s = 2, #width of signal pulse
    A = 4, #amplitude of signal pulse
    w = 3*2*pi, #frequency of signal pulse 
    dA = 3, #center 
    sA = 2, #width 
    AA = 2, #amplitude 
    wA = 2*2*pi, #frequency 
    dB = 8, #center 
    sB = 2, #width 
    AB = 2, #amplitude 
    wB = 4*2*pi #frequency 
    )
all_qfi(p_nominal,eps)
```

## Sanity check - second filter off
```@example readout
p = merge(p_nominal,(AB=0,))
all_qfi(p,eps)
```

## sanity check - late signal
```@example readout
p = merge(p_nominal,(d=50,))
all_qfi(p,eps)
```

## increase signal amplitude
```@example readout
p = merge(p_nominal,(A=p_nominal.A*4,))
all_qfi(p,eps)
```


