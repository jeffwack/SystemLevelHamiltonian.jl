"""
    cavity(name)

Create a basic optical cavity SLH system.

Creates a single-mode cavity with detuning and decay. The cavity has one input 
and one output port with direct transmission (S=1).

# Arguments
- `name`: Symbol identifying the cavity (used for operator and parameter naming)

# Returns
- `SLH`: System with Hamiltonian H = Δ·a†a and coupling L = [√κ·a]

# Parameters
- `κ`: Cavity decay rate
- `Δ`: Cavity detuning from driving field
"""
function cavity(name)
    #Define the filter cavity 
    hilb = FockSpace(name)
    a = Destroy(hilb,Symbol(:a,:_,name))

    (κ, Δ) = cnumbers(Symbol(:κ,:_,name),Symbol(:Δ,:_,name))

    return SLH(name,
                [:In],
                [:Out],
                [1],
                [sqrt(κ)*a],
                Δ*adjoint(a)*a)
end

"""
    qed_cavity(name)

Create a cavity QED system with a two-level atom.

Creates a cavity containing a two-level atom with Jaynes-Cummings coupling
and an external driving field.

# Arguments  
- `name`: Symbol identifying the system (used for operator and parameter naming)

# Returns
- `SLH`: System with cavity-atom interaction and driving

# Parameters
- `Δ`: Cavity detuning
- `g`: Atom-cavity coupling strength  
- `κ`: Cavity decay rate
- `h`: External driving amplitude
"""
function qed_cavity(name)
    hcav = FockSpace(Symbol(:cav,:_,name))
    hspin = NLevelSpace(Symbol(:cav,:_,name),2)
    hilb = SecondQuantizedAlgebra.tensor(hcav,hspin)

    # Define the operators
    a = Destroy(hilb,Symbol(:a,:_,name))

    σ(i,j) = Transition(hilb, Symbol(:σ,:_,name), i, j, 2)

    (Δ, g, κ, h,) = cnumbers(Symbol(:Δ,:_,name),Symbol(:g,:_,name),Symbol(:κ,:_,name),Symbol(:h,:_,name))

    # Hamiltonian
    H = Δ*a'*a + g*(a'*σ(1,2) + a*σ(2,1)) + h*(a + a')
    L = [sqrt(κ)*a]

    return SLH(name,
                [:In],
                [:Out],
                [1],
                L,
                H)
end


"""
    squeezing_cavity(name)

Create a squeezing cavity SLH system.

Creates a cavity that generates squeezed light through a parametric 
interaction (two-mode squeezing Hamiltonian).

# Arguments
- `name`: Symbol identifying the cavity (used for operator and parameter naming)

# Returns  
- `SLH`: System with squeezing Hamiltonian H = iϵ(a†² - a²) and coupling L = [√κ·a]

# Parameters
- `κ`: Cavity decay rate
- `ϵ`: Squeezing strength 
"""
function squeezing_cavity(name)
    hilb = FockSpace(:squeezer)
    a = Destroy(hilb, Symbol(:a,:_,name))

    (κ,ϵ) = cnumbers(Symbol(:κ,:_,name),Symbol(:ϵ,:_,name))

    return SLH(name,
                [:In],
                [:Out],
                [1],
                [sqrt(κ)*a],
                1im*ϵ*(adjoint(a)^2- a^2))
end

"""
    radiation_pressure_cavity(name)

Create an optomechanical cavity with radiation pressure coupling.

Creates a two-mode system with an optical cavity mode coupled to a 
mechanical oscillator through radiation pressure.

# Arguments
- `name`: Symbol identifying the system (used for operator and parameter naming)

# Returns
- `SLH`: Optomechanical system with optical and mechanical modes

# Parameters  
- `κ`: Cavity decay rate
- `Δ`: Cavity detuning
- `Ω`: Mechanical frequency
- `g`: Optomechanical coupling strength
"""
function radiation_pressure_cavity(name)
    hilb = SecondQuantizedAlgebra.tensor(FockSpace(Symbol(name,:_,:optical)), FockSpace(Symbol(name,:_,:mechanical)))
    a = Destroy(hilb,Symbol(:a,:_,name),1) 
    b = Destroy(hilb,Symbol(:b,:_,name),2) 

    (κ, Δ, Ω, g) = cnumbers(Symbol(:κ,:_,name),Symbol(:Δ,:_,name),Symbol(:Ω,:_,name),Symbol(:g,:_,name))

    return SLH(name,
                [:In],
                [:Out],
                [1],
                [sqrt(κ)*a],
                Δ*a'*a+Ω*b'*b - g*a'*a*(b'+b))
end
