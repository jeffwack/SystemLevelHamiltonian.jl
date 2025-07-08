using SystemLevelHamiltonian
using QuantumCumulants

# Example: Linear optomechanical system
# H = b†*b*ω + g*(a†*b + a*b†)  # Hamiltonian for linear optomechanical system  
# L = [√κ * a]                   # Optical cavity loss
# S = [1]                        # Trivial scattering

# Create Hilbert spaces
hilb = FockSpace(:optical) ⊗ FockSpace(:mechanical)

# Create operators
a = Destroy(hilb,:a,1)
b = Destroy(hilb,:b,2)

# Define parameters
@cnumbers ω κ g

# Build the system
H = b'*b*ω + g*(a'*b + a*b')
L = [√κ * a]
S = [1]

# Create SLH system
sys = SLH(:optomech, [:in], [:out], S, L, H)

println("System created:")
println("Name: ", sys.name)
println("Inputs: ", sys.inputs)
println("Outputs: ", sys.outputs)
println("Hamiltonian: ", sys.H)
println("Dissipation: ", sys.L)
println("Scattering: ", sys.S)

# Convert to ABCD representation
try
    (A, B, C, D) = SLH2ABCD(sys)
    
    println("\nABCD matrices:")
    println("A (drift): ", A)
    println("B (input coupling): ", B)
    println("C (output coupling): ", C)
    println("D (feedthrough): ", D)
catch e
    println("Error during conversion: ", e)
end
