# SystemLevelHamiltonian.jl

SystemLevelHamiltonian.jl is a Julia package for creating and composing open
quantum systems using the SLH framework. 

## Quick Start

```@example quick
using SystemLevelHamiltonian
using SecondQuantizedAlgebra

# Create a Hilbert space and operators
hilbA = FockSpace(:cavityA)
hilbB = FockSpace(:cavityB)

@qnumbers a::Destroy(hilbA) b::Destroy(hilbB)
@cnumbers ω1 ω2 κ1 κ2

# Define system components
H1 = ω1 * a' * a
L1 = [sqrt(κ1) * a]            
S1 = [1]                    

H2 = ω2 * b' * b
L2 = [sqrt(κ2) * b]
S2 = [1]

# Create SLH systems
cavityA = SLH(:A, [:in], [:out], S1, L1, H1)
cavityB = SLH(:B, [:in], [:out], S2, L2, H2)

sys = concatenate(:chain, [cavityA, cavityB])
sys = feedbackreduce(sys, :out_A, :in_B)

println(sys.H)

# Extract information
println(operators(sys))        # Get quantum operators
println(parameters(sys))    # Get symbolic parameters
```
## Overview of SLH systems

The SLH framework represents each open quantum systems with three components:
- **S**: Scattering matrix describing direct input-output coupling of external
  (bath) modes
- **L**: Coupling vector describing the interaction of the internal modes with
  the external modes 
- **H**: System Hamiltonian describing internal dynamics

### Component Library
The SLH framework enables you to create complicated quantum systems by combining
simple, reusable components
- Pre-built quantum components including:
  - Basic cavities
  - Squeezing cavities  
  - Radiation pressure cavities
  - Jaynes-Cummings QED cavity


## Dependencies

- [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) provides the symbolic algebra system for quantum operators

## References

- [The SLH framework for modeling quantum input-output networks](https://arxiv.org/pdf/1611.00375)
