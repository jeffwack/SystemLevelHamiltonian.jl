# SystemLevelHamiltonian.jl

SystemLevelHamiltonian.jl is a Julia package for creating and composing open
quantum systems using the SLH framework. 

## Quick Start

```@example quick
using SystemLevelHamiltonian
using SecondQuantizedAlgebra

# Create a Hilbert space and operators
hilb = FockSpace(:cavity)

@qnumbers a::Destroy(hilb)
@cnumbers ω κ

# Define system components 
H = ω * a' * a
L = [sqrt(κ) * a]
S = [1]

# Create SLH systems
cavityA = SLH(:A, [:in], [:out], S, L, H)
cavityB = SLH(:B, [:in], [:out], S, L, H)

sys = concatenate([cavityA, cavityB],:chain)
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
