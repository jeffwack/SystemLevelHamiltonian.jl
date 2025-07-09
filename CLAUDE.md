# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Architecture Overview

SystemLevelHamiltonian.jl is a Julia package for creating and combining open quantum systems using the SLH (System-Level Hamiltonian) framework. The package provides tools for quantum optics simulations and quantum information processing.

### Core Components

- **SLH Framework**: The main abstraction for open quantum systems, defined in `src/slh.jl`. Each SLH system has inputs, outputs, scattering matrix (S), coupling vector (L), and Hamiltonian (H). 

- **Quantum Symbols**: Defined in `src/qsymbols.jl`, provides utilities for working with quantum operators and symbolic expressions defined by QuantumCumulants.jl

- **Component Library**: Pre-built quantum components in `src/componentlibrary.jl` including cavities, squeezing cavities, radiation pressure cavities, and QED cavities.

- **Fisher Information**: Quantum Fisher Information calculations in `src/fisherinfo.jl` for parameter estimation and metrology applications.

- **Linear Systems**: State-space representation and transfer functions in `src/linearsystems.jl` for control theory applications.

### Key Dependencies

- **QuantumCumulants**: For symbolic quantum calculations
- **QuantumToolbox**: For simulating quantum system dynamics with a master
  equation solver.
- **ModelingToolkit**: For symbolic modeling
- **Symbolics**: For symbolic mathematics
- **OrdinaryDiffEq**: For differential equation solving

### System Operations

Two primary ways to combine SLH systems:
1. `concatenate()`: Series connection of systems
2. `feedbackreduce()`: Feedback connections between systems

### Examples Structure

The `examples/` directory contains demonstrations organized by physics:
- `Atom_IFO/`: Atom-interferometer systems
- `Jaynes_Cummings/`: Jaynes-Cummings model implementations
- `empty_cavity/`: Cavity quantum electrodynamics
- `lambda_in_cavity/`: Photon storage systems
- `linear_optomechanical/`: Optomechanical systems
- `single_spin/`: Single spin systems

### Documentation

Documentation is built using Documenter.jl and includes:
- API documentation
- Examples for IFO models (both QuantumCumulants and QuantumToolbox implementations)
- Frequency-dependent squeezing examples
- Quantum Fisher Information examples
- SLH to ABCD conversion examples

#### Documentation Structure
```
docs/
├── make.jl              # Documenter.jl configuration
├── Project.toml         # Documentation dependencies
├── Manifest.toml        # Locked dependencies
├── build/               # Generated documentation (gitignored)
└── src/                 # Source markdown files
    ├── index.md         # Main documentation page
    ├── api.md           # API documentation
    ├── cascadedoutputfilters.md
    ├── freqdepsqz.md
    ├── ifoQC.md
    ├── ifoQT.md
    ├── jcfisher.md
    ├── oameyeQFI.md
    └── slh_to_abcd.md   # SLH to ABCD conversion examples
```

#### Adding Documentation Pages
When adding new documentation pages:
1. Create the markdown file in `docs/src/`
2. Use `@example` blocks for executable Julia code
3. Update `docs/make.jl` to include the new page in the `pages` array
4. The `pages` array controls both the navigation structure and build order

### Testing

Minimal test suite in `test/runtests.jl`. The package is in development with version 1.0.0-DEV.

## Julia Best Practices

This section outlines key Julia programming best practices based on the Julia manual style guide, tailored for scientific computing package development.

### Code Organization

- **Write functions, not scripts**: Functions are more reusable, testable, and clarify what steps are being done
- **Functions should take arguments**: Avoid operating on global variables
- **Design for reusability**: Aim for generic, flexible function designs that work with multiple types

### Naming Conventions

- **Modules and types**: Use capitalized camel case (`SparseArrays`, `UnitRange`, `SystemLevelHamiltonian`)
- **Functions**: Use lowercase (`maximum`, `haskey`, `concatenate`)
- **Multi-word functions**: Use underscores when needed (`feedback_reduce`)
- **Mutating functions**: Append `!` to functions that modify arguments (`push!`, `pop!`)

### Type Design and Performance

- **Avoid overly-specific types**: Prefer abstract types like `Integer` over concrete types like `Int32`
- **Let Julia specialize**: The compiler will automatically optimize generic code
- **Handle type conversions at caller level**: Don't force conversions inside functions
- **Use `Int` literals in generic code**: When writing generic numerical algorithms

### Function Design Principles

Follow Base library argument ordering:
1. Function argument
2. IO stream  
3. Input being mutated
4. Type
5. Immutable inputs

### Best Practices for Scientific Computing

- **Prefer exported methods**: Don't directly access fields of types you don't own
- **Avoid type piracy**: Don't extend types from other packages inappropriately
- **Use multiple dispatch**: Design flexible implementations that work with different types
- **Minimize type conversions**: Let Julia's type system work efficiently
- **Use proper type checking**: Use `isa` and `<:` for type checking, not `==`

### Performance Guidelines

- **Create type-generic numerical algorithms**: Design functions to work with multiple numeric types
- **Let the compiler optimize**: Trust Julia's compiler to specialize generic code
- **Avoid unnecessary static type parameters**: Don't over-constrain your functions
- **Don't expose unsafe operations**: Keep unsafe operations internal to your implementation

## Julia Documentation Best Practices

This section outlines best practices for writing clear, comprehensive documentation in Julia, based on the Julia manual and Documenter.jl guidelines.

### Docstring Structure

Every function should have a docstring following this template:

```julia
"""
    function_name(arg1[, arg2])

Concise imperative description of function purpose.

Optional detailed explanation of function behavior.

# Arguments
- `arg1`: Description of first argument
- `arg2`: Description of second argument (if needed)

# Examples
```jldoctest
julia> example_usage()
expected_output
```

See also: [`related_function1`](@ref), [`related_function2`](@ref)
"""
```

### Docstring Guidelines

- **Function signature**: Always show with 4-space indent at the top
- **One-line description**: Use imperative form ("Create", "Calculate", "Return")
- **Avoid redundant language**: Don't use "The function..." or similar
- **Use Markdown formatting**: Backticks for code, proper headers
- **Line length**: Respect 92-character limit for readability
- **Self-contained examples**: Include executable doctests that work independently

### Documentation Organization

#### Package Structure
- Create `docs/` directory in package root
- Use `docs/src/` for markdown files
- Include `docs/make.jl` for Documenter.jl configuration
- Use `DocumenterTools.generate()` to set up initial structure

#### Content Organization
- **index.md**: Main documentation page with overview
- **API sections**: Use `@docs` blocks to include function docstrings
- **Navigation**: Use `@contents` and `@index` blocks for automatic navigation
- **Cross-references**: Use `[@ref]` syntax to link between sections

### Doctests and Examples

- Use `# Examples` section with executable doctests
- Format as ` ```jldoctest ` code blocks
- Make examples mimic REPL interactions
- Avoid `rand()` for consistent test outputs
- Include realistic use cases from the domain

### Advanced Documentation Techniques

- **Generic functions**: Document both generic functions and specific methods when behavior differs
- **Type documentation**: Use `@doc` macro for flexible documentation attachment
- **Dynamic documentation**: Consider for complex parametric types
- **Module context**: Use `@meta` block to specify current module

### Documentation Configuration

Using Documenter.jl:
- Configure `makedocs()` with `sitename` and `modules`
- Control sidebar navigation with `pages` argument
- Support both light and dark theme logos in `src/assets/`
- Use LiveServer for local preview during development

### Scientific Computing Specific Guidelines

- **Mathematical notation**: Use LaTeX math mode when appropriate
- **Physical units**: Always specify units in physics/engineering contexts
- **Parameter ranges**: Document valid input ranges and constraints
- **References**: Include citations to relevant papers or algorithms
- **Performance notes**: Document computational complexity when relevant
