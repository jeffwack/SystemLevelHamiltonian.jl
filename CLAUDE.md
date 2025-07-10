# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Architecture Overview

SystemLevelHamiltonian.jl is a focused Julia package providing core SLH (System-Level Hamiltonian) framework functionality for creating and combining open quantum systems. After recent refactoring, the package now provides minimal, reusable components for quantum system modeling.

### Core Components

- **SLH Framework**: The main abstraction for open quantum systems, defined in `src/slh_core.jl`. Each SLH system has inputs, outputs, scattering matrix (S), coupling vector (L), and Hamiltonian (H).

- **Quantum Symbol Utilities**: Defined in `src/qsymbols_core.jl`, provides utilities for extracting quantum operators and symbolic parameters from expressions.

- **Component Library**: Pre-built quantum components in `src/componentlibrary.jl` including basic cavities, squeezing cavities, radiation pressure cavities, and QED cavities.

### Key Dependencies

- **SecondQuantizedAlgebra.jl**: For quantum operators, Hilbert spaces, and symbolic algebra
- **Symbolics.jl**: For symbolic mathematics
- **LinearAlgebra.jl**: For matrix operations

### System Operations

Two primary ways to combine SLH systems:
1. `concatenate(name, systems_list)`: Parallel composition of systems
2. `feedbackreduce(system, output, input)`: Feedback connections between ports

### Package Structure

```
SystemLevelHamiltonian/
├── src/
│   ├── SystemLevelHamiltonian.jl  # Main module file
│   ├── slh_core.jl               # SLH struct and operations
│   ├── qsymbols_core.jl          # Quantum symbol utilities
│   └── componentlibrary.jl       # Pre-built components
├── test/
│   └── runtests.jl               # Test suite
├── docs/
│   ├── make.jl                   # Documentation builder
│   ├── Project.toml              # Doc dependencies
│   └── src/
│       ├── index.md              # Main documentation
│       └── api.md                # API reference
└── CLAUDE.md                     # This file
```

### Documentation

Documentation is built using Documenter.jl and includes:
- Package overview and quick start guide
- Complete API reference with docstrings
- Examples of system composition and feedback

### Testing

Test suite in `test/runtests.jl` covers:
- Quantum symbol utilities (`get_qsymbols`, `get_numsymbols`, `get_additive_terms`)
- SLH struct operations (`operators`, `parameters`)
- System composition (`concatenate`)
- Component library functions

## Development Guidelines

When adding or modifying package functionality, ensure consistency across these areas:

### 1. Function Implementation
- Add the function to the appropriate file (`slh_core.jl`, `qsymbols_core.jl`, or `componentlibrary.jl`)
- Follow Julia naming conventions (lowercase, underscores for multi-word)
- Include proper type annotations where helpful

### 2. Docstrings
- Add comprehensive docstrings following the template below
- Include function signature, description, arguments, and returns
- Keep examples simple and focused (avoid complex doctests)

```julia
"""
    function_name(arg1, arg2)

Brief description of what the function does.

More detailed explanation if needed.

# Arguments
- `arg1`: Description of first argument
- `arg2`: Description of second argument

# Returns
- Description of return value and type
"""
```

### 3. Exports
- Add new public functions to the appropriate export list in `src/SystemLevelHamiltonian.jl`
- Group exports logically (SLH operations, utilities, components)

### 4. Tests
- Add tests for new functionality in `test/runtests.jl`
- Use `@testset` to group related tests
- Test both expected behavior and edge cases
- Ensure tests are self-contained and don't rely on external state

### 5. Documentation
- Update API reference in `docs/src/api.md` to include new functions in appropriate `@docs` blocks
- Add examples to `docs/src/index.md` if the functionality warrants it
- Consider whether new functionality needs its own documentation page

### 6. Consistency Checklist

When adding/changing functionality, verify:
- [ ] Function is implemented with proper docstring
- [ ] Function is exported in `SystemLevelHamiltonian.jl`
- [ ] Function is included in `@docs` block in `docs/src/api.md`
- [ ] Tests are added to `test/runtests.jl`
- [ ] Documentation builds without errors (`julia --project=docs docs/make.jl`)
- [ ] All tests pass (`julia --project=. -e "using Pkg; Pkg.test()"`)

### Notes on Removed Functionality

The package was recently refactored to focus on core SLH functionality:
- Removed redundant `hilbert()` function (use `SecondQuantizedAlgebra.hilbert()` directly)
- Removed `check_hilberts()` function (SecondQuantizedAlgebra handles this automatically)
- Examples and complex applications moved to separate Breadboard package
- Fixed `concatenate()` function to properly flatten coupling vectors

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
