using Documenter
using SystemLevelHamiltonian

makedocs(
    sitename = "SystemLevelHamiltonian.jl",
    format = Documenter.HTML(),
    modules = [SystemLevelHamiltonian],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)