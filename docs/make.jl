using Documenter, SystemLevelHamiltonian

push!(LOAD_PATH,"../src/")

makedocs()

deploydocs(
    repo = "github.com/jeffwack/SystemLevelHamiltonian.jl.git",
)