using Documenter, SystemLevelHamiltonian

push!(LOAD_PATH,"../src/")

makedocs(sitename="SystemLevelHamiltonian", remotes = nothing)

deploydocs(
    repo = "github.com/jeffwack/SystemLevelHamiltonian.jl.git",
)
