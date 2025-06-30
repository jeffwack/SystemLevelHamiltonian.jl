using Documenter
using SystemLevelHamiltonian

push!(LOAD_PATH,"../src/")

pages = ["Introduction" => "index.md",
         "General readout" => "cascadedoutputfilters.md",
         "Linear IFO" => "linearsystems.md",
         "Examples" => ["freqdepsqz.md", "jcfisher.md","oameyeQFI.md"],
         "API" => "api.md"] 

makedocs(sitename="SystemLevelHamiltonian", 
        pages = pages,
        remotes = nothing)  

deploydocs(
    repo = "github.com/jeffwack/SystemLevelHamiltonian.jl.git"
)
