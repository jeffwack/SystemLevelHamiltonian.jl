using SystemLevelHamiltonian
using Documenter

DocMeta.setdocmeta!(SystemLevelHamiltonian, :DocTestSetup, :(using SystemLevelHamiltonian); recursive=true)

makedocs(;
    modules=[SystemLevelHamiltonian],
    authors="Jeffrey Wack <jeffwack111@gmail.com> and contributors",
    sitename="SystemLevelHamiltonian.jl",
    format=Documenter.HTML(;
        canonical="https://jeffwack.github.io/SystemLevelHamiltonian.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/jeffwack/SystemLevelHamiltonian.jl",
    devbranch="main",
)
