using FDTDSolver
using Documenter

DocMeta.setdocmeta!(FDTDSolver, :DocTestSetup, :(using FDTDSolver); recursive=true)

makedocs(;
    modules=[FDTDSolver],
    authors="Mohamed Kamal AbdElrahman",
    repo="https://github.com/MKAbdElrahman/FDTDSolver.jl/blob/{commit}{path}#{line}",
    sitename="FDTDSolver.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MKAbdElrahman.github.io/FDTDSolver.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MKAbdElrahman/FDTDSolver.jl",
    devbranch="main",
)
