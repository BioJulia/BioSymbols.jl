using Documenter, BiologicalSymbols

makedocs()
deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-biojulia"),
    repo = "github.com/BioJulia/BiologicalSymbols.jl.git",
    julia = "0.5",
    osname = "linux",
)
