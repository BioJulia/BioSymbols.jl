using Documenter, BioSymbols

makedocs()
deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-biojulia"),
    repo = "github.com/BioJulia/BioSymbols.jl.git",
    julia = "0.5",
    osname = "linux",
)
