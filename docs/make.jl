using Documenter, BioSymbols

makedocs()
deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-material"),
    repo = "github.com/BioJulia/BioSymbols.jl.git",
    julia = "0.6",
    osname = "linux",
)
