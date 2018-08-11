using Documenter, BioSymbols

makedocs(
    format = :html,
    sitename = "BioSymbols",
    pages = [
        "Home" => "index.md",
        "Nucleic Acids" => "nucleicacids.md",
        "Amino Acids" => "aminoacids.md",
        "Alphabets" => "alphabets.md",
        "Sequences" => "sequences.md",
        "References" => "references.md"
    ],
    authors = "Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-material"),
    repo = "github.com/BioJulia/BioSymbols.jl.git",
    julia = "1.0",
    osname = "linux",
)
