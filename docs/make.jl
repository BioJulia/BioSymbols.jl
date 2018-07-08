using Documenter, BioSymbols

makedocs(
    format = :html,
    sitename = "BioSymbols",
    pages = [
        "Home" => "index.md",
        "Nucleic Acids" => "nucleicacids.md",
        "Amino Acids" => "aminoacids.md",
        "Sequences" => "sequences.md",
        "References" => "references.md"
    ],
    authors = "Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/BioSymbols.jl.git",
    julia = "0.7",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
