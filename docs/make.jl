using Documenter, BioSymbols

makedocs(
    format = Documenter.HTML(),
    sitename = "BioSymbols",
    pages = [
        "Home" => "index.md",
        "Nucleic Acids" => "nucleicacids.md",
        "Amino Acids" => "aminoacids.md",
        "Sequences" => "sequences.md",
        #"References" => "references.md",
        "Library" => [
            "Public" => "lib/public.md",
            hide("Internals" => "lib/internals.md")
        ]
    ],
    authors = "Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/BioSymbols.jl.git",
    deps = nothing,
    make = nothing
)
