# printing.jl
# ===========
#
# Print functions for NucleicAcids.
#
# This file is a part of the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/NucleicAcids.jl/blob/master/LICENSE.md

# NucleicAcid
# ===========

function Base.show(io::IO, nt::DNA)
    if isvalid(nt)
        if nt == DNA_Gap
            write(io, "DNA_Gap")
        else
            write(io, "DNA_", Char(nt))
        end
    else
        write(io, "Invalid DNA")
    end
    return
end

function Base.show(io::IO, nt::RNA)
    if isvalid(nt)
        if nt == RNA_Gap
            write(io, "RNA_Gap")
        else
            write(io, "RNA_", Char(nt))
        end
    else
        write(io, "Invalid RNA")
    end
    return
end

function Base.print(io::IO, nt::NucleicAcid)
    if !isvalid(nt)
        throw(ArgumentError("nucleic acid is invalid"))
    end
    write(io, Char(nt))
    return
end

# AminoAcid
# =========

function Base.show(io::IO, aa::AminoAcid)
    if isvalid(aa)
        if aa == AA_Term
            write(io, "AA_Term")
        elseif aa == AA_Gap
            write(io, "AA_Gap")
        else
            write(io, "AA_", Char(aa))
        end
    else
        write(io, "Invalid Amino Acid")
    end
    return
end

function Base.print(io::IO, aa::AminoAcid)
    if !isvalid(aa)
        throw(ArgumentError("invalid amino acid"))
    end
    write(io, Char(aa))
    return
end
