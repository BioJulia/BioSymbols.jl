# BioSymbols.jl
# =============
#
# Types representing Biological Symbols.
#
# This file is a part of the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/NucleicAcids.jl/blob/master/LICENSE.md

__precompile__()

module BioSymbols

export
    # Types
    NucleicAcid,
    DNA,
    RNA,
    AminoAcid,

    # Constants
    DNA_A,
    DNA_C,
    DNA_G,
    DNA_T,
    DNA_M,
    DNA_R,
    DNA_W,
    DNA_S,
    DNA_Y,
    DNA_K,
    DNA_V,
    DNA_H,
    DNA_D,
    DNA_B,
    DNA_N,
    DNA_Gap,
    ACGT,
    ACGTN,
    RNA_A,
    RNA_C,
    RNA_G,
    RNA_U,
    RNA_M,
    RNA_R,
    RNA_W,
    RNA_S,
    RNA_Y,
    RNA_K,
    RNA_V,
    RNA_H,
    RNA_D,
    RNA_B,
    RNA_N,
    RNA_Gap,
    ACGU,
    ACGUN,
    AA_A,
    AA_R,
    AA_N,
    AA_D,
    AA_C,
    AA_Q,
    AA_E,
    AA_G,
    AA_H,
    AA_I,
    AA_L,
    AA_K,
    AA_M,
    AA_F,
    AA_P,
    AA_S,
    AA_T,
    AA_W,
    AA_Y,
    AA_V,
    AA_O,
    AA_U,
    AA_B,
    AA_J,
    AA_Z,
    AA_X,
    AA_Term,
    AA_Gap,

    # Bit Operations
    gap,
    isGC,
    ispurine,
    ispyrimidine,
    isambiguous,
    iscertain,
    isgap,
    complement,
    iscompatible,
    compatbits,
    alphabet

import Automa
import Automa.RegExp: @re_str

abstract type BioSymbol end

include("nucleicacid.jl")
include("aminoacid.jl")

"""
    isgap(aa::AminoAcid)

Test if `aa` is a gap.
"""
isgap(symbol::BioSymbol) = symbol == gap(typeof(symbol))

isterm(symbol::NucleicAcid) = false
isterm(symbol::AminoAcid) = symbol == AA_Term


# Printing BioSymbols
# -------------------

prefix(::DNA) = "DNA"
prefix(::RNA) = "RNA"
prefix(::AminoAcid) = "AA"
type_text(::AminoAcid) = "Amino Acid"
type_text(x::NucleicAcid) = prefix(x)

function suffix(symbol::BioSymbol)
    if isterm(symbol)
        return "_Term"
    end
    if isgap(symbol)
        return "_Gap"
    end
    return "_$(convert(Char, symbol))"
end

function Base.show(io::IO, symbol::BioSymbol)
    if isvalid(symbol)
        write(io, prefix(symbol), suffix(symbol))
    else
        write(io, "Invalid ", type_text(symbol))
    end
    return nothing
end

function Base.print(io::IO, symbol::BioSymbol)
    if !isvalid(symbol)
        throw(ArgumentError("Invalid $(type_text(symbol))"))
    end
    write(io, convert(Char, symbol))
    return nothing
end

end
