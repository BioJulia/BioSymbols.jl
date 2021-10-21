# BioSymbols.jl
# =============
#
# Types representing Biological Symbols.
#
# This file is a part of the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/NucleicAcids.jl/blob/master/LICENSE.md

__precompile__()

"""
A package to define types for amino acids and nucleic acids and to perform operations on these types.

Amino Acid Symbols:
- `AA_A`: Alanine 
- `AA_B`: Aspartic Acid or Asparagine 
- `AA_C`: Cysteine 
- `AA_D`: Aspartic Acid 
- `AA_E`: Glutamic Acid 
- `AA_F`: Phenylalanine 
- `AA_G`: Glycine 
- `AA_Gap`: Amino Acid Gap 
- `AA_H`: Histidine 
- `AA_I`: Isoleucine 
- `AA_J`: Leucine or Isoleucine 
- `AA_K`: Lysine 
- `AA_L`: Leucine 
- `AA_M`: Methionine 
- `AA_N`: Asparagine 
- `AA_O`: Pyrrolysine 
- `AA_P`: Proline 
- `AA_Q`: Glutamine 
- `AA_R`: Arginine 
- `AA_S`: Serine 
- `AA_T`: Threonine 
- `AA_Term`: Terminal 
- `AA_U`: Selenocysteine 
- `AA_V`: Valine 
- `AA_W`: Tryptophan 
- `AA_X`: Unspecified or Unknown Amino Acid 
- `AA_Y`: Tyrosine 
- `AA_Z`: Glutamine or Glutamic Acid 

- `ACGT`: Unambiguous DNA.
- `ACGTN`: Unambiguous DNA and `DNA_N`.
- `ACGU`: Unambiguous RNA.
- `ACGUN`: Unambiguous RNA and `RNA_N`.
- `AminoAcid`: An amino acid type.
- `BioSymbol`: The BioSymbol type is an abstract type that represents any kind of biological symbol that may appear in data such as biological sequences, SNP datasets and more.

DNA Symbols:
- `DNA`: A deoxyribonucleic acid type.
- `DNA_A`: DNA Adenine 
- `DNA_B`: DNA Cytosine, Guanine or Thymine 
- `DNA_C`: DNA Cytosine 
- `DNA_D`: DNA Adenine, Guanine or Thymine 
- `DNA_G`: DNA Guanine 
- `DNA_Gap`: DNA Gap 
- `DNA_H`: DNA Adenine, Cytosine or Thymine 
- `DNA_K`: DNA Guanine or Thymine 
- `DNA_M`: DNA Adenine or Cytosine 
- `DNA_N`: DNA Adenine, Cytosine, Guanine or Thymine 
- `DNA_R`: DNA Adenine or Guanine 
- `DNA_S`: DNA Cytosine or Guanine 
- `DNA_T`: DNA Thymine 
- `DNA_V`: DNA Adenine, Cytosine or Guanine 
- `DNA_W`: DNA Adenine or Thymine 
- `DNA_Y`: DNA Cytosine or Thymine 
- `NucleicAcid`: An abstract nucleic acid type.

RNA Symbols:
- `RNA`: A ribonucleic acid type.
- `RNA_A`: RNA Adenine 
- `RNA_B`: RNA Cytosine, Guanine or Uracil 
- `RNA_C`: RNA Cytosine 
- `RNA_D`: RNA Adenine, Guanine or Uracil 
- `RNA_G`: RNA Guanine 
- `RNA_Gap`: RNA Gap 
- `RNA_H`: RNA Adenine, Cytosine or Uracil 
- `RNA_K`: RNA Guanine or Uracil 
- `RNA_M`: RNA Adenine or Cytosine 
- `RNA_N`: RNA Adenine, Cytosine, Guanine or Uracil 
- `RNA_R`: RNA Adenine or Guanine 
- `RNA_S`: RNA Cytosine or Guanine 
- `RNA_U`: RNA Uracil 
- `RNA_V`: RNA Adenine, Cytosine or Guanine 
- `RNA_W`: RNA Adenine or Uracil 
- `RNA_Y`: RNA Cytosine or Uracil 

Functions:
- `alphabet`:
  + Get all symbols of `AminoAcid` in sorted order.
  + Get all symbols of `RNA` in sorted order.
  + Get all symbols of `DNA` in sorted order.

- `compatbits`:
  + Return the compatibility bits of `aa` as `UInt32`.
  + Return the compatibility bits of `nt` as `UInt8`.
  - `complement`: Return the complementary nucleotide of `nt`.

- `gap`:
  + Return `AA_Gap`.
  + Return `RNA_Gap`.
  + Return `DNA_Gap`.

- `isGC`: Test if `nt` is surely either guanine or cytosine.

- `isambiguous`:
  + Test if `aa` is an ambiguous amino acid.
  + Test if `nt` is an ambiguous nucleotide.

- `iscertain`:
  + Test if `aa` is a non-ambiguous amino acid.
  + Test if `nt` is a non-ambiguous nucleotide e.g.

- `iscompatible`: Test if `x` and `y` are compatible with each other.

- `isgap`: Test if `symbol` is a gap.

- `ispurine`: Test if `nt` is surely a purine.

- `ispyrimidine`: Test if `nt` is surely a pyrimidine.

- `BioSymbols.AA_INVALID`: Invalid Amino Acid 

- `BioSymbols.char_to_dna`: Lookup table used for converting characters to DNA symbol values The provided `convert` method should be used rather than this table, but you can use it if you insist and know what your are doing.

- `BioSymbols.char_to_rna`: Lookup table used for converting characters to RNA symbol values The provided `convert` method should be used rather than this table, but you can use it if you insist and know what your are doing.

- `BioSymbols.dna_to_char`: Lookup table for converting DNA symbol values to characters The provided `convert` method should be used rather than this table, but you can use it if you insist and know what your are doing.

- `BioSymbols.rna_to_char`: Lookup table for converting RNA symbol values to characters The provided `convert` method should be used rather than this table, but you can use it if you insist and know what your are doing.
"""
module BioSymbols

export
    # Types
    BioSymbol,
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
    alphabet,
    encoded_data,
    encode


"""
The BioSymbol type is an abstract type that represents
any kind of biological symbol that may appear in data
such as biological sequences, SNP datasets and more.

Every abstract or concrete subtype of BioSymbol is
expected to have the following methods defined:

`isterm`
`bytemask`
`prefix`
`type_text`

"""
abstract type BioSymbol end

encoded_data(x::BioSymbol) = reinterpret(encoded_data_eltype(typeof(x)), x)
function encode(::Type{T}, x) where T <: BioSymbol
    return reinterpret(T, convert(encoded_data_eltype(T), x))
end

# Enable broadcasting
Base.broadcastable(x::BioSymbol) = (x,)

include("nucleicacid.jl")
include("aminoacid.jl")

"""
    isgap(symbol::BioSymbol)

Test if `symbol` is a gap.
"""
isgap(symbol::BioSymbol) = symbol === gap(typeof(symbol))


# Arithmetic and Order
# --------------------

# These methods are necessary when deriving some algorithims
# like iteration, sort, comparison, and so on.
Base.isless(x::S, y::S) where S <: BioSymbol = isless(encoded_data(x), encoded_data(y))

@inline function Base.count_ones(symbol::BioSymbol)
    return count_ones(encoded_data(symbol))
end

@inline function Base.trailing_zeros(symbol::BioSymbol)
    return trailing_zeros(encoded_data(symbol))
end

# Printing BioSymbols
# -------------------

function suffix(symbol::BioSymbol)
    if isterm(symbol)
        return "_Term"
    end
    if isgap(symbol)
        return "_Gap"
    end
    return string('_', convert(Char, symbol))
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

Base.write(io::IO, symbol::BioSymbol) = write(io, encoded_data(symbol))
Base.read(io::IO, ::Type{T}) where T<:BioSymbol = encode(T, read(io, encoded_data_eltype(T)))

"""
    iscompatible(x::S, y::S) where S <: BioSymbol

Test if `x` and `y` are compatible with each other.

Examples
--------

```jldoctest
julia> iscompatible(AA_A, AA_R)
false

julia> iscompatible(AA_A, AA_X)
true

julia> iscompatible(DNA_A, DNA_A)
true

julia> iscompatible(DNA_C, DNA_N)  # DNA_N can be DNA_C
true

julia> iscompatible(DNA_C, DNA_R)  # DNA_R (A or G) cannot be DNA_C
false

```
"""
@inline function iscompatible(x::S, y::S) where S <: BioSymbol
    return compatbits(x) & compatbits(y) != 0
end

end
