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
    alphabet

import Automa
import Automa.RegExp: @re_str

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

Base.length(::BioSymbol) = 1
Base.iterate(sym::BioSymbol) = (sym, nothing)
Base.iterate(sym::BioSymbol, state) = nothing

include("nucleicacid.jl")
include("aminoacid.jl")

"""
    isgap(symbol::BioSymbol)

Test if `symbol` is a gap.
"""
isgap(symbol::BioSymbol) = symbol == gap(typeof(symbol))


# Arithmetic and Order
# --------------------

# These methods are necessary when deriving some algorithims
# like iteration, sort, comparison, and so on.

Base.:-(x::S, y::S) where S <: BioSymbol = convert(Int, x) - convert(Int, y)

Base.:+(x::Integer, y::BioSymbol) = y + x

Base.isless(x::S, y::S) where S <: BioSymbol = isless(convert(Integer, x), convert(Integer, y))

function Base.:~(x::BioSymbol)
    return reinterpret(typeof(x), ~reinterpret(UInt8, x) & bytemask(x))
end

function Base.:&(x::S, y::S) where S <: BioSymbol
    return reinterpret(S, convert(Integer, x) & convert(Integer, y))
end


@inline function Base.count_ones(symbol::BioSymbol)
    return count_ones(convert(Integer, symbol))
end

@inline function Base.trailing_zeros(symbol::BioSymbol)
    return trailing_zeros(convert(Integer, symbol))
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

Base.write(io::IO, symbol::BioSymbol) = write(io, convert(Integer, symbol))
Base.read(io::IO, ::Type{T}) where T<:BioSymbol = reinterpret(T, read(io, UInt8))    


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
