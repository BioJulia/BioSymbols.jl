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
For more information see:
- Amino Acid Symbols: https://github.com/BioJulia/BioSymbols.jl/blob/master/docs/src/aminoacids.md
- Nucleic Acid Symbols: https://github.com/BioJulia/BioSymbols.jl/blob/master/docs/src/nucleicacids.md
- Creating Sequences: https://github.com/BioJulia/BioSymbols.jl/blob/master/docs/src/sequences.md

Functions:
- `alphabet`
- `compatbits`
- `gap`
- `isGC`
- `isambiguous`
- `iscertain`
- `iscompatible`
- `isgap`
- `ispurine`
- `ispyrimidine`
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
    encode,
    stringbyte


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

# Less efficient fallback. Should only be called for symbols of AsciiAlphabet
"""
	stringbyte(::BioSymbol)::UInt8

For biosymbol types that can be represented as ASCII characters, `stringbyte(x)`
returns the printable ASCII byte that represents the character in a string.

# Examples
```julia
julia> stringbyte(DNA_A) == UInt8('A')
true

julia> stringbyte(AA_Gap) == UInt8('-')
true
```
"""
function stringbyte end

# Create a lookup table from biosymbol to the UInt8 for the character that would
# represent it in a string, e.g. DNA_G -> UInt8('G')
for alphabettype in ("DNA", "RNA", "AminoAcid")
    tablename = Symbol(uppercase(alphabettype), "_TO_BYTE")
    typ = Symbol(alphabettype)
    @eval begin
        const $(tablename) = let
            alph = alphabet($(typ))
            bytes = zeros(UInt8, length(alph))
            @inbounds for letter in alph
                bytes[reinterpret(UInt8, letter) + 1] = UInt8(Char(letter))
            end
            Tuple(bytes)
        end
        stringbyte(x::$(typ)) = @inbounds $(tablename)[reinterpret(UInt8, x) + 1]
    end
end

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
