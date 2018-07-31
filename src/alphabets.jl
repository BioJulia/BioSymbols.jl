# Alphabet
# ========
#
# Alphabet of biological symbols.
#
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

abstract type Alphabet end

"""
Alphabet of nucleic acids.
"""
abstract type NucleicAcidAlphabet <: Alphabet end

"""
Abstract type representing any DNA nucleotide alphabet.
"""
abstract type DNAAlphabet <: NucleicAcidAlphabet end

struct AmbiguousDNA <: DNAAlphabet end
struct UnambiguousDNA <: DNAAlphabet end

"""
Abstract type representing any RNA nucleotide alphabet.
"""
abstract type RNAAlphabet <: NucleicAcidAlphabet end

struct AmbiguousRNA <: RNAAlphabet end
struct UnambiguousRNA <: RNAAlphabet end

"""
Amino acid alphabet.
"""
struct AminoAcidAlphabet <: Alphabet end


"""
General character alphabet.
"""
struct CharAlphabet <: Alphabet end


"""
Nothing alphabet (internal use only).
"""
struct NothingAlphabet <: Alphabet end

Base.eltype(::Type{A}) where A <: DNAAlphabet = DNA
Base.eltype(::Type{A}) where A <: RNAAlphabet = RNA
Base.eltype(::Type{AminoAcidAlphabet}) = AminoAcid
Base.eltype(::Type{CharAlphabet}) = Char
Base.eltype(::Type{VoidAlphabet}) = Nothing
Base.eltype(::A) where A <: Alphabet = eltype(A)

@inline function Base.iterate(x::Union{UnambiguousDNA, UnambiguousRNA}, state = 0x01)
    state > 0x08 ? nothing : (reinterpret(eltype(x), state), state << 0x01)
end

@inline function Base.iterate(x::Union{AmbiguousDNA, AmbiguousRNA}, state = 0x00)
    state > 0x0F ? nothing : (reinterpret(eltype(x), state), state + 0x01)
end

Base.length(x::Union{UnambiguousDNA, UnambiguousRNA}) = 4
Base.length(x::Union{AmbiguousDNA, AmbiguousRNA}) = 16

function Base.getindex(x::Union{UnambiguousDNA, UnambiguousRNA}, i::Int)
    return reinterpret(eltype(x), 0x01 << (i - 1))
end
function Base.getindex(x::Union{AmbiguousDNA, AmbiguousRNA}, i::Int)
    return reinterpret(eltype(x), i - 1)
end
firstindex(x::Alphabet) = 1
lastindex(x::Alphabet) = length(x)

symbols(::Type{UnambiguousDNA}) = ACGT
symbols(::Type{UnambiguousRNA}) = ACGU
symbols(::Type{AmbiguousDNA}) = alphabet(DNA)
symbols(::Type{AmbiguousRNA}) = alphabet(RNA)
symbols(::Type{AminoAcidAlphabet}) = alphabet(AminoAcid)
symbols(::Type{CharAlphabet}) = typemin(Char):typemax(Char)
symbols(::Type{VoidAlphabet}) = nothing
symbols(::A) where A <: Alphabet = eltype(A)

# Promotion of Alphabets
# ----------------------

for alph in (DNAAlphabet, RNAAlphabet)
    @eval function Base.promote_rule(::Type{A}, ::Type{B}) where {A<:$alph,B<:$alph}
        return length(symbols(A)) >= length(symbols(B)) ? A : B
    end
end