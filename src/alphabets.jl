# Alphabet
# ========
#
# Alphabet of biological symbols.
#
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
# Alphabets of biological symbols.

`Alphabet` is perhaps the most important type trait for biological sequences in
BioSequences.jl.

An `Alphabet` represents a domain of biological symbols.

For example, `DNAAlphabet{2}` has a domain of unambiguous nucleotides
(i.e. A, C, G, and T).

Alphabet types restrict and define the set of biological symbols,
that can be encoded in a given biological sequence type.
They ALSO define *HOW* that encoding is done.

An `Alphabet` type defines the encoding of biological symbols with a pair
of associated `encoder` and `decoder` methods. These paired methods map
between biological symbol values and a binary representation of the symbol.

Any type A <: Alphabet, is expected to implement the `Base.eltype` method
for itself.
It is also expected to implement the `BitsPerSymbol` method.


"""
abstract type Alphabet end

"""
Alphabet of nucleic acids.
"""
abstract type NucleicAcidAlphabet <: Alphabet end

"""
DNA nucleotide alphabet.
"""
abstract type DNAAlphabet <: NucleicAcidAlphabet end
Base.eltype(::Type{A}) where A <: DNAAlphabet = DNA



struct AmbiguousDNA <: DNAAlphabet end
struct UnambiguousDNA <: DNAAlphabet end

"""
RNA nucleotide alphabet.
"""
abstract type RNAAlphabet <: NucleicAcidAlphabet end
Base.eltype(::Type{A}) where A <: RNAAlphabet = RNA

struct AmbiguousRNA <: RNAAlphabet end
struct UnambiguousRNA <: RNAAlphabet end

"""
Amino acid alphabet.
"""
struct AminoAcidAlphabet <: Alphabet end
Base.eltype(::Type{AminoAcidAlphabet}) = AminoAcid

"""
General character alphabet.
"""
struct CharAlphabet <: Alphabet end
Base.eltype(::Type{CharAlphabet}) = Char

"""
Void alphabet (internal use only).
"""
struct VoidAlphabet <: Alphabet end
Base.eltype(::Type{VoidAlphabet}) = Void




Base.eltype(::A) where A <: Alphabet = eltype(A)

@inline function Base.iterate(x::Union{UnambiguousDNA, UnambiguousRNA}, state = 0x01)
    state > 0x08 ? nothing : (reinterpret(eltype(x), state), state << 0x01)
end

@inline function Base.iterate(x::Union{AmbiguousDNA, AmbiguousRNA}, state = 0x00)
    state > 0x0F ? nothing : (reinterpret(eltype(x), state), state + 0x01)
end



symbols(::UnambiguousDNA) = ACGT
symbols(::UnambiguousRNA) = ACGU
symbols(::AmbiguousDNA) = alphabet(DNA)
symbols(::AmbiguousRNA) = alphabet(RNA)
symbols(::AminoAcidAlphabet) = alphabet(AminoAcid)
# TODO: this alphabet includes invalid Unicode scalar values
symbols(::CharAlphabet) = typemin(Char):typemax(Char)
symbols(::VoidAlphabet) = nothing


# Promotion of Alphabets
# ----------------------

for alph in (DNAAlphabet, RNAAlphabet)
    @eval function Base.promote_rule(::Type{A}, ::Type{B}) where {A<:$alph,B<:$alph}
        # TODO: Resolve this use of bits_per_symbol.
        return $alph{max(bits_per_symbol(A()),bits_per_symbol(B()))}
    end
end