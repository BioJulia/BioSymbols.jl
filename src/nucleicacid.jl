# Nucleic Acid
# ============
#
# DNA and RNA types.
#
# This file is a part of the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/NucleicAcids.jl/blob/master/LICENSE.md

# NucleicAcid Encoding
# -------------------
#
# Unambiguous nucleotides are represented in one-hot encoding as follows:
#
#   | NucleicAcid | Bits |
#   | ----------- | ---- |
#   |     A       | 0001 |
#   |     C       | 0010 |
#   |     G       | 0100 |
#   |    T/U      | 1000 |
#
# Ambiguous nucleotides are bitwise OR of these four nucleotides. For example, R
# , A or G, is represented as 0101 (= A: 0001 | G: 0100). The gap symbol is
# always 0000.  The meaningful four bits are stored in the least significant
# bits of a byte.

"""
An abstract nucleic acid type.
"""
abstract type NucleicAcid <: BioSymbol end

"""
A deoxyribonucleic acid type.
"""
primitive type DNA <: NucleicAcid 8 end

"""
A ribonucleic acid type.
"""
primitive type RNA <: NucleicAcid 8 end

prefix(::DNA) = "DNA"
prefix(::RNA) = "RNA"
type_text(x::NucleicAcid) = prefix(x)
isterm(symbol::NucleicAcid) = false
encoded_data_eltype(::Type{DNA}) = UInt8
encoded_data_eltype(::Type{RNA}) = UInt8


# Conversion from/to compatible types
# ---------------------------

Base.convert(::Type{DNA}, nt::RNA) = reinterpret(DNA, nt)
DNA(nt::RNA) = convert(DNA, nt)
Base.convert(::Type{RNA}, nt::DNA) = reinterpret(RNA, nt)
RNA(nt::DNA) = convert(RNA, nt)


# Conversion from/to characters
# -----------------------------

function Base.convert(::Type{DNA}, c::Char)
    if c > '\uff'
        throw(InexactError(:convert, DNA, repr(c)))
    end
    @inbounds dna = char_to_dna[convert(Int, c) + 1]
    if !isvalid(DNA, dna)
        throw(InexactError(:convert, DNA, repr(c)))
    end
    return encode(DNA, dna)
end
DNA(c::Char) = convert(DNA, c)

function Base.convert(::Type{RNA}, c::Char)
    if c > '\uff'
        throw(InexactError(:convert, RNA, repr(c)))
    end
    @inbounds rna = char_to_rna[convert(Int, c) + 1]
    if !isvalid(RNA, rna)
        throw(InexactError(:convert, RNA, repr(c)))
    end
    return encode(RNA, rna)
end
RNA(c::Char) = convert(RNA, c)

function Base.convert(::Type{Char}, nt::DNA)
    return dna_to_char[encoded_data(nt) + 1]
end

Base.Char(nt::DNA) = convert(Char, nt)

function Base.convert(::Type{Char}, nt::RNA)
    return rna_to_char[encoded_data(nt) + 1]
end

Base.Char(nt::RNA) = convert(Char, nt)


# Encoding of DNA and RNA NucleicAcids
# ------------------------------------

"""
Lookup table used for converting characters to DNA symbol values

The provided `convert` method should be used rather than this table, but you can
use it if you insist and know what your are doing.

!!! note
    The array is indexed by converting a character to an integer. When indexed, it
    returns a UInt8 with the bit pattern on the corresponding nucleic acid.
    The `convert(DNA, x)` method does this for you.

!!! warning
    If you index this array with a character that is greater than '\uff', then
    you will get a bounds error. The `convert(DNA, x)` method checks such things
    to avoid this for you.

!!! warning
    If you index this array with a character that does not have a corresonding
    DNA symbol, then you get a byte with the bit pattern `0x80`, which is an
    invalid DNA symbol and will be of no use to you. The `convert(DNA, x)`
    checks such things for you and throws an exception gracefully if such a
    situation arises.
"""
const char_to_dna = [0x80 for _ in 0x00:0xff]

"""
Lookup table for converting DNA symbol values to characters

The provided `convert` method should be used rather than this table, but you can
use it if you insist and know what your are doing.

!!! note
    The array is indexed by reinterpreting a DNA symbol value as an UInt8.
    When indexed, it returns the character corresponding to the symbol.
    The `convert(Char, x::DNA)` method does this for you.

!!! warning
    If you index this array with an invalid DNA symbol, then you will hit a
    bounds error. If you construct DNA symbols properly, then this scenario
    should never occur.
"""
dna_to_char

# Derived from "The DDBJ/ENA/GenBank Feature Table Definition"
# §7.4.1 Nucleotide base code (IUPAC)
# http://www.insdc.org/documents/feature_table.html#7.4.1
const dna_to_char = let
    chararray = Vector{Char}(undef, 16)
    for (char, doc, bits) in [
        ('-', "DNA Gap",                                   0b0000),
        ('A', "DNA Adenine",                               0b0001),
        ('C', "DNA Cytosine",                              0b0010),
        ('G', "DNA Guanine",                               0b0100),
        ('T', "DNA Thymine",                               0b1000),
        ('M', "DNA Adenine or Cytosine",                   0b0011),
        ('R', "DNA Adenine or Guanine",                    0b0101),
        ('W', "DNA Adenine or Thymine",                    0b1001),
        ('S', "DNA Cytosine or Guanine",                   0b0110),
        ('Y', "DNA Cytosine or Thymine",                   0b1010),
        ('K', "DNA Guanine or Thymine",                    0b1100),
        ('V', "DNA Adenine, Cytosine or Guanine",          0b0111),
        ('H', "DNA Adenine, Cytosine or Thymine",          0b1011),
        ('D', "DNA Adenine, Guanine or Thymine",           0b1101),
        ('B', "DNA Cytosine, Guanine or Thymine",          0b1110),
        ('N', "DNA Adenine, Cytosine, Guanine or Thymine", 0b1111)]
        var = Symbol("DNA_", char != '-' ? char : "Gap")
        @eval begin
            @doc $(doc) const $(var) = encode(DNA, $(bits))
            char_to_dna[$(convert(Int, char) + 1)] = char_to_dna[$(convert(Int, lowercase(char)) + 1)] = $(bits)
            $(chararray)[$(convert(Int, bits) + 1)] = $(char)
        end
    end
    Tuple(chararray)
end

@eval function alphabet(::Type{DNA})
    return $(tuple([encode(DNA, x) for x in 0b0000:0b1111]...))
end

"""
    alphabet(DNA)

Get all symbols of `DNA` in sorted order.

Examples
--------

```jldoctest
julia> alphabet(DNA)
(DNA_Gap, DNA_A, DNA_C, DNA_M, DNA_G, DNA_R, DNA_S, DNA_V, DNA_T, DNA_W, DNA_Y, DNA_H, DNA_K, DNA_D, DNA_B, DNA_N)

julia> issorted(alphabet(DNA))
true

```
"""
alphabet(::Type{DNA})

"""
    ACGT

Unambiguous DNA.

Examples
--------

```jldoctest
julia> ACGT
(DNA_A, DNA_C, DNA_G, DNA_T)

```
"""
const ACGT = (DNA_A, DNA_C, DNA_G, DNA_T)

"""
    ACGTN

Unambiguous DNA and `DNA_N`.

Examples
--------

```jldoctest
julia> ACGTN
(DNA_A, DNA_C, DNA_G, DNA_T, DNA_N)

```
"""
const ACGTN = (DNA_A, DNA_C, DNA_G, DNA_T, DNA_N)

"""
Lookup table used for converting characters to RNA symbol values

The provided `convert` method should be used rather than this table, but you can
use it if you insist and know what your are doing.

!!! note
    The array is indexed by converting a character to an integer. When indexed, it
    returns a UInt8 with the bit pattern on the corresponding nucleic acid.
    The `convert(RNA, x)` method does this for you.

!!! warning
    If you index this array with a character that is greater than '\uff', then
    you will get a bounds error. The `convert(RNA, x)` method checks such things
    to avoid this for you.

!!! warning
    If you index this array with a character that does not have a corresonding
    RNA symbol, then you get a byte with the bit pattern `0x80`, which is an
    invalid RNA symbol and will be of no use to you. The `convert(RNA, x)`
    checks such things for you and throws an exception gracefully if such a
    situation arises.
"""
const char_to_rna = [0x80 for _ in 0x00:0xff]

"""
Lookup table for converting RNA symbol values to characters

The provided `convert` method should be used rather than this table, but you can
use it if you insist and know what your are doing.

!!! note
    The array is indexed by reinterpreting a RNA symbol value as an UInt8.
    When indexed, it returns the character corresponding to the symbol.
    The `convert(Char, x::RNA)` method does this for you.

!!! warning
    If you index this array with an invalid RNA symbol, then you will hit a
    bounds error. If you construct RNA symbols properly, then this scenario
    should never occur.
"""
rna_to_char


const rna_to_char = let
    chararray = Vector{Char}(undef, 16)
    for (char, doc, dna) in [
        ('-', "RNA Gap",                                  DNA_Gap),
        ('A', "RNA Adenine",                              DNA_A  ),
        ('C', "RNA Cytosine",                             DNA_C  ),
        ('G', "RNA Guanine",                              DNA_G  ),
        ('U', "RNA Uracil",                               DNA_T  ),
        ('M', "RNA Adenine or Cytosine",                  DNA_M  ),
        ('R', "RNA Adenine or Guanine",                   DNA_R  ),
        ('W', "RNA Adenine or Uracil",                    DNA_W  ),
        ('S', "RNA Cytosine or Guanine",                  DNA_S  ),
        ('Y', "RNA Cytosine or Uracil",                   DNA_Y  ),
        ('K', "RNA Guanine or Uracil",                    DNA_K  ),
        ('V', "RNA Adenine, Cytosine or Guanine",         DNA_V  ),
        ('H', "RNA Adenine, Cytosine or Uracil",          DNA_H  ),
        ('D', "RNA Adenine, Guanine or Uracil",           DNA_D  ),
        ('B', "RNA Cytosine, Guanine or Uracil",          DNA_B  ),
        ('N', "RNA Adenine, Cytosine, Guanine or Uracil", DNA_N  )]
        var = Symbol("RNA_", char != '-' ? char : "Gap")
        @eval begin
            @doc $(doc) const $(var) = reinterpret(RNA, $(dna))
            char_to_rna[$(convert(Int, char) + 1)] = char_to_rna[$(convert(Int, lowercase(char) + 1))] = reinterpret(UInt8, $(dna))
            $(chararray)[$(encoded_data(dna) + 1)] = $(char)
        end
    end
    Tuple(chararray)
end

@eval function alphabet(::Type{RNA})
    return $(tuple([encode(RNA, x) for x in 0b0000:0b1111]...))
end

"""
    alphabet(RNA)

Get all symbols of `RNA` in sorted order.

Examples
--------

```jldoctest
julia> alphabet(RNA)
(RNA_Gap, RNA_A, RNA_C, RNA_M, RNA_G, RNA_R, RNA_S, RNA_V, RNA_U, RNA_W, RNA_Y, RNA_H, RNA_K, RNA_D, RNA_B, RNA_N)

julia> issorted(alphabet(RNA))
true

```
"""
alphabet(::Type{RNA})

"""
    ACGU

Unambiguous RNA.

Examples
--------

```jldoctest
julia> ACGU
(RNA_A, RNA_C, RNA_G, RNA_U)

```
"""
const ACGU = (RNA_A, RNA_C, RNA_G, RNA_U)

"""
    ACGUN

Unambiguous RNA and `RNA_N`.

Examples
--------

```jldoctest
julia> ACGUN
(RNA_A, RNA_C, RNA_G, RNA_U, RNA_N)

```
"""
const ACGUN = (RNA_A, RNA_C, RNA_G, RNA_U, RNA_N)

gap(::Type{DNA}) = DNA_Gap
gap(::Type{RNA}) = RNA_Gap

"""
    isGC(nt::NucleicAcid)

Test if `nt` is surely either guanine or cytosine.
"""
function isGC(nt::NucleicAcid)
    bits = encoded_data(nt)
    return bits != 0 && (bits & 0b1001) == 0
end

"""
    ispurine(nt::NucleicAcid)

Test if `nt` is surely a purine.
"""
@inline function ispurine(nt::NucleicAcid)
    bits = encoded_data(nt)
    return bits != 0 && (bits & 0b1010) == 0
end

"""
    ispyrimidine(nt::NucleicAcid)

Test if `nt` is surely a pyrimidine.
"""
@inline function ispyrimidine(nt::NucleicAcid)
    bits = encoded_data(nt)
    return bits != 0 && (bits & 0b0101) == 0
end

"""
    isambiguous(nt::NucleicAcid)

Test if `nt` is an ambiguous nucleotide.
"""
@inline function isambiguous(nt::NucleicAcid)
    return count_ones(nt) > 1
end

"""
    iscertain(nt::NucleicAcid)

Test if `nt` is a non-ambiguous nucleotide e.g. ACGT.
"""
@inline function iscertain(nt::NucleicAcid)
    return count_ones(nt) == 1
end

"""
    complement(nt::NucleicAcid)

Return the complementary nucleotide of `nt`.

This function returns the union of all possible complementary nucleotides.

Examples
--------

```jldoctest
julia> complement(DNA_A)
DNA_T

julia> complement(DNA_N)
DNA_N

julia> complement(RNA_U)
RNA_A

```
"""
function complement(nt::Union{DNA, RNA})
    # This is essentially a lookup table of 16 x 4 bits.
    # It's the concatenation of the bitpatterns of the nucleotides,
    # in order, complemented.
    u64 = 0xf7b3d591e6a2c480 >>> ((4 * encoded_data(nt)) & 63)
    reinterpret(typeof(nt), (u64 % UInt8) & 0x0f)
end

function Base.isvalid(::Type{T}, x::Integer) where T <: NucleicAcid
    return 0 ≤ x < 16
end

function Base.isvalid(nt::NucleicAcid)
    return encoded_data(nt) ≤ 0b1111
end

function Base.:~(x::N) where N <: NucleicAcid
    return encode(N, encoded_data(x) ⊻ 0b1111)
end

function Base.:|(x::N, y::N) where N <: NucleicAcid
    return encode(N, encoded_data(x) | encoded_data(y))
end

function Base.:&(x::N, y::N) where N <: NucleicAcid
    return encode(N, encoded_data(x) & encoded_data(y))
end


"""
    compatbits(nt::NucleicAcid)

Return the compatibility bits of `nt` as `UInt8`.
The resulting `UInt8` has the lower four bits set
if `nt` is compatible with `A`, `C`, `G` and `T/U`, respectively.

Hence, `RNA_Gap` is `0x00` (not compatible with any nucleotide),
and `DNA_W` is `0x09` (compatible with `A` and `T`)

Examples
--------

```jldoctest
julia> compatbits(DNA_A)
0x01

julia> compatbits(DNA_C)
0x02

julia> compatbits(DNA_N)
0x0f

julia> compatbits(DNA_W)
0x09

julia> compatbits(RNA_Gap)
0x00
```
"""
@inline function compatbits(nt::NucleicAcid)
    return encoded_data(nt)
end
