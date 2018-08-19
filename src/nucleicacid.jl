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
abstract type NucleicAcid end

"""
A deoxyribonucleic acid type.
"""
primitive type DNA <: NucleicAcid 8 end

"""
A ribonucleic acid type.
"""
primitive type RNA <: NucleicAcid 8 end


Base.length(::NucleicAcid) = 1
Base.iterate(nt::NucleicAcid) = (nt, nothing)
Base.iterate(nt::NucleicAcid, state) = nothing


# Conversion from/to integers
# ---------------------------

Base.convert(::Type{T}, nt::UInt8) where T <: NucleicAcid = reinterpret(T, nt)
Base.convert(::Type{UInt8}, nt::T) where T <: NucleicAcid = reinterpret(UInt8, nt)
Base.convert(::Type{T}, nt::S) where {T <: Number, S <: NucleicAcid} = convert(T, convert(UInt8, nt))
Base.convert(::Type{S}, nt::T) where {T <: Number, S <: NucleicAcid} = convert(S, convert(UInt8, nt))
DNA(nt::Integer) = convert(DNA, nt)
RNA(nt::Integer) = convert(RNA, nt)

Base.convert(::Type{DNA}, nt::RNA) = DNA(convert(UInt8, nt))
DNA(nt::RNA) = convert(DNA, nt)
Base.convert(::Type{RNA}, nt::DNA) = RNA(convert(UInt8, nt))
RNA(nt::DNA) = convert(RNA, nt)


# Conversion from/to characters
# -----------------------------

function Base.convert(::Type{DNA}, c::Char)
    if c > '\uff'
        throw(InexactError(:convert, DNA, c))
    end
    @inbounds dna = char_to_dna[convert(Int, c) + 1]
    if !isvalid(DNA, dna)
        throw(InexactError(:convert, DNA, c))
    end
    return reinterpret(DNA, dna)
end
DNA(c::Char) = convert(DNA, c)

function Base.convert(::Type{RNA}, c::Char)
    if c > '\uff'
        throw(InexactError(:convert, RNA, c))
    end
    @inbounds rna = char_to_rna[convert(Int, c) + 1]
    if !isvalid(RNA, rna)
        throw(InexactError(:convert, RNA, c))
    end
    return reinterpret(RNA, rna)
end
RNA(c::Char) = convert(RNA, c)

function Base.convert(::Type{Char}, nt::DNA)
    return dna_to_char[convert(UInt8, nt) + 1]
end
Char(nt::DNA) = convert(Char, nt)

function Base.convert(::Type{Char}, nt::RNA)
    return rna_to_char[convert(UInt8, nt) + 1]
end
Char(nt::RNA) = convert(Char, nt)

# Print
# -----

prefix(::Type{DNA}) = "DNA"
prefix(::Type{RNA}) = "RNA"

function Base.show(io::IO, nt::T) where T <: NucleicAcid
    if isvalid(nt)
        if nt == gap(T)
            write(io, prefix(T), "_Gap")
        else
            write(io, prefix(T), "_", convert(Char, nt))
        end
    else
        write(io, "Invalid ", prefix(T))
    end
    return
end

function Base.print(io::IO, nt::NucleicAcid)
    if !isvalid(nt)
        throw(ArgumentError("nucleic acid is invalid"))
    end
    write(io, convert(Char, nt))
    return
end

Base.write(io::IO, na::NucleicAcid) = write(io, reinterpret(UInt8, na))
Base.read(io::IO, ::Type{T}) where T<:NucleicAcid = reinterpret(T, read(io, UInt8))

# Encoding of DNA and RNA NucleicAcids
# ------------------------------------

# lookup table for characters
const char_to_dna = [0x80 for _ in 0x00:0xff]
const dna_to_char = Vector{Char}(undef, 16)

# derived from "The DDBJ/ENA/GenBank Feature Table Definition"
# §7.4.1 Nucleotide base code (IUPAC)
# http://www.insdc.org/documents/feature_table.html#7.4.1
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
        @doc $(doc) const $(var) = reinterpret(DNA, $(bits))
        char_to_dna[$(convert(Int, char) + 1)] = char_to_dna[$(convert(Int, lowercase(char)) + 1)] = $(bits)
        dna_to_char[$(convert(Int, bits) + 1)] = $(char)
    end
end

@eval function alphabet(::Type{DNA})
    return $(tuple([reinterpret(DNA, x) for x in 0b0000:0b1111]...))
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

# lookup table for characters
const char_to_rna = [0x80 for _ in 0x00:0xff]
const rna_to_char = Vector{Char}(undef, 16)

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
        rna_to_char[$(convert(Int, dna) + 1)] = $(char)
    end
end

@eval function alphabet(::Type{RNA})
    return $(tuple([reinterpret(RNA, x) for x in 0b0000:0b1111]...))
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

function Base.:~(x::N) where N <: NucleicAcid
    return reinterpret(N, ~reinterpret(UInt8, x) & 0b1111)
end

function Base.:|(x::N, y::N) where N <: NucleicAcid
    return reinterpret(N, reinterpret(UInt8, x) | reinterpret(UInt8, y))
end

function Base.:&(x::N, y::N) where N <: NucleicAcid
    return reinterpret(N, reinterpret(UInt8, x) & reinterpret(UInt8, y))
end

function Base.:-(x::N, y::N) where N <: NucleicAcid
    return convert(Int, x) - convert(Int, y)
end

function Base.:-(x::N, y::Integer) where N <: NucleicAcid
    return x + (-y)
end

function Base.:+(x::N, y::Integer) where N <: NucleicAcid
    return reinterpret(N, (convert(UInt8, x) + y % UInt8) & 0b1111)
end

function Base.isless(x::N, y::N) where N <: NucleicAcid
    return isless(reinterpret(UInt8, x), reinterpret(UInt8, y))
end

@inline function Base.count_ones(nt::NucleicAcid)
    return count_ones(reinterpret(UInt8, nt))
end

function Base.trailing_zeros(nt::NucleicAcid)
    return trailing_zeros(reinterpret(UInt8, nt))
end

"""
    gap(DNA)

Return `DNA_Gap`.
"""
function gap(::Type{DNA})
    return DNA_Gap
end

"""
    gap(RNA)

Return `RNA_Gap`.
"""
function gap(::Type{RNA})
    return RNA_Gap
end

"""
    isGC(nt::NucleicAcid)

Test if `nt` is surely either guanine or cytosine.
"""
function isGC(nt::NucleicAcid)
    bits = reinterpret(UInt8, nt)
    return bits != 0 && (bits & 0b1001) == 0
end

"""
    ispurine(nt::NucleicAcid)

Test if `nt` is surely a purine.
"""
@inline function ispurine(nt::NucleicAcid)
    bits = reinterpret(UInt8, nt)
    return bits != 0 && (bits & 0b1010) == 0
end

"""
    ispyrimidine(nt::NucleicAcid)

Test if `nt` is surely a pyrimidine.
"""
@inline function ispyrimidine(nt::NucleicAcid)
    bits = reinterpret(UInt8, nt)
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
    isgap(nt::NucleicAcid)

Test if `nt` is a gap.
"""
@inline function isgap(nt::NucleicAcid)
    return count_ones(nt) == 0
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
function complement(nt::NucleicAcid)
    bits = compatbits(nt)
    return reinterpret(
        typeof(nt),
        (bits & 0x01) << 3 | (bits & 0x08) >> 3 |
        (bits & 0x02) << 1 | (bits & 0x04) >> 1)
end

function Base.isvalid(::Type{T}, x::Integer) where T <: NucleicAcid
    return 0 ≤ x < 16
end

function Base.isvalid(nt::NucleicAcid)
    return reinterpret(UInt8, nt) ≤ 0b1111
end

"""
    iscompatible(x::T, y::T) where T <: NucleicAcid

Test if `x` and `y` are compatible with each other (i.e. `x` and `y` can be the same symbol).

`x` and `y` must be the same type.

Examples
--------

```jldoctest
julia> iscompatible(DNA_A, DNA_A)
true

julia> iscompatible(DNA_C, DNA_N)  # DNA_N can be DNA_C
true

julia> iscompatible(DNA_C, DNA_R)  # DNA_R (A or G) cannot be DNA_C
false

```
"""
@inline function iscompatible(x::T, y::T) where T <: NucleicAcid
    return compatbits(x) & compatbits(y) != 0
end

"""
    compatbits(nt::NucleicAcid)

Return the compatibility bits of `nt` as `UInt8`.

Examples
--------

```jldoctest
julia> compatbits(DNA_A)
0x01

julia> compatbits(DNA_C)
0x02

julia> compatbits(DNA_N)
0x0f

```
"""
@inline function compatbits(nt::NucleicAcid)
    return reinterpret(UInt8, nt)
end

#=
# TODO - extra boundschecking which can be turned off using @inbounds.

@inline function Base.getindex(arr::SVector{4}, idx::NucleicAcid)
    @inbounds return arr[trailing_zeros(idx) + 1]
end

@inline function Base.getindex(arr::SVector{16}, idx::NucleicAcid)
    @inbounds return arr[reinterpret(UInt8, idx) + 1]
end

@inline function Base.getindex(arr::SMatrix{4,4}, i::NucleicAcid, j::NucleicAcid)
    @inbounds return arr[trailing_zeros(i) + 1, trailing_zeros(j) + 1]
end

@inline function Base.getindex(arr::SMatrix{16,16}, i::NucleicAcid, j::NucleicAcid)
    @inbounds return arr[reinterpret(UInt8, i) + 1, reinterpret(UInt8, j) + 1]
end
=#
