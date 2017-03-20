# bitops.jl
# =========
#
# Bit Operations for nucleic acid types.
#
# This file is a part of the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/NucleicAcids.jl/blob/master/LICENSE.md

# NucleicAcids
# ============

function Base.:~{N<:NucleicAcid}(x::N)
    return reinterpret(N, ~reinterpret(UInt8, x) & 0b1111)
end

function Base.:|{N<:NucleicAcid}(x::N, y::N)
    return reinterpret(N, reinterpret(UInt8, x) | reinterpret(UInt8, y))
end

function Base.:&{N<:NucleicAcid}(x::N, y::N)
    return reinterpret(N, reinterpret(UInt8, x) & reinterpret(UInt8, y))
end

function Base.:-{N<:NucleicAcid}(x::N, y::N)
    return Int(x) - Int(y)
end

function Base.:-{N<:NucleicAcid}(x::N, y::Integer)
    return x + (-y)
end

function Base.:+{N<:NucleicAcid}(x::N, y::Integer)
    return reinterpret(N, (UInt8(x) + y % UInt8) & 0b1111)
end

function Base.isless{N<:NucleicAcid}(x::N, y::N)
    return isless(reinterpret(UInt8, x), reinterpret(UInt8, y))
end

@inline function Base.count_ones(nt::NucleicAcid)
    return count_ones(reinterpret(UInt8, nt))
end

function Base.trailing_zeros(nt::NucleicAcid)
    return trailing_zeros(reinterpret(UInt8, nt))
end

function gap{N<:NucleicAcid}(::Type{N})
    return reinterpret(N, 0b0000)
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
Test if nucleotide is surely a purine.
"""
@inline function ispurine(nt::NucleicAcid)
    bits = reinterpret(UInt8, nt)
    return bits != 0 && (bits & 0b1010) == 0
end

"""
    ispyrimidine(nt::NucleicAcid)
Test if nucleotide is surely a pyrimidine.
"""
@inline function ispyrimidine(nt::NucleicAcid)
    bits = reinterpret(UInt8, nt)
    return bits != 0 && (bits & 0b0101) == 0
end

"""
    isambiguous(nt::NucleicAcid)
Test if `nt` is ambiguous nucleotide.
"""
@inline function isambiguous(nt::NucleicAcid)
    return count_ones(nt) > 1
end

"""
    iscertain(nt::NucleicAcid)
Test if `nt` is a non-ambiguous nucleotide e.g. ATCG.
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
"""
function complement(nt::NucleicAcid)
    bits = compatbits(nt)
    return reinterpret(
        typeof(nt),
        (bits & 0x01) << 3 | (bits & 0x08) >> 3 |
        (bits & 0x02) << 1 | (bits & 0x04) >> 1)
end

function Base.isvalid{T<:NucleicAcid}(::Type{T}, x::Integer)
    return 0 ≤ x < 16
end

function Base.isvalid(nt::NucleicAcid)
    return reinterpret(UInt8, nt) ≤ 0b1111
end

"""
    iscompatible(x, y)
Return `true` if and only if `x` and `y` are compatible with each other (i.e.
`x` and `y` can be the same symbol).
`x` and `y` must be the same type (`DNA`, `RNA` or `AminoAcid`).
# Examples
```julia
julia> iscompatible(DNA_A, DNA_A)
true
julia> iscompatible(DNA_C, DNA_N)  # DNA_N can be DNA_C
true
julia> iscompatible(DNA_C, DNA_R)  # DNA_R (A or G) cannot be DNA_C
false
julia> iscompatible(AA_A, AA_X)    # AA_X can be AA_A
true
```
"""
@inline function iscompatible{T<:NucleicAcid}(x::T, y::T)
    return compatbits(x) & compatbits(y) != 0
end

# Return the compatibility bits of `nt`.
@inline function compatbits(nt::NucleicAcid)
    return reinterpret(UInt8, nt)
end

# AminoAcids
# ==========

# Arithmetic and Order
# --------------------

# These methods are necessary when deriving some algorithims
# like iteration, sort, comparison, and so on.
Base.:-(x::AminoAcid, y::AminoAcid) = Int(x) - Int(y)
# 0x1c is the size of the amino acid alphabet
Base.:-(x::AminoAcid, y::Integer) = x + mod(-y, 0x1c)
Base.:+(x::AminoAcid, y::Integer) = reinterpret(AminoAcid, mod((UInt8(x) + y) % UInt8, 0x1c))
Base.:+(x::Integer, y::AminoAcid) = y + x
Base.isless(x::AminoAcid, y::AminoAcid) = isless(UInt8(x), UInt8(y))

Base.isvalid(::Type{AminoAcid}, x::Integer) = 0 ≤ x ≤ 0x1b
Base.isvalid(aa::AminoAcid) = aa ≤ AA_Gap
isambiguous(aa::AminoAcid) = AA_B ≤ aa ≤ AA_X
gap(::Type{AminoAcid}) = AA_Gap

compatbits(aa::AminoAcid) = compatbits_aa[reinterpret(UInt8, aa)+1]

function iscompatible(x::AminoAcid, y::AminoAcid)
    return compatbits(x) & compatbits(y) != 0
end
