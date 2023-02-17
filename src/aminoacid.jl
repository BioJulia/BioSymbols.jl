# Amino Acid
# ==========
#
# Amino acid type.
#
# This file is a part of the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/NucleicAcids.jl/blob/master/LICENSE.md

"""
An amino acid type.
"""
primitive type AminoAcid <: BioSymbol 8 end

prefix(::AminoAcid) = "AA"
type_text(::AminoAcid) = "Amino Acid"
isterm(symbol::AminoAcid) = symbol === AA_Term
bytemask(symbol::AminoAcid) = 0b11111

encoded_data_eltype(::Type{AminoAcid}) = UInt8

# Conversion from/to Char
# -----------------------

function Base.convert(::Type{AminoAcid}, c::Char)
    aa = tryparse(AminoAcid, c)
    if aa === nothing
        throw(InexactError(:convert, AminoAcid, repr(c)))
    end
    return aa
end
AminoAcid(c::Char) = convert(AminoAcid, c)

Base.convert(::Type{Char}, aa::AminoAcid) = aa_to_char[encoded_data(aa) + 1]
Char(aa::AminoAcid) = convert(Char, aa)


# Amino acid encoding definition
# ------------------------------

"Invalid Amino Acid"
const AA_INVALID = encode(AminoAcid, 0x1c)  # Used during conversion from strings

# lookup table for characters
const char_to_aa = [AA_INVALID for _ in 0x00:0x7f]

# This set of amino acids is defined by IUPAC-IUB Joint Commission on Biochemical Nomenclature.
# Reference: http://www.insdc.org/documents/feature_table.html#7.4.3

const (aa_to_char, compatbits_aa) = let
    aatochar = Vector{Char}(undef, 0x1c)
    compatbitsaa = Vector{UInt32}(undef, 28)
    for (aa, doc, code) in [
        ('A', "Alanine",                           0x00),
        ('R', "Arginine",                          0x01),
        ('N', "Asparagine",                        0x02),
        ('D', "Aspartic Acid",                     0x03),
        ('C', "Cysteine",                          0x04),
        ('Q', "Glutamine",                         0x05),
        ('E', "Glutamic Acid",                     0x06),
        ('G', "Glycine",                           0x07),
        ('H', "Histidine",                         0x08),
        ('I', "Isoleucine",                        0x09),
        ('L', "Leucine",                           0x0a),
        ('K', "Lysine",                            0x0b),
        ('M', "Methionine",                        0x0c),
        ('F', "Phenylalanine",                     0x0d),
        ('P', "Proline",                           0x0e),
        ('S', "Serine",                            0x0f),
        ('T', "Threonine",                         0x10),
        ('W', "Tryptophan",                        0x11),
        ('Y', "Tyrosine",                          0x12),
        ('V', "Valine",                            0x13),
        ('O', "Pyrrolysine",                       0x14),  # non-standard
        ('U', "Selenocysteine",                    0x15),  # non-standard
        ('B', "Aspartic Acid or Asparagine",       0x16),  # ambiguous
        ('J', "Leucine or Isoleucine",             0x17),  # ambiguous
        ('Z', "Glutamine or Glutamic Acid",        0x18),  # ambiguous
        ('X', "Unspecified or Unknown Amino Acid", 0x19)]  # ambiguous
        var = Symbol("AA_", aa)
        @eval begin
            @doc $doc const $var = encode(AminoAcid, $code)
            char_to_aa[$(Int(aa)+1)] = char_to_aa[$(Int(lowercase(aa))+1)] = $var
            $(aatochar)[$(code)+1] = $aa
        end
        if code ≤ 0x15
            compatbitsaa[code+1] = 1 << code
        elseif code == 0x16
            compatbitsaa[code+1] = 1 << 0x02 | 1 << 0x03
        elseif code == 0x17
            compatbitsaa[code+1] = 1 << 0x09 | 1 << 0x0a
        elseif code == 0x18
            compatbitsaa[code+1] = 1 << 0x05 | 1 << 0x06
        elseif code == 0x19
            compatbitsaa[code+1] = (1 << 0x16) - 1
        end
    end

    @eval begin
        "Terminal"
        const AA_Term = encode(AminoAcid, 0x1a)
        char_to_aa[Int('*')+1] = AA_Term
        $(aatochar)[0x1a+1] = '*'
        $(compatbitsaa)[0x1a+1] = 1 << 0x1a

        "Amino Acid Gap"
        const AA_Gap = encode(AminoAcid, 0x1b)
        char_to_aa[Int('-') + 1] = AA_Gap
        $(aatochar)[0x1b+1] = '-'
        $(compatbitsaa)[0x1b+1] = 0
    end
    (Tuple(aatochar), Tuple(compatbitsaa))
end

@eval alphabet(::Type{AminoAcid}) = $(tuple([encode(AminoAcid, x) for x in 0x00:0x1b]...))

"""
    alphabet(AminoAcid)

Get all symbols of `AminoAcid` in sorted order.

Examples
--------

```jldoctest
julia> alphabet(AminoAcid)
(AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V, AA_O, AA_U, AA_B, AA_J, AA_Z, AA_X, AA_Term, AA_Gap)

julia> issorted(alphabet(AminoAcid))
true

```
"""
alphabet(::Type{AminoAcid})

# lookup table of 20 standard amino acids
const threeletter_to_aa = Dict(
    "ALA" => AA_A, "ARG" => AA_R, "ASN" => AA_N, "ASP" => AA_D, "CYS" => AA_C,
    "GLN" => AA_Q, "GLU" => AA_E, "GLY" => AA_G, "HIS" => AA_H, "ILE" => AA_I,
    "LEU" => AA_L, "LYS" => AA_K, "MET" => AA_M, "PHE" => AA_F, "PRO" => AA_P,
    "SER" => AA_S, "THR" => AA_T, "TRP" => AA_W, "TYR" => AA_Y, "VAL" => AA_V,
    "ASX" => AA_B, "XLE" => AA_J, "GLX" => AA_Z, "XAA" => AA_X,
    "PYL" => AA_O, "SEC" => AA_U,
)

function Base.tryparse(::Type{AminoAcid}, c::Char)
    @inbounds aa = c <= '\x7f' ? char_to_aa[Int(c)+1] : AA_INVALID
    if aa == AA_INVALID
        return nothing
    else
        return aa
    end
end

function Base.tryparse(::Type{AminoAcid}, s::AbstractString)
    if sizeof(s) == 1
        tryparse(AminoAcid, first(s))
    else
        stripped = strip(s)
        sizeof(stripped) == 1 && return tryparse(AminoAcid, first(stripped))
        sizeof(stripped) == 3 || return nothing
        return get(threeletter_to_aa, uppercase(stripped), nothing)
    end
end

function Base.parse(::Type{AminoAcid}, c::Union{AbstractString,Char})
    aa = tryparse(AminoAcid, c)
    if aa === nothing
        throw(ArgumentError("invalid amino acid"))
    end
    return aa
end

# Arithmetic and Order
# --------------------

Base.isvalid(::Type{AminoAcid}, x::Integer) = 0 ≤ x ≤ 0x1b
Base.isvalid(aa::AminoAcid) = aa ≤ AA_Gap

"""
    isambiguous(aa::AminoAcid)

Test if `aa` is an ambiguous amino acid.
"""
isambiguous(aa::AminoAcid) = AA_B ≤ aa ≤ AA_X

"""
    iscertain(aa::AminoAcid)

Test if `aa` is a non-ambiguous amino acid.
"""
function iscertain(aa::AminoAcid)
    return AA_A ≤ aa ≤ AA_U || aa == AA_Term
end

gap(::Type{AminoAcid}) = AA_Gap

"""
    compatbits(aa::AminoAcid)

Return the compatibility bits of `aa` as `UInt32`.

Examples
--------

```jldoctest
julia> compatbits(AA_A)
0x00000001

julia> compatbits(AA_J)
0x00000600

```
"""
compatbits(aa::AminoAcid) = @inbounds compatbits_aa[encoded_data(aa) + 1]
