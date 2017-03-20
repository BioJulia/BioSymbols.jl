# typedefs.jl
# ===========
#
# Definitions of NucleicAcid and AminoAcid types.
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

@compat abstract type NucleicAcid end
@compat primitive type DNA <: NucleicAcid 8 end
@compat primitive type RNA <: NucleicAcid 8 end

# Conversion from/to integers
# ---------------------------

Base.convert(::Type{DNA}, nt::UInt8) = reinterpret(DNA, nt)
Base.convert(::Type{RNA}, nt::UInt8) = reinterpret(RNA, nt)
Base.convert(::Type{UInt8}, nt::DNA) = reinterpret(UInt8, nt)
Base.convert(::Type{UInt8}, nt::RNA) = reinterpret(UInt8, nt)
Base.convert{T<:Number,S<:NucleicAcid}(::Type{T}, nt::S) = convert(T, UInt8(nt))
Base.convert{T<:Number,S<:NucleicAcid}(::Type{S}, nt::T) = convert(S, UInt8(nt))

# Conversion from/to characters
# -----------------------------

function Base.convert(::Type{DNA}, c::Char)
    if c > '\xff'
        throw(InexactError())
    end
    @inbounds dna = char_to_dna[Int(c) + 1]
    if !isvalid(DNA, dna)
        throw(InexactError())
    end
    return reinterpret(DNA, dna)
end

function Base.convert(::Type{RNA}, c::Char)
    if c > '\xff'
        throw(InexactError())
    end
    @inbounds rna = char_to_rna[Int(c) + 1]
    if !isvalid(RNA, rna)
        throw(InexactError())
    end
    return reinterpret(RNA, rna)
end

function Base.convert(::Type{Char}, nt::DNA)
    return dna_to_char[convert(UInt8, nt) + 1]
end

function Base.convert(::Type{Char}, nt::RNA)
    return rna_to_char[convert(UInt8, nt) + 1]
end

# Encoding of DNA and RNA NucleicAcids
# ------------------------------------
# DNA

# lookup table for characters
const char_to_dna = [0x80 for _ in 0x00:0xff]
const dna_to_char = Vector{Char}(16)

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
        char_to_dna[$(Int(char)+1)] = char_to_dna[$(Int(lowercase(char))+1)] = $(bits)
        dna_to_char[$(Int(bits)+1)] = $(char)
    end
end

@eval alphabet(::Type{DNA}) = $(tuple([reinterpret(DNA, x)
                                                 for x in 0b0000:0b1111]...))

const ACGT = (DNA_A, DNA_C, DNA_G, DNA_T)
const ACGTN = (DNA_A, DNA_C, DNA_G, DNA_T, DNA_N)

# RNA

# lookup table for characters
const char_to_rna = [0x80 for _ in 0x00:0xff]
const rna_to_char = Vector{Char}(16)

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
        char_to_rna[$(Int(char)+1)] = char_to_rna[$(Int(lowercase(char)+1))] = reinterpret(UInt8, $(dna))
        rna_to_char[$(Int(dna)+1)] = $(char)
    end
end

@eval alphabet(::Type{RNA}) = $(tuple([reinterpret(RNA, x)
                                                 for x in 0b0000:0b1111]...))

const ACGU = (RNA_A, RNA_C, RNA_G, RNA_U)
const ACGUN = (RNA_A, RNA_C, RNA_G, RNA_U, RNA_N)


# Aminoacids
# ==========
#
# The amino acid type.

"Type representing AminoAcids"
@compat primitive type AminoAcid 8 end

# Conversion from/to integers
# ---------------------------

Base.convert(::Type{AminoAcid}, aa::UInt8) = reinterpret(AminoAcid, aa)
Base.convert(::Type{UInt8}, aa::AminoAcid) = reinterpret(UInt8, aa)
Base.convert{T<:Number}(::Type{T}, aa::AminoAcid) = convert(T, UInt8(aa))
Base.convert{T<:Number}(::Type{AminoAcid}, aa::T) = convert(AminoAcid, UInt8(aa))

# Conversion from/to Char
# -----------------------

function Base.convert(::Type{AminoAcid}, c::Char)
    @inbounds aa = c <= '\x7f' ? char_to_aa[Int(c)+1] : AA_INVALID
    @assert aa != AA_INVALID error("$(c) is not a valid amino acid")
    return aa
end

Base.convert(::Type{Char}, aa::AminoAcid) = aa_to_char[convert(UInt8, aa) + 1]

# Amino acid encoding definition
# ------------------------------

"Invalid Amino Acid"
const AA_INVALID = convert(AminoAcid, 0x1c)  # Used during conversion from strings

# lookup table for characters
const char_to_aa = [AA_INVALID for _ in 0x00:0x7f]
const aa_to_char = Vector{Char}(0x1c)

# compatibility bits
const compatbits_aa = Vector{UInt32}(28)

# This set of amino acids is defined by IUPAC-IUB Joint Commission on Biochemical Nomenclature.
# Reference: http://www.insdc.org/documents/feature_table.html#7.4.3

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
        @doc $doc const $var = convert(AminoAcid, $code)
        char_to_aa[$(Int(aa)+1)] = char_to_aa[$(Int(lowercase(aa))+1)] = $var
        aa_to_char[$(code)+1] = $aa
    end
    if code ≤ 0x15
        compatbits_aa[code+1] = 1 << code
    elseif code == 0x16
        compatbits_aa[code+1] = 1 << 0x02 | 1 << 0x03
    elseif code == 0x17
        compatbits_aa[code+1] = 1 << 0x09 | 1 << 0x0a
    elseif code == 0x18
        compatbits_aa[code+1] = 1 << 0x05 | 1 << 0x06
    elseif code == 0x19
        compatbits_aa[code+1] = (1 << 0x16) - 1
    end
end

"Terminal"
const AA_Term = convert(AminoAcid, 0x1a)
char_to_aa[Int('*')+1] = AA_Term
aa_to_char[0x1a+1] = '*'
compatbits_aa[0x1a+1] = 1 << 0x1a

"Amino Acid Gap"
const AA_Gap = convert(AminoAcid, 0x1b)
char_to_aa[Int('-') + 1] = AA_Gap
aa_to_char[0x1b+1] = '-'
compatbits_aa[0x1b+1] = 0

@eval alphabet(::Type{AminoAcid}) = $(tuple([reinterpret(AminoAcid, x) for x in 0x00:0x1b]...))

# lookup table of 20 standard amino acids
const threeletter_to_aa = Dict(
    "ALA" => AA_A, "ARG" => AA_R, "ASN" => AA_N, "ASP" => AA_D, "CYS" => AA_C,
    "GLN" => AA_Q, "GLU" => AA_E, "GLY" => AA_G, "HIS" => AA_H, "ILE" => AA_I,
    "LEU" => AA_L, "LYS" => AA_K, "MET" => AA_M, "PHE" => AA_F, "PRO" => AA_P,
    "SER" => AA_S, "THR" => AA_T, "TRP" => AA_W, "TYR" => AA_Y, "VAL" => AA_V,
    "ASX" => AA_B, "XLE" => AA_J, "GLX" => AA_Z, "XAA" => AA_X,
    "PYL" => AA_O, "SEC" => AA_U,
)

function Base.parse(::Type{AminoAcid}, s::AbstractString)
    s′ = strip(s)
    if length(s′) == 1
        return convert(AminoAcid, s′[1])
    end
    try
        return threeletter_to_aa[uppercase(s′)]
    catch ex
        if isa(ex, KeyError)
            error("invalid amino acid string: \"$s\" ")
        end
        rethrow()
    end
end
