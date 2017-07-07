```@meta
CurrentModule = BioSymbols
DocTestSetup = quote
    using BioSymbols
end
```

Nucleic Acids
=============

Type definitions
----------------

BioSymbols provides two types of nucleic acids:

| Type  | Meaning               |
| :---- | :-------------------- |
| `DNA` | deoxyribonucleic acid |
| `RNA` | ribonucleic acid      |

These two are an 8-bit primitive type and a subtype of `NucleicAcid`.

```jldoctest
julia> sizeof(DNA)
1

julia> DNA <: NucleicAcid
true

```

The set of nucleotide symbols in BioSymbols.jl covers the IUPAC nucleotides
as well as a GAP (`-`) symbol.

| Symbol | Constant              | Meaning                    |
| :----- | :-------------------- | :------------------------- |
| `'A'`  | `DNA_A` / `RNA_A`     | A; Adenine                 |
| `'C'`  | `DNA_C` / `RNA_C`     | C; Cytosine                |
| `'G'`  | `DNA_G` / `RNA_G`     | G; Guanine                 |
| `'T'`  | `DNA_T`               | T; Thymine (DNA only)      |
| `'U'`  | `RNA_U`               | U; Uracil (RNA only)       |
| `'M'`  | `DNA_M` / `RNA_M`     | A or C                     |
| `'R'`  | `DNA_R` / `RNA_R`     | A or G                     |
| `'W'`  | `DNA_W` / `RNA_W`     | A or T/U                   |
| `'S'`  | `DNA_S` / `RNA_S`     | C or G                     |
| `'Y'`  | `DNA_Y` / `RNA_Y`     | C or T/U                   |
| `'K'`  | `DNA_K` / `RNA_K`     | G or T/U                   |
| `'V'`  | `DNA_V` / `RNA_V`     | A or C or G; not T/U       |
| `'H'`  | `DNA_H` / `RNA_H`     | A or C or T; not G         |
| `'D'`  | `DNA_D` / `RNA_D`     | A or G or T/U; not C       |
| `'B'`  | `DNA_B` / `RNA_B`     | C or G or T/U; not A       |
| `'N'`  | `DNA_N` / `RNA_N`     | A or C or G or T/U         |
| `'-'`  | `DNA_Gap` / `RNA_Gap` | Gap (none of the above)    |

<http://www.insdc.org/documents/feature_table.html#7.4.1>

These are accessible as constants with `DNA_` or `RNA_` prefix:
```jldoctest
julia> DNA_A
DNA_A

julia> DNA_T
DNA_T

julia> RNA_U
RNA_U

julia> DNA_Gap
DNA_Gap

julia> typeof(DNA_A)
BioSymbols.DNA

julia> typeof(RNA_A)
BioSymbols.RNA

```

Symbols can be constructed by converting regular characters:
```jldoctest
julia> convert(DNA, 'C')
DNA_C

julia> convert(DNA, 'C') === DNA_C
true

julia> convert(DNA, 'c') === convert(DNA, 'C')  # convertion is not case-sensitive
true

```

`print` and `show` methods are defined to output the text representation of a symbol:
```jldoctest
julia> print(DNA_A)  # un-decorated text
A
julia> show(DNA_A)   # informative text
DNA_A
```


Bit encoding
------------

Every nucleotide is encoded using the lower 4 bits of a byte. An unambiguous
nucleotide has only one set bit and the other bits are unset. The table below
summarizes all unambiguous nucleotides and their corresponding bits. An
ambiguous nucleotide is the bitwise OR of unambiguous nucleotides that the
ambiguous nucleotide can take. For example, `DNA_R` (meaning the nucleotide is
either `DNA_A` or `DNA_G`) is encoded as `0101` because `0101` is the bitwise OR
of `0001` (`DNA_A`) and `0100` (`DNA_G`). The gap symbol is always `0000`.

```jldoctest
julia> bits(reinterpret(UInt8, DNA_A))
"00000001"

julia> bits(reinterpret(UInt8, DNA_G))
"00000100"

julia> bits(reinterpret(UInt8, DNA_R))
"00000101"

```

This bit encoding enables efficient bit operations:

```jldoctest
julia> DNA_A | DNA_G  # A or G
DNA_R

julia> DNA_A & DNA_G  # A and G
DNA_Gap

julia> DNA_A | ~DNA_A  # A or not A
DNA_N

julia> DNA_A | DNA_C | DNA_G | DNA_T  # any DNA
DNA_N

```
