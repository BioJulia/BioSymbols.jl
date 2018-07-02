```@meta
CurrentModule = BioSymbols
DocTestSetup = quote
    using BioSymbols
end
```

Amino Acids
===========

Set of amino acid symbols also covers IUPAC amino acid symbols plus a gap
symbol:

| Symbol | Constant  | Meaning                     |
| :----- | :-------- | :-------------------------- |
| `'A'`  | `AA_A`    | Alanine                     |
| `'R'`  | `AA_R`    | Arginine                    |
| `'N'`  | `AA_N`    | Asparagine                  |
| `'D'`  | `AA_D`    | Aspartic acid (Aspartate)   |
| `'C'`  | `AA_C`    | Cysteine                    |
| `'Q'`  | `AA_Q`    | Glutamine                   |
| `'E'`  | `AA_E`    | Glutamic acid (Glutamate)   |
| `'G'`  | `AA_G`    | Glycine                     |
| `'H'`  | `AA_H`    | Histidine                   |
| `'I'`  | `AA_I`    | Isoleucine                  |
| `'L'`  | `AA_L`    | Leucine                     |
| `'K'`  | `AA_K`    | Lysine                      |
| `'M'`  | `AA_M`    | Methionine                  |
| `'F'`  | `AA_F`    | Phenylalanine               |
| `'P'`  | `AA_P`    | Proline                     |
| `'S'`  | `AA_S`    | Serine                      |
| `'T'`  | `AA_T`    | Threonine                   |
| `'W'`  | `AA_W`    | Tryptophan                  |
| `'Y'`  | `AA_Y`    | Tyrosine                    |
| `'V'`  | `AA_V`    | Valine                      |
| `'O'`  | `AA_O`    | Pyrrolysine                 |
| `'U'`  | `AA_U`    | Selenocysteine              |
| `'B'`  | `AA_B`    | Aspartic acid or Asparagine |
| `'J'`  | `AA_J`    | Leucine or Isoleucine       |
| `'Z'`  | `AA_Z`    | Glutamine or Glutamic acid  |
| `'X'`  | `AA_X`    | Any amino acid              |
| `'*'`  | `AA_Term` | Termination codon           |
| `'-'`  | `AA_Gap`  | Gap (none of the above)     |

<http://www.insdc.org/documents/feature_table.html#7.4.3>

Symbols are accessible as constants with `AA_` prefix:
```jldoctest
julia> AA_A
AA_A

julia> AA_Q
AA_Q

julia> AA_Term
AA_Term

julia> typeof(AA_A)
BioSymbols.AminoAcid

```

Symbols can be constructed by converting regular characters:
```jldoctest
julia> convert(AminoAcid, 'A')
AA_A

julia> convert(AminoAcid, 'P') === AA_P
true

julia> convert(AminoAcid, 'a') === convert(AminoAcid, 'A')
true

```

3-letter and 1-letter abbreviations can be parsed using `parse` in a
case-insensitive way:
```jldoctest
julia> parse(AminoAcid, "Pro")  # 3-letter abbreviation
AA_P

julia> parse(AminoAcid, "P")    # 1-letter abbreviation
AA_P

julia> parse(AminoAcid, "Pro") == parse(AminoAcid, "pRo")
true

julia> tryparse(AminoAcid, "Pro")  # tryparse returns either an amino acid or nothing
AA_P

julia> tryparse(AminoAcid, "Pr")
nothing

```
