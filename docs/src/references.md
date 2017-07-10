References
==========

```@meta
CurrentModule = BioSymbols
DocTestSetup = quote
    using BioSymbols
end
```

Nucleic acids
-------------

```@docs
NucleicAcid
```

### DNA

```@docs
DNA
DNA_A
DNA_C
DNA_G
DNA_T
DNA_M
DNA_R
DNA_W
DNA_S
DNA_Y
DNA_K
DNA_V
DNA_H
DNA_D
DNA_B
DNA_N
DNA_Gap
ACGT
ACGTN
```

### RNA

```@docs
RNA
RNA_A
RNA_C
RNA_G
RNA_U
RNA_M
RNA_R
RNA_W
RNA_S
RNA_Y
RNA_K
RNA_V
RNA_H
RNA_D
RNA_B
RNA_N
RNA_Gap
ACGU
ACGUN
```

### Functions

```@docs
alphabet(::Type{DNA})
alphabet(::Type{RNA})
gap(::Type{DNA})
gap(::Type{RNA})
complement(::NucleicAcid)
isgap(::NucleicAcid)
compatbits(::NucleicAcid)
iscompatible{T<:NucleicAcid}(::T, ::T)
isambiguous(::NucleicAcid)
iscertain(::NucleicAcid)
isGC(::NucleicAcid)
ispurine(::NucleicAcid)
ispyrimidine(::NucleicAcid)
```


Amino acids
-----------

### Amino acids

```@docs
AminoAcid
AA_A
AA_R
AA_N
AA_D
AA_C
AA_Q
AA_E
AA_G
AA_H
AA_I
AA_L
AA_K
AA_M
AA_F
AA_P
AA_S
AA_T
AA_W
AA_Y
AA_V
AA_O
AA_U
AA_B
AA_J
AA_Z
AA_X
AA_Term
AA_Gap
```

### Functions

```@docs
alphabet(::Type{AminoAcid})
gap(::Type{AminoAcid})
isgap(::AminoAcid)
compatbits(::AminoAcid)
iscompatible(::AminoAcid, ::AminoAcid)
isambiguous(::AminoAcid)
iscertain(::AminoAcid)
```
