# Bit encoding of nucleic acid types

Unambiguous nucleotides are represented in one-hot encoding as follows:

| NucleicAcid | Bits |
| ----------- | ---- |
|     A       | 0001 |
|     C       | 0010 |
|     G       | 0100 |
|    T/U      | 1000 |

Ambiguous nucleotides are the bitwise OR of these four nucleotides.
For example, R, A or G, is represented as 0101 (= A: 0001 | G: 0100).
The gap symbol is always 0000.
The meaningful four bits are stored in the least significant bits of a byte.

This encoding applies to both the `DNA` and `RNA` types.