```@meta
CurrentModule = BioSymbols
DocTestSetup = quote
    using BioSymbols
end
```

Sequences
=========

Using `Vector`
--------------

A quick way to create a DNA/RNA sequence is storing symbols in a vector.

```jldoctest
julia> seq = [DNA_A, DNA_C, DNA_G, DNA_T]
4-element Array{DNA,1}:
 DNA_A
 DNA_C
 DNA_G
 DNA_T

julia> [convert(DNA, x) for x in "ACGT"]  # from a string
4-element Array{DNA,1}:
 DNA_A
 DNA_C
 DNA_G
 DNA_T

```


Using `Tuple`
-------------

Julia offers a tuple type to represent multiple values in a value. It is similar
to `Vector` but is significantly different in some points. First, `Tuple` is
immutable while `Vector` is mutable. So you cannot update elements in a tuple
once created. Second, a tuple type is parameterized by its length. That means it
is inefficient to represent variable-length sequences in tuple due to type
instability problem.

```jldoctest
julia> (RNA_A, RNA_U, RNA_C)  # RNA triplet (or codon)
(RNA_A, RNA_U, RNA_C)

```


Using the BioSequences package
------------------------------

Using `Vector` or `Tuple` is simple, however, BioSymbols does not offer useful
operations for these representations. So you need to use built-in operations of
Julia or other packages. Moreover, these representations are not necessarily
efficient. For example, `DNA` is an 8-bit primitive but it only uses 4 bits,
which means 50% of a `Vector{DNA}`'s space is not used at all.

For the purpose of representing sequences as efficient as possible BioJulia has
developed [BioSequences](https://github.com/BioJulia/BioSequences.jl)
package. The `BioSequence` type is able to represent a DNA/RNA sequence in 2 or
4 bits per symbol. It also offers many efficient algorithms and I/O tools for
common file formats such as FASTA.
