
using BioSymbols

using Test

function round_trip(x)
    io = IOBuffer()
    write(io, x)
    io = IOBuffer(take!(io))
    y = read(io, typeof(x))
    return x == y
end

@testset "NucleicAcids" begin
    @testset "Conversions" begin
        @testset "UInt8" begin
            @testset "DNA conversions from UInt8" begin
                @test convert(DNA, 0b0000) === DNA_Gap
                @test convert(DNA, 0b0001) === DNA_A
                @test convert(DNA, 0b0010) === DNA_C
                @test convert(DNA, 0b0011) === DNA_M
                @test convert(DNA, 0b0100) === DNA_G
                @test convert(DNA, 0b0101) === DNA_R
                @test convert(DNA, 0b0110) === DNA_S
                @test convert(DNA, 0b0111) === DNA_V
                @test convert(DNA, 0b1000) === DNA_T
                @test convert(DNA, 0b1001) === DNA_W
                @test convert(DNA, 0b1010) === DNA_Y
                @test convert(DNA, 0b1011) === DNA_H
                @test convert(DNA, 0b1100) === DNA_K
                @test convert(DNA, 0b1101) === DNA_D
                @test convert(DNA, 0b1110) === DNA_B
                @test convert(DNA, 0b1111) === DNA_N
            end

            @testset "RNA conversions from UInt8" begin
                @test convert(RNA, 0b0000) === RNA_Gap
                @test convert(RNA, 0b0001) === RNA_A
                @test convert(RNA, 0b0010) === RNA_C
                @test convert(RNA, 0b0011) === RNA_M
                @test convert(RNA, 0b0100) === RNA_G
                @test convert(RNA, 0b0101) === RNA_R
                @test convert(RNA, 0b0110) === RNA_S
                @test convert(RNA, 0b0111) === RNA_V
                @test convert(RNA, 0b1000) === RNA_U
                @test convert(RNA, 0b1001) === RNA_W
                @test convert(RNA, 0b1010) === RNA_Y
                @test convert(RNA, 0b1011) === RNA_H
                @test convert(RNA, 0b1100) === RNA_K
                @test convert(RNA, 0b1101) === RNA_D
                @test convert(RNA, 0b1110) === RNA_B
                @test convert(RNA, 0b1111) === RNA_N
            end

            @testset "DNA conversions to UInt8" begin
                @test convert(UInt8, DNA_Gap) === 0b0000
                @test convert(UInt8, DNA_A)   === 0b0001
                @test convert(UInt8, DNA_C)   === 0b0010
                @test convert(UInt8, DNA_M)   === 0b0011
                @test convert(UInt8, DNA_G)   === 0b0100
                @test convert(UInt8, DNA_R)   === 0b0101
                @test convert(UInt8, DNA_S)   === 0b0110
                @test convert(UInt8, DNA_V)   === 0b0111
                @test convert(UInt8, DNA_T)   === 0b1000
                @test convert(UInt8, DNA_W)   === 0b1001
                @test convert(UInt8, DNA_Y)   === 0b1010
                @test convert(UInt8, DNA_H)   === 0b1011
                @test convert(UInt8, DNA_K)   === 0b1100
                @test convert(UInt8, DNA_D)   === 0b1101
                @test convert(UInt8, DNA_B)   === 0b1110
                @test convert(UInt8, DNA_N)   === 0b1111
            end

            @testset "RNA conversions to UInt8" begin
                @test convert(UInt8, RNA_Gap) === 0b0000
                @test convert(UInt8, RNA_A)   === 0b0001
                @test convert(UInt8, RNA_C)   === 0b0010
                @test convert(UInt8, RNA_M)   === 0b0011
                @test convert(UInt8, RNA_G)   === 0b0100
                @test convert(UInt8, RNA_R)   === 0b0101
                @test convert(UInt8, RNA_S)   === 0b0110
                @test convert(UInt8, RNA_V)   === 0b0111
                @test convert(UInt8, RNA_U)   === 0b1000
                @test convert(UInt8, RNA_W)   === 0b1001
                @test convert(UInt8, RNA_Y)   === 0b1010
                @test convert(UInt8, RNA_H)   === 0b1011
                @test convert(UInt8, RNA_K)   === 0b1100
                @test convert(UInt8, RNA_D)   === 0b1101
                @test convert(UInt8, RNA_B)   === 0b1110
                @test convert(UInt8, RNA_N)   === 0b1111
            end
        end

        @testset "UInt64" begin
            @testset "DNA conversions from UInt64" begin
                @test convert(DNA, UInt64(0b0000)) === DNA_Gap
                @test convert(DNA, UInt64(0b0001)) === DNA_A
                @test convert(DNA, UInt64(0b0010)) === DNA_C
                @test convert(DNA, UInt64(0b0100)) === DNA_G
                @test convert(DNA, UInt64(0b1000)) === DNA_T
                @test convert(DNA, UInt64(0b1111)) === DNA_N
            end

            @testset "RNA conversions from UInt64" begin
                @test convert(RNA, UInt64(0b0000)) === RNA_Gap
                @test convert(RNA, UInt64(0b0001)) === RNA_A
                @test convert(RNA, UInt64(0b0010)) === RNA_C
                @test convert(RNA, UInt64(0b0100)) === RNA_G
                @test convert(RNA, UInt64(0b1000)) === RNA_U
                @test convert(RNA, UInt64(0b1111)) === RNA_N
            end

            @testset "DNA conversions to UInt64" begin
                @test convert(UInt64, DNA_Gap) === UInt64(0b0000)
                @test convert(UInt64, DNA_A)   === UInt64(0b0001)
                @test convert(UInt64, DNA_C)   === UInt64(0b0010)
                @test convert(UInt64, DNA_G)   === UInt64(0b0100)
                @test convert(UInt64, DNA_T)   === UInt64(0b1000)
                @test convert(UInt64, DNA_N)   === UInt64(0b1111)
            end

            @testset "RNA conversions to UInt64" begin
                @test convert(UInt64, RNA_Gap) === UInt64(0b0000)
                @test convert(UInt64, RNA_A)   === UInt64(0b0001)
                @test convert(UInt64, RNA_C)   === UInt64(0b0010)
                @test convert(UInt64, RNA_G)   === UInt64(0b0100)
                @test convert(UInt64, RNA_U)   === UInt64(0b1000)
                @test convert(UInt64, RNA_N)   === UInt64(0b1111)
            end
        end

        @testset "Char" begin
            @testset "DNA conversions from Char" begin
                @test convert(DNA, 'A') === DNA('A') === DNA_A
                @test convert(DNA, 'C') === DNA('C') === DNA_C
                @test convert(DNA, 'G') === DNA('G') === DNA_G
                @test convert(DNA, 'T') === DNA('T') === DNA_T
                @test convert(DNA, 'N') === DNA('N') === DNA_N
                @test_throws InexactError convert(DNA, 'Z')
                @test_throws InexactError convert(DNA, '核')
            end

            @testset "RNA conversions from Char" begin
                @test convert(RNA, 'A') === RNA('A') === RNA_A
                @test convert(RNA, 'C') === RNA('C') === RNA_C
                @test convert(RNA, 'G') === RNA('G') === RNA_G
                @test convert(RNA, 'U') === RNA('U') === RNA_U
                @test convert(RNA, 'N') === RNA('N') === RNA_N
                @test_throws InexactError convert(RNA, 'Z')
                @test_throws InexactError convert(RNA, '核')
            end

            @testset "DNA conversions to Char" begin
                @test convert(Char, DNA_A) === Char(DNA_A) === 'A'
                @test convert(Char, DNA_C) === Char(DNA_C) === 'C'
                @test convert(Char, DNA_G) === Char(DNA_G) === 'G'
                @test convert(Char, DNA_T) === Char(DNA_T) === 'T'
                @test convert(Char, DNA_N) === Char(DNA_N) === 'N'
            end

            @testset "RNA conversions to Char" begin
                @test convert(Char, RNA_A) === Char(RNA_A) === 'A'
                @test convert(Char, RNA_C) === Char(RNA_C) === 'C'
                @test convert(Char, RNA_G) === Char(RNA_G) === 'G'
                @test convert(Char, RNA_U) === Char(RNA_U) === 'U'
                @test convert(Char, RNA_N) === Char(RNA_N) === 'N'
            end
        end

        @testset "Other numeric types" begin
            @test convert(Int, DNA_A) === 1
            @test convert(Int, DNA_C) === 2
            @test convert(Int, DNA_G) === 4
            @test convert(Int, DNA_T) === 8
            @test convert(Int, DNA_N) === 15
            @test convert(DNA,  1) === DNA_A
            @test convert(DNA,  2) === DNA_C
            @test convert(DNA,  4) === DNA_G
            @test convert(DNA,  8) === DNA_T
            @test convert(DNA, 15) === DNA_N

            @test convert(Int, RNA_A) === 1
            @test convert(Int, RNA_C) === 2
            @test convert(Int, RNA_G) === 4
            @test convert(Int, RNA_U) === 8
            @test convert(Int, RNA_N) === 15
            @test convert(RNA,  1) === RNA_A
            @test convert(RNA,  2) === RNA_C
            @test convert(RNA,  4) === RNA_G
            @test convert(RNA,  8) === RNA_U
            @test convert(RNA, 15) === RNA_N
        end
        
        @testset "Nucleic acid types" begin
            fromto = [(DNA_Gap, RNA_Gap), (DNA_A, RNA_A), (DNA_C, RNA_C),
                      (DNA_M, RNA_M), (DNA_G, RNA_G), (DNA_R, RNA_R),
                      (DNA_S, RNA_S), (DNA_V, RNA_V), (DNA_T, RNA_U),
                      (DNA_W, RNA_W), (DNA_Y, RNA_Y), (DNA_H, RNA_H),
                      (DNA_K, RNA_K), (DNA_D, RNA_D), (DNA_B, RNA_B), (DNA_N, RNA_N)]
            
            for (from, to) in fromto
                @test convert(RNA, from) === RNA(from) === to
                @test convert(DNA, to) === DNA(to) === from
            end
        end
    end

    @testset "iscompatible" begin
        @test  iscompatible(DNA_A, DNA_A)
        @test  iscompatible(DNA_A, DNA_R)
        @test !iscompatible(DNA_C, DNA_A)
        @test !iscompatible(DNA_C, DNA_R)

        for x in alphabet(DNA)
            @test iscompatible(x, DNA_N) == (x != DNA_Gap)
            @test iscompatible(DNA_N, x) == (x != DNA_Gap)
        end

        @test  iscompatible(RNA_A, RNA_A)
        @test  iscompatible(RNA_A, RNA_R)
        @test !iscompatible(RNA_C, RNA_A)
        @test !iscompatible(RNA_C, RNA_R)

        for x in alphabet(RNA)
            @test iscompatible(x, RNA_N) == (x != RNA_Gap)
            @test iscompatible(RNA_N, x) == (x != RNA_Gap)
        end
    end

    @testset "isambiguous" begin
        for nt in alphabet(DNA)
            @test isambiguous(nt) == (nt ∉ (DNA_A, DNA_C, DNA_G, DNA_T, DNA_Gap))
        end
        for nt in alphabet(RNA)
            @test isambiguous(nt) == (nt ∉ (RNA_A, RNA_C, RNA_G, RNA_U, RNA_Gap))
        end
    end

    @testset "isGC" begin
        for nt in alphabet(DNA)
            @test isGC(nt) == (nt ∈ (DNA_G, DNA_C, DNA_S))
        end
        for nt in alphabet(RNA)
            @test isGC(nt) == (nt ∈ (RNA_G, RNA_C, RNA_S))
        end
    end

    @testset "ispurine" begin
        for nt in alphabet(DNA)
            @test ispurine(nt) == (nt == DNA_A || nt == DNA_G || nt == DNA_R)
        end
        for nt in alphabet(RNA)
            @test ispurine(nt) == (nt == RNA_A || nt == RNA_G || nt == RNA_R)
        end
    end

    @testset "ispyrimidine" begin
        for nt in alphabet(DNA)
            @test ispyrimidine(nt) == (nt == DNA_T || nt == DNA_C || nt == DNA_Y)
        end
        for nt in alphabet(RNA)
            @test ispyrimidine(nt) == (nt == RNA_U || nt == RNA_C || nt == RNA_Y)
        end
    end

    @testset "iscertain" begin
        for nt in alphabet(DNA)
            @test iscertain(nt) == (nt ∈ ACGT)
        end
        for nt in alphabet(RNA)
            @test iscertain(nt) == (nt ∈ ACGU)
        end
    end

    @testset "isgap" begin
        for nt in alphabet(DNA)
            @test isgap(nt) == (nt === DNA_Gap)
        end
        for nt in alphabet(RNA)
            @test isgap(nt) == (nt === RNA_Gap)
        end
    end
    
    @testset "isterm" begin
        for nt in alphabet(DNA)
            @test BioSymbols.isterm(nt) === false
        end
        for nt in alphabet(RNA)
            @test BioSymbols.isterm(nt) === false
        end
    end

    @testset "complement" begin
        @test complement(DNA_A) === DNA_T
        @test complement(DNA_C) === DNA_G
        @test complement(DNA_G) === DNA_C
        @test complement(DNA_T) === DNA_A
        @test complement(DNA_Gap) === DNA_Gap
        @test complement(DNA_N) === DNA_N

        @test complement(RNA_A) === RNA_U
        @test complement(RNA_C) === RNA_G
        @test complement(RNA_G) === RNA_C
        @test complement(RNA_U) === RNA_A
        @test complement(RNA_Gap) === RNA_Gap
        @test complement(RNA_N) === RNA_N
    end

    @testset "Arithmetic and Order" begin
        @testset "DNA" begin
            @test ~DNA_Gap === DNA_N
            @test ~DNA_N   === DNA_Gap
            @test DNA_A | DNA_C === DNA_M
            @test DNA_A & DNA_C === DNA_Gap
            @test DNA_Gap - DNA_A   === -1
            @test DNA_A   - DNA_Gap === +1
            @test DNA_Gap + 1 === DNA_Gap + 17 === DNA_A
            @test DNA_Gap - 1 === DNA_Gap - 17 === DNA_N
            @test DNA_Gap < DNA_A < DNA_C < DNA_G < DNA_T < DNA_N
            @test !(DNA_A > DNA_G)
            @test trailing_zeros(DNA_A) === 0
            @test trailing_zeros(DNA_C) === 1
            @test trailing_zeros(DNA_G) === 2
            @test trailing_zeros(DNA_T) === 3
            @test gap(DNA) === DNA_Gap
            @test collect(alphabet(DNA)) == sort([
                DNA_A, DNA_C, DNA_G, DNA_T,
                DNA_M, DNA_R, DNA_W, DNA_S,
                DNA_Y, DNA_K, DNA_V, DNA_H,
                DNA_D, DNA_B, DNA_N, DNA_Gap])
        end
        @testset "RNA" begin
            @test ~RNA_Gap === RNA_N
            @test ~RNA_N   === RNA_Gap
            @test RNA_A | RNA_C === RNA_M
            @test RNA_A & RNA_C === RNA_Gap
            @test RNA_Gap - RNA_A   === -1
            @test RNA_A   - RNA_Gap === +1
            @test RNA_Gap + 1 === RNA_Gap + 17 === RNA_A
            @test RNA_Gap - 1 === RNA_Gap - 17 === RNA_N
            @test RNA_Gap < RNA_A < RNA_C < RNA_G < RNA_U < RNA_N
            @test !(RNA_A > RNA_G)
            @test trailing_zeros(RNA_A) === 0
            @test trailing_zeros(RNA_C) === 1
            @test trailing_zeros(RNA_G) === 2
            @test trailing_zeros(RNA_U) === 3
            @test gap(RNA) === RNA_Gap
            @test collect(alphabet(RNA)) == sort([
                RNA_A, RNA_C, RNA_G, RNA_U,
                RNA_M, RNA_R, RNA_W, RNA_S,
                RNA_Y, RNA_K, RNA_V, RNA_H,
                RNA_D, RNA_B, RNA_N, RNA_Gap])
        end
    end

    @testset "Show DNA" begin
        dnas = [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N, DNA_Gap]
        @testset "print" begin
            buf = IOBuffer()
            for nt in dnas
                print(buf, nt)
            end
            @test String(take!(buf)) == "ACGTN-"
            @test_throws ArgumentError print(reinterpret(DNA, 0xf0))
        end

        @testset "show" begin
            buf = IOBuffer()
            for nt in dnas
                @test BioSymbols.prefix(nt) === BioSymbols.type_text(nt) === "DNA"
            end
            for nt in dnas
                show(buf, nt)
                write(buf, ' ')
            end
            @test String(take!(buf)) == "DNA_A DNA_C DNA_G DNA_T DNA_N DNA_Gap "
            @test sprint(show, reinterpret(DNA, 0xf0)) == "Invalid DNA"
        end
    end

    @testset "Show RNA" begin
        rnas = [RNA_A, RNA_C, RNA_G, RNA_U, RNA_N, RNA_Gap]
        @testset "print" begin
            buf = IOBuffer()
            for nt in rnas
                print(buf, nt)
            end
            @test String(take!(buf)) == "ACGUN-"
            @test_throws ArgumentError print(reinterpret(RNA, 0xf0))
        end

        @testset "show" begin
            for nt in rnas
                @test BioSymbols.prefix(nt) === BioSymbols.type_text(nt) === "RNA"
            end
            buf = IOBuffer()
            for nt in rnas
                show(buf, nt)
                write(buf, ' ')
            end
            @test String(take!(buf)) == "RNA_A RNA_C RNA_G RNA_U RNA_N RNA_Gap "
            @test sprint(show, reinterpret(RNA, 0xf0)) == "Invalid RNA"
        end
        
        @testset "read and write" begin
            @test round_trip(AA_X)
            @test round_trip(AA_Y)
        end
    end

    @testset "Sets" begin
        @test length(ACGT) == 4
        @test ACGT[1] === DNA_A
        @test ACGT[2] === DNA_C
        @test ACGT[3] === DNA_G
        @test ACGT[4] === DNA_T
        @test collect(ACGT) == [DNA_A, DNA_C, DNA_G, DNA_T]

        @test length(ACGTN) == 5
        @test ACGTN[1] === DNA_A
        @test ACGTN[2] === DNA_C
        @test ACGTN[3] === DNA_G
        @test ACGTN[4] === DNA_T
        @test ACGTN[5] === DNA_N
        @test collect(ACGTN) == [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]

        @test length(ACGU) == 4
        @test ACGU[1] === RNA_A
        @test ACGU[2] === RNA_C
        @test ACGU[3] === RNA_G
        @test ACGU[4] === RNA_U
        @test collect(ACGU) == [RNA_A, RNA_C, RNA_G, RNA_U]

        @test length(ACGUN) == 5
        @test ACGUN[1] === RNA_A
        @test ACGUN[2] === RNA_C
        @test ACGUN[3] === RNA_G
        @test ACGUN[4] === RNA_U
        @test ACGUN[5] === RNA_N
        @test collect(ACGUN) == [RNA_A, RNA_C, RNA_G, RNA_U, RNA_N]
    end
end

@testset "Aminoacids" begin
    @testset "conversion" begin
        @test BioSymbols.bytemask(AA_A) === 0b11111
        
        for (int, aa) in [
            (0x00, AA_A), (0x01, AA_R), (0x02, AA_N), (0x03, AA_D), (0x04, AA_C),
            (0x05, AA_Q), (0x06, AA_E), (0x07, AA_G), (0x08, AA_H), (0x09, AA_I),
            (0x0a, AA_L), (0x0b, AA_K), (0x0c, AA_M), (0x0d, AA_F), (0x0e, AA_P),
            (0x0f, AA_S), (0x10, AA_T), (0x11, AA_W), (0x12, AA_Y), (0x13, AA_V),
            (0x14, AA_O), (0x15, AA_U), (0x16, AA_B), (0x17, AA_J), (0x18, AA_Z),
            (0x19, AA_X), (0x1a, AA_Term), (0x1b, AA_Gap)]
            @test convert(AminoAcid, int) === AminoAcid(int) === aa
        end

        @test convert(AminoAcid, 0) === AA_A
        @test convert(AminoAcid, 10) === AA_L

        for (c, aa) in [
                ('A', AA_A), ('R', AA_R), ('N', AA_N), ('D', AA_D), ('C', AA_C),
                ('Q', AA_Q), ('E', AA_E), ('G', AA_G), ('H', AA_H), ('I', AA_I),
                ('L', AA_L), ('K', AA_K), ('M', AA_M), ('F', AA_F), ('P', AA_P),
                ('S', AA_S), ('T', AA_T), ('W', AA_W), ('Y', AA_Y), ('V', AA_V),
                ('O', AA_O), ('U', AA_U), ('B', AA_B), ('J', AA_J), ('Z', AA_Z),
                ('X', AA_X), ('*', AA_Term), ('-', AA_Gap)]
            @test convert(AminoAcid, c) === convert(AminoAcid, lowercase(c)) == AminoAcid(c) === aa
            @test Char(aa) === c
        end
        @test_throws InexactError convert(AminoAcid, '\0')
        @test_throws InexactError convert(AminoAcid, '@')
        @test_throws InexactError convert(AminoAcid, '亜')
    end

    @testset "isvalid" begin
        for aa in alphabet(AminoAcid)
            @test isvalid(aa)
        end
        @test !isvalid(reinterpret(AminoAcid, 0x1c))
        @test !isvalid(reinterpret(AminoAcid, 0xff))
        @test  isvalid(AminoAcid, 0x1b)
        @test !isvalid(AminoAcid, 0x1c)
    end

    @testset "Arithmetic and Order" begin
        @test AA_A + 1 == 1 + AA_A == AA_R
        @test AA_R + 1 == 1 + AA_R == AA_N
        @test AA_A + 2 == 2 + AA_A == AA_N
        @test AA_A + 28 == 28 + AA_A == AA_A
        @test AA_R - 1 == AA_A
        @test AA_N - 2 == AA_A
        @test AA_A - 28 == AA_A
        @test AA_D - AA_A ==  3
        @test AA_A - AA_D == -3
        @test (AA_A < AA_R < AA_N < AA_V < AA_O < AA_U <
               AA_B < AA_J < AA_Z < AA_X < AA_Term < AA_Gap)
        @test !(AA_J < AA_B)

        @test gap(AminoAcid) === AA_Gap
        @test length(alphabet(AminoAcid)) == 28
        @test AA_A in alphabet(AminoAcid)
        @test AA_I in alphabet(AminoAcid)
        @test AA_U in alphabet(AminoAcid)
    end

    @testset "Range" begin
        @test !(AA_C in AA_Q:AA_H)
        @test   AA_Q in AA_Q:AA_H
        @test   AA_E in AA_Q:AA_H
        @test   AA_G in AA_Q:AA_H
        @test   AA_H in AA_Q:AA_H
        @test !(AA_I in AA_Q:AA_H)

        @test collect(AA_W:AA_V) == [AA_W, AA_Y, AA_V]
    end

    @testset "iscompatible" begin
        @test  iscompatible(AA_A, AA_A)
        @test !iscompatible(AA_A, AA_R)

        for x in alphabet(AminoAcid)
            @test iscompatible(x, AA_B) == (x ∈ (AA_D, AA_N, AA_B, AA_X))
            @test iscompatible(x, AA_J) == (x ∈ (AA_I, AA_L, AA_J, AA_X))
            @test iscompatible(x, AA_Z) == (x ∈ (AA_E, AA_Q, AA_Z, AA_X))
            @test iscompatible(x, AA_X) == (x ∉ (AA_Term, AA_Gap))
        end
    end

    @testset "isambiguous" begin
        for x in alphabet(AminoAcid)
            @test isambiguous(x) == (AA_B <= x <= AA_X)
        end
    end

    @testset "iscertain" begin
        for x in alphabet(AminoAcid)
            @test iscertain(x) == (!isambiguous(x) && x != AA_Gap)
        end
    end

    @testset "isgap" begin
        for x in alphabet(AminoAcid)
            @test isgap(x) == (x == AA_Gap)
        end
    end

    @testset "Show amino acid" begin
        aas = [AA_A, AA_D, AA_B, AA_X, AA_Term, AA_Gap]
        @testset "print" begin
            buf = IOBuffer()
            for aa in aas
                print(buf, aa)
            end
            @test String(take!(buf)) == "ADBX*-"
            @test_throws ArgumentError print(BioSymbols.AA_INVALID)
        end

        @testset "show" begin
            for aa in aas
                @test BioSymbols.prefix(aa) === "AA"
                @test BioSymbols.type_text(aa) === "Amino Acid"
            end
            buf = IOBuffer()
            for aa in aas
                show(buf, aa)
                write(buf, ' ')
            end
            @test String(take!(buf)) == "AA_A AA_D AA_B AA_X AA_Term AA_Gap "
            @test sprint(show, BioSymbols.AA_INVALID) == "Invalid Amino Acid"
        end
        
        @testset "read and write" begin
            @test round_trip(DNA_G)
            @test round_trip(RNA_U)
        end
    end

    @testset "Parsers" begin
        @testset "Valid Cases" begin
            # case-insensitive and ignores spaces
            @test parse(AminoAcid, "a") == AA_A
            @test parse(AminoAcid, "Ala") == AA_A
            @test parse(AminoAcid, "aLa ") == AA_A
            @test parse(AminoAcid, " alA ") == AA_A
            @test parse(AminoAcid, "\tAlA\n") == AA_A
            @test parse(AminoAcid, "x") == AA_X
            @test parse(AminoAcid, "X") == AA_X
            aas = [
                ("A", "ALA", AA_A),
                ("R", "ARG", AA_R),
                ("N", "ASN", AA_N),
                ("D", "ASP", AA_D),
                ("C", "CYS", AA_C),
                ("E", "GLU", AA_E),
                ("Q", "GLN", AA_Q),
                ("G", "GLY", AA_G),
                ("H", "HIS", AA_H),
                ("I", "ILE", AA_I),
                ("L", "LEU", AA_L),
                ("K", "LYS", AA_K),
                ("M", "MET", AA_M),
                ("F", "PHE", AA_F),
                ("P", "PRO", AA_P),
                ("S", "SER", AA_S),
                ("T", "THR", AA_T),
                ("W", "TRP", AA_W),
                ("Y", "TYR", AA_Y),
                ("V", "VAL", AA_V),
                ("O", "PYL", AA_O),
                ("U", "SEC", AA_U),
                ("B", "ASX", AA_B),
                ("J", "XLE", AA_J),
                ("Z", "GLX", AA_Z),
                ("X", "XAA", AA_X),
            ]
            @test length(aas) == 26
            for (one, three, aa) in aas
                @test parse(AminoAcid, one) == aa
                @test parse(AminoAcid, three) == aa
                @test parse(AminoAcid, Char(one[1])) == aa
                @test tryparse(AminoAcid, one) === aa
                @test tryparse(AminoAcid, three) === aa
                @test tryparse(AminoAcid, Char(one[1])) === aa
            end
            @test parse(AminoAcid, "*") == AA_Term
            @test parse(AminoAcid, "-") == AA_Gap
            @test tryparse(AminoAcid, "*") === AA_Term
            @test tryparse(AminoAcid, "-") === AA_Gap
        end

        @testset "Invalid Cases" begin
            @test_throws ArgumentError parse(AminoAcid, "")
            @test_throws ArgumentError parse(AminoAcid, "AL")
            @test_throws ArgumentError parse(AminoAcid, "LA")
            @test_throws ArgumentError parse(AminoAcid, "ALAA")
            @test_throws ArgumentError parse(AminoAcid, '\0')
            @test_throws ArgumentError parse(AminoAcid, '@')
            @test_throws ArgumentError parse(AminoAcid, '亜')
            @test tryparse(AminoAcid, "") == nothing
            @test tryparse(AminoAcid, "AL") == nothing
            @test tryparse(AminoAcid, "LA") == nothing
            @test tryparse(AminoAcid, "ALAA") == nothing
            @test tryparse(AminoAcid, '\0') == nothing
            @test tryparse(AminoAcid, '@') == nothing
            @test tryparse(AminoAcid, '亜') == nothing
        end
    end
end
