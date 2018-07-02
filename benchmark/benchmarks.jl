using BenchmarkTools, BioSymbols

sorted_DNA_symbols = alphabet(DNA)
n_nucleic_acids = length(sorted_DNA_symbols)

function symbol_combos(alphabet)
    combos = eltype(alphabet)
    for i = 1:length(alphabet)
        for j = (i + 1):length(alphabet)
            
        end
    end
end
    





SUITE = BenchmarkGroup()
SUITE["NucleicAcids"] = BenchmarkGroup()

SUITE["NucleicAcids"]["iscompatible"] = BenchmarkGroup()
SUITE["NucleicAcids"]["iscompatible"]["DNA"] = @benchmarkable iscompatible(DNA_A, DNA_R)
SUITE["NucleicAcids"]["iscompatible"]["RNA"] = @benchmarkable iscompatible(RNA_A, RNA_R)

SUITE["NucleicAcids"]["isambiguous"] = BenchmarkGroup()

@benchmarkable isambiguous(DNA_M)



SUITE["NucleicAcids"]["isGC"] = @benchmarkable isGC(DNA_G)
SUITE["NucleicAcids"]["ispurine"] = @benchmarkable ispurine(DNA_A)
SUITE["NucleicAcids"]["iscertain"] = @benchmarkable ispyrimidine(DNA_A)
SUITE["NucleicAcids"]["isgap"] = @benchmarkable isgap(DNA_Gap)
SUITE["NucleicAcids"]["complement"] = @benchmarkable complement(DNA_A)
SUITE["NucleicAcids"]["arithmetic"] = @benchmarkable ispyrimidine(DNA_A)
SUITE["NucleicAcids"]["ispyrimidine"] = @benchmarkable ispyrimidine(DNA_A)
