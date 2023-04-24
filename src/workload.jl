using PrecompileTools

@compile_workload begin
    for sym in [DNA_W, RNA_U, AA_Y]
        print(IOBuffer(), sym)
        Char(sym)

        alphabet(typeof(sym))
        encoded_data(sym)

        isgap(sym)
        iscertain(sym)
        isambiguous(sym)
    end

    for sym in [DNA_W, RNA_U]
        isGC(sym)
        ispurine(sym)
        ispyrimidine(sym)
        complement(sym)
    end
end