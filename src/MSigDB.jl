module MSigDB

#using RCall
using CSV
using DataFrames
using NoLongerProblems

function HumanEntrezEnsemblID()
    folder = "MSigDB"
    file = normpath(folder, "Human_EnsemblID_.csv")
    tb = DataFrames.DataFrame(DelimitedFiles.readdlm(file)[2:end, 2:end])
    rename!(tb, :x1 => :EnsemblID_human, :x2 => :Entrez)
end

function mouse_human_ortologs()
    ortologs = "MSigDB/mart_mouseortologstohuman.csv"
    return DataFrame!(CSV.File(ortologs))
end



function get_MSigDB(pat; mouseGeneSymbol = true)
    MSigDB_folder = "MSigDB"
    matrix = DelimitedFiles.readdlm(normpath(MSigDB_folder, pat))
    
    new_dict = Dict{String, Array}()
    
    for ii in 1:size(matrix)[1]
        key = matrix[ii, 1] 
        bool = [isa(jj, Int) for jj in matrix[ii, :]]
        valu = matrix[ii, bool]
        new_dict[key] = valu
    end
    
    if mouseGeneSymbol
        
         hinf = HumanEntrezEnsemblID(); hinf = hinf[hinf[!,:Entrez].!="NA", :]
         hinf
        
         mort = mouse_human_ortologs(); rename!(mort, 
            Symbol("Gene stable ID") => :EnsemblID_human, 
            Symbol("Mouse gene stable ID") => :EnsemblID,
            Symbol("Mouse orthology confidence [0 low, 1 high]") => :Ortology,
            Symbol("Mouse gene name") => :GeneSymbol,)
        
         mort = innerjoin(hinf, mort, on =:EnsemblID_human)
        # Only genes with high Ortology confidence and Ortology ortholog_one2one       
        mort =  mort[mort[!,:Ortology].==1, :]
        mort =  mort[mort[!,Symbol("Mouse homology type")].=="ortholog_one2one", :]
        mort_dict = Dict(mort[!,:Entrez], mort[!,:GeneSymbol])
        # Translate Entrez to gene symbols
        
        new_new_dict = Dict()
        
        for key in keys(new_dict)
            new_new_dict[key] = []
            for ii in  new_dict[key]
                try 
                    push!(new_new_dict[key], mort_dict[ii])
                catch;
                end
            end
            
        end
        
        return new_new_dict
      
    end
     return new_dict
    
    
end

function HALLMARK_OXIDATIVE_PHOSPHORYLATION()
    
    hallmarkdatabase = MSigDB.get_MSigDB("h.all.v6.2.entrez.gmt", mouseGeneSymbol = true)
    df = DataFrames.DataFrame(GeneSymbol = hallmarkdatabase["HALLMARK_OXIDATIVE_PHOSPHORYLATION"])
end


end