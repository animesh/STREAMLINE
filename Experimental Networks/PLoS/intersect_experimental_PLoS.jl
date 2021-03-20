using CSV
using DataFrames
using LinearAlgebra
using DelimitedFiles
using Polyhedra
using NumericalIntegration
using StatsPlots
using LightGraphs
using Statistics
using GraphPlot
using SpecialFunctions
using GLM
using DataFrames
using CSV

Exp = readdlm("C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Data/PLoS/C_elegans_expression_normalized.txt")
#Exp = DataFrame(CSV.File("C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Data/PLoS/A_thaliana_expression_notnormalized.csv"))
inter = readdlm("C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Data/PLoS/C_elegans_interactions.txt",String)
inter2 = inter

for i in 1:size(inter)[1]
    inter2[i,1] = string("ENSMUSG000",inter[i,1])
    inter2[i,2] = string("ENSMUSG000",inter[i,2])
end


exp_genes = Exp[2:end,1]
inter_genes = intersect(unique(inter2[:,1]),unique(inter2[:,2]))
genes = intersect(inter_genes,exp_genes)
writedlm("Cthaliana_genes.txt",genes)

tfs = genes
tf_network = inter

s=length(tfs)

degrees_in=zeros(length(tfs))
degrees_out=zeros(length(tfs))
for i in 1:length(tfs)
    degrees_in[i]=length(findall(x->x==tfs[i],tf_network[:,2]))
    degrees_out[i]=length(findall(x->x==tfs[i],tf_network[:,1]))
end


mA=zeros(s,s)

for i in 1:s
    if degrees_out[i] != 0
        connection_ids = findall(x->x==tfs[i],tf_network[:,1])
        target_genes = tf_network[connection_ids,2]
        for c in target_genes
            target_id = findall(x->x==c,tfs)
            if length(target_id)==1
                mA[i,target_id[1]] = 1
            end
        end
    end
end

sum(mA)

writedlm("Celegans_adjcency.txt",mA)
