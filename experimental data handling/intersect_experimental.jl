using CSV
using DataFrames
using LinearAlgebra
using DelimitedFiles
using LightGraphs
using Statistics
using GraphPlot
using SpecialFunctions
using DataFrames


datasets=["... put the name of the reference database containing the intersections here"]

for set in datasets
    ## Set input
    gene_expression = readdlm("... put the pathe to the expression data here",',')
    ref_int = string("... put the path to the folder containing the reference interactions here",set,".txt")
    ref_interactions = readdlm(ref_int,'\t')
    network_name = string("... put the path to the output folder here",set,".csv")

    a1 = ref_interactions[:,1]
    a2 = ref_interactions[:,2]
    exp_genes = gene_expression[:,1]

    ## Perform he intersection
    open(network_name, "w") do io
        print(io, "Gene1,Gene2")
        for i in 1:size(ref_interactions)[1]
            if (a1[i] in exp_genes) && (a2[i] in exp_genes)
                intersection_name = string(string(a1[i]),",",string(a2[i]))
                print(io, "\n")
                print(io, intersection_name)
            end
        end
    end
    ## Compute adjacency matrixes
    ref_genes_unique = unique(vcat(unique(a1),unique(a2)))
    A = zeros(length(ref_genes_unique),length(ref_genes_unique))

    for i in 1:size(ref_interactions)[1]
        if (a1[i] in ref_genes_unique) && (a2[i] in ref_genes_unique)
            pos1 = findall(x->x==a1[i],ref_genes_unique)
            pos2 = findall(x->x==a2[i],ref_genes_unique)
            A[pos1[1],pos2[1]]=1
            A[pos2[1],pos1[1]]=1
        end
    end

    out_name = string("... put the path to the output folder here",set,"_adjacency.csv")
    writedlm(out_name,A,',')
end
