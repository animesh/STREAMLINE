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


tfs = readdlm("expdata/human-tfs.csv",',')[:,1]
tf_network = readdlm("C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Data/BEELINE-Networks/Networks/human/hESC-ChIP-seq-network.csv",',')

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


writedlm("expdata_networks_correct/hESC_ChIP.txt",mA)
