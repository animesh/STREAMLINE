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

# mouse_tfs = readdlm("expdata/mouse-tfs.csv",',')[:,1]
# mESC_ordering = readdlm("expdata/GeneOrdering.csv",',')[2:end,:]
# mESC_ChIP = readdlm("expdata/mESC-ChIP-seq-network.csv",',')
# mESC_lofgof = readdlm("expdata/mESC-lofgof-network.csv",',')
#
# degrees_in=zeros(length(mouse_tfs))
# degrees_out=zeros(length(mouse_tfs))
# for i in 1:1480
#     degrees_in[i]=length(findall(x->x==mouse_tfs[i],mESC_ChIP[:,2]))
#     degrees_out[i]=length(findall(x->x==mouse_tfs[i],mESC_ChIP[:,1]))
# end
#
#
# mA=zeros(1480,1480)
#
# for i in 1:1480
#     if degrees_out[i] != 0
#         connection_ids = findall(x->x==mouse_tfs[i],mESC_ChIP[:,1])
#         target_genes = mESC_ChIP[connection_ids,2]
#         print("\n")
#         print(target_genes)
#         for c in target_genes
#             target_id = findall(x->x==c,mouse_tfs)
#             if length(target_id)==1
#                 mA[i,target_id[1]] = 1
#             end
#         end
#     end
# end

m_degrees = readdlm("mESC_ChIP_degrees.txt")
tf_mouse = h_degrees[:,1]
m_in_degrees = h_degrees[:,2]
m_out_degrees = h_degrees[:,3]

### Estimate distribution and plot histogram
m_dist_est = zeros(50)
for i in 1:50
    m_dist_est[i]=sum(degrees_out[degrees_out.>0][degrees_out[degrees_out.>0].<360*i].>360*(i-1))
end
mESCbar = bar(360*[1:50],m_dist_est,
    label="degrees grouped by 360",
    xlabel="Degree",
    ylabel="Number of TF with deg>0",
    title="mESC network out-degree distribution")
savefig(mESCbar,"images/mESC_degree_hist.png")

### Plot log-log distribution
hESCscatter = scatter(log.(360*vec(1:50)),log.(m_dist_est),
        label="degrees grouped by 360",
        xlabel="Degree",
        ylabel="Number of TF with deg>0",
        title="mESC network out-degree distribution")
data_m = DataFrame(X=log.(360*vec(1:50)), Y=log.(m_dist_est))
ols_m = lm(@formula(Y ~ X), data_m)
slope = -1
intercept = 7.5
ms=5.8:0.001:9.8
plot!(hESCscatter,ms,slope*ms.+intercept, color="red", label="linear regression fit (significant)")
savefig(mESCscatter,"images/mESC_degree_scatter.png")

# mESC_ChIP_network=SimpleDiGraph(mA)
# global_clustering_coefficient(mESC_ChIP_network)
# diameter(mESC_ChIP_network)
# connected_components(mESC_ChIP_network)
# weakly_connected_components(mESC_ChIP_network)[2]
