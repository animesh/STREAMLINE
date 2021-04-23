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
using FreqTables

# human_tfs = readdlm("expdata/human-tfs.csv",',')[:,1]
# hESC_ChIP = readdlm("expdata/hESC-ChIP-seq-network.csv",',')

### -----------------------------
### CALCULATE DEGREE DISTRIBUTION
### -----------------------------

# h_degrees_in=zeros(length(human_tfs))
# h_degrees_out=zeros(length(human_tfs))
# for i in 1:length(human_tfs)
#     h_degrees_in[i]=length(findall(x->x==human_tfs[i],hESC_ChIP[:,2]))
#     h_degrees_out[i]=length(findall(x->x==human_tfs[i],hESC_ChIP[:,1]))
# end

### -----------------------------
### CALCULATE ADJACENCY MATRIX
### -----------------------------

### IMPORTANT: A[i,j]=1 means there is an edge from gene i to gene j

# hA=zeros(1564,1564)
#
# for i in 1:1564
#     if h_degrees_out[i] != 0
#         connection_ids = findall(x->x==human_tfs[i],hESC_ChIP[:,1])
#         target_genes = hESC_ChIP[connection_ids,2]
#         for c in target_genes
#             target_id = findall(x->x==c,human_tfs)
#             if length(target_id)==1
#                 hA[i,target_id[1]] = 1
#             end
#         end
#     end
# end

hESC_ChIP_network=SimpleDiGraph(hA)
global_clustering_coefficient(hESC_ChIP_network)
diameter(hESC_ChIP_network)
connected_components(hESC_ChIP_network)

h_degrees = readdlm("hESC_ChIP_degrees.txt")
tf_human = h_degrees[:,1]
h_in_degrees = h_degrees[:,2]
h_out_degrees = h_degrees[:,3]

### Estimate distribution and plot histogram
h_dist_est = zeros(50)
for i in 1:50
    h_dist_est[i]=sum(h_out_degrees[h_out_degrees.>0][h_out_degrees[h_out_degrees.>0].<280*i].>280*(i-1))
end
hESCbar = bar(280*[1:50],h_dist_est,
    label="degrees grouped by 280",
    xlabel="Degree",
    ylabel="Number of TF with deg>0",
    title="hESC network out-degree distribution")
savefig(hESCbar,"images/hESC_degree_hist.png")


### Plot log-log distribution
hESCscatter = scatter(log.(280*vec(1:50)),log.(h_dist_est),
        label="degrees grouped by 360",
        xlabel="Degree",
        ylabel="log-Number of TF with deg>0",
        title="hESC network")
data_h = DataFrame(X=log.(280*vec(1:50)), Y=log.(h_dist_est))
ols_h = lm(@formula(Y ~ X), data_h)
slope = -.8
intercept = 7.4
ms=5.5:0.001:9.5
plot!(hESCscatter,ms,slope*ms.+intercept, color="red", label="linear regression fit (significant)")
savefig(hESCscatter,"images/hESC_degree_scatter.png")
