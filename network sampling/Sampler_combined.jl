using LightGraphs
using Plots
#using GraphPlot
using Random
using Distributions
using DelimitedFiles

function ER_sample(N_nodes,N_edges,M)
    for u in 1:M
        for n in N_nodes
            for k in N_edges
                p=k/(n*(n-1))
                network=erdos_renyi(n, p, is_directed=true)
                name = string("samples/ER_",string(n),"nodes_",string(k),"edges_",string(u),".txt")
                A=adjacency_matrix(network)
                writedlm(name,Matrix(A))
            end
        end
    end
end

function scale_free_sample(N_nodes,N_edges,α_ins::Array{Float64,1},α_out::Array{Float64,1},M)
    u=1
    while u<=M
        for n in N_nodes
            for k in N_edges
                for ai in α_ins
                    for ao in α_out
                        network = static_scale_free(n, k, ao, ai)
                        maxdeg = maximum(indegree(network))
                        mindeg = minimum(indegree(network))
                        name = string("samples/scale_free_",string(n),"nodes_",string(k),"edges_",string(ao),"outexp_",string(ai),"inexp_",string(u),".txt")
                        A=adjacency_matrix(network)
                        if mindeg>0 && maxdeg<=16
                            writedlm(name,Matrix(A))
                            u+=1
                        end
                    end
                end
            end
        end
    end
end

function small_world_sample(N_nodes,N_sdegrees,betas,M)
    for u in 1:M
        for n in N_nodes
            for k in N_sdegrees
                for b in betas
                    network = watts_strogatz(n, Int(k), b)
                    name = string("samples/smallworld_",string(n),"nodes_",string(k),"sdegree_",string(b),"beta_",string(u),".txt")
                    A=adjacency_matrix(network)
                    writedlm(name,Matrix(A))
                end
            end
        end
    end
end

function semi_scale_free(N_nodes,N_edges,α_out,t,M)
    for u in 1:M
        for n in N_nodes
            for k in N_edges
                for ao in α_out
                    network=static_scale_free(n, Int.(ceil(k/t)), ao, Inf)
                    A=adjacency_matrix(network)
                    delete_nodes=setdiff(1:n,rand(1:n,Int(ceil(n*t))))
                    A[delete_nodes,:].=0
                    name = string("samples/semi_scale_free_",string(n),"nodes_",string(k),"edges_",string(ao),"outexp_",string(t),"inper_",string(u),".txt")
                    writedlm(name,Matrix(A))
                end
            end
        end
    end
end


####### Parameter ranges:
# Set them as desired
# N_nodes = 
# N_edges = 

### Scale-free
# α_ins= 
# α_out= 2:.5:10

### Small-world
# betas = 
# N_sdegrees = 

### Using the intervals given above
#ER_sample(,,)
#scale_free_sample(N_nodes,N_edges,α_ins,α_out)
#semi_scale_free(,,,,)
#small_world_sample(,,,)
#scale_free_sample(,,,,)
