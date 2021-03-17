using LightGraphs
using Plots
using GraphPlot
using Random
using Distributions

function ER_sample(N_nodes=Array{Any},N_edges=Array{Any})
    for n in N_nodes
        for k in N_edges
            p=k/(n*(n-1))
            network=erdos_renyi(n, p, is_directed=true)
            name = string("samples/ER_",string(n),"nodes_",string(k),"edges.txt")
            A=adjacency_matrix(network)
            writedlm(name,Matrix(A))
        end
    end
end

function scale_free_sample(N_nodes::Array{Any},N_edges::Array{Any},α_ins::Array{Any},α_out::Array{Any})
    for n in N_nodes
        for k in N_edges
            for ai in α_ins
                for ao in α_out
                    network=static_scale_free(n, k, ao, ai)
                    name = string("samples/scale_free_",string(n),"nodes_",string(k),"edges_",string(ao),"outexp_",string(ai),"inexp.txt")
                    A=adjacency_matrix(network)
                    writedlm(name,Matrix(A))
                end
            end
        end
    end
end

function small_world_sample(N_nodes::Array{Any},N_sdegrees::Array{Any},betas::Array{Any})
    for n in N_nodes
        for k in N_sdegrees
            for b in betas
                network = watts_strogatz(n, k, b)
                name = string("samples/smallworld_",string(n),"nodes_",string(k),"sdegree_",string(b),"beta.txt")
                A=adjacency_matrix(network)
                writedlm(name,Matrix(A))
            end
        end
    end
end

function semi_scale_free(N_nodes::Array{Int64,1},N_edges::Array{Int64,1},α_out::Array{Float64,1},t::Float64)
    for n in N_nodes
        for k in N_edges
            for ao in α_out
                network=static_scale_free(n, Int(ceil(k/t)), ao, Inf)
                A=adjacency_matrix(network)
                delete_nodes=setdiff(1:n,rand(1:n,Int(ceil(n*t))))
                A[delete_nodes,:].=0
                name = string("samples/semi_scale_free_",string(n),"nodes_",string(k),"edges_",string(ao),"outexp_",string(t),"inper.txt")
                writedlm(name,Matrix(A))
            end
        end
    end
end


####### Parameter ranges:
### General:
# N_nodes = 500:500:6000
# N_edges = 10000:10000:100000

### Scale-free
# α_ins= 2:.5:10
# α_out= 2:.5:10

### Scale-free
# betas = 0.1:0.1:0.9
# N_sdegrees = N_nodes.*[0.1:0.1:0.9]
