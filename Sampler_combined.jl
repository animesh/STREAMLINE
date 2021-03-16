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

function scale_free_sample(N_nodes=Array{Any},N_edges=Array{Any},α_ins=Array{Any},α_out=Array{Any})
    for n in N_nodes
        for k in N_edges
            for ao in α_ins
                for ai in α_out
                    network=static_scale_free(n, k, ao, ai)
                    name = string("samples/scale_free_",string(n),"nodes_",string(k),"edges_",string(ao),"outexp_",string(ai),"inexp.txt")
                    A=adjacency_matrix(network)
                    writedlm(name,Matrix(A))
                end
            end
        end
    end
end

function small_world_sample(N_nodes=Array{Any},N_edges=Array{Any},betas=Array{Any})
    for n in N_nodes
        for k in N_edges
            for b in betas
                network = watts_strogatz(n, k, b)
                name = string("samples/smallworld_",string(n),"nodes_",string(k),"degree_",string(b),"beta.txt")
                A=adjacency_matrix(network)
                writedlm(name,Matrix(A))
            end
        end
    end
end
