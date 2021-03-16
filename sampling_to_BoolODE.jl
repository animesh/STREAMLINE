using LightGraphs
using Plots
using GraphPlot
using Random
using Distributions

function BoolODE_textfiles(directory)
    txts=readdir(directory)
    for txt in txts
            # set network name
            network_name_BO = string("BoolODE_",txt)

            # read network
            A_path=string("samples/",txt)
            A=readdlm(A_path)
            N=size(A)[1]

            #matrix storing the positive connections
            B=zeros(N,N)
            p1=0
            p2=1

            # write txt file for BoolODE
            open(network_name_BO, "w") do io
                print(io, "Gene  Rule \n")
                for i in 1:N
                    gene_name = string("g",string(i))
                    print(io, gene_name)
                    print(io, "    ")
                    connections_pos=[]
                    connections_neg=[]
                    for j in 1:N
                        p=rand(1)[1]
                        if A[j,i]!=0
                            if p>p1
                                append!(connections_pos,j)
                                B[j,i]=1
                            else
                                append!(connections_neg,j)
                            end
                        end
                    end
                    if connections_pos==[]
                        if connections_neg==[]
                            print(io, "\n")
                        else
                            print(io,"not ")
                            print(io, "( ")
                            for g in 1:length(connections_neg)
                                if g-1>0
                                    print(io, "or ")
                                end
                                gene_name = string("g",connections_neg[g]," ")
                                print(io, gene_name)
                            end
                            print(io, ")")
                            print(io, "\n")
                        end
                    else
                        if connections_neg==[]
                            print(io, "( ")
                            for g in 1:length(connections_pos)
                                if g-1>0
                                    print(io, "or ")
                                end
                                gene_name = string("g",connections_pos[g]," ")
                                print(io, gene_name)
                            end
                            print(io, ")")
                            print(io, "\n")
                        else
                            print(io, "( ")
                            for g in 1:length(connections_pos)
                                if g-1>0
                                    print(io, "or ")
                                end
                                gene_name = string("g",connections_pos[g]," ")
                                print(io, gene_name)
                            end
                            print(io, ")")
                            print(io," and not ")
                            print(io, "( ")
                            for g in 1:length(connections_neg)
                                if g-1>0
                                    print(io, "or ")
                                end
                                gene_name = string("g",connections_neg[g]," ")
                                print(io, gene_name)
                            end
                            print(io, ")")
                            print(io, "\n")
                        end
                    end
                end
            end
    end
end
