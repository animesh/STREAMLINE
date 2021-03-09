using LightGraphs
using Plots
using GraphPlot
using Random
using Distributions



function ER_textfiles(M,N,p1,p2)
    for k in 1:M
        # set network name
        network_name_BO = string("ER_BoolODE_",string(k),".txt")
        network_name_BE = string("ER_BEELINE_",string(k),".txt")

        #matrix storing the positive connections
        B=zeros(N,N)

        # sample network
        network=erdos_renyi(N, p1, is_directed=true)
        A=adjacency_matrix(network)

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

        # write txt file for BEELINE
        open(network_name_BE, "w") do jo
            print(jo, "Gene1,Gene2,Type")
            for i in 1:N
                gene_name_2 = string("g",string(i))
                for j in 1:N
                    if A[j,i]!=0
                        print(jo, "\n")
                        gene_name_1 = string("g",string(j))
                        if B[j,i]==1
                            print(jo,gene_name_1,",",gene_name_2,",","+")
                        else
                            print(jo,gene_name_1,",",gene_name_2,",","-")
                        end
                    end
                end
            end
        end

    end
end

ER_textfiles(5,100,.1,.9)
