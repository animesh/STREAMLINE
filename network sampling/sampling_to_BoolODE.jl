using LightGraphs
using Plots
#using GraphPlot
using Random
using Distributions
using DelimitedFiles

function BoolODE_textfiles_old(directory)
    txts=readdir(directory)
    for txt in txts
            # set network name
            network_name_BO = string("BoolODE_",txt)

            # read network
            A_path=string(directory,txt)
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
                    if i in 1:9
                        print(io, "    ")
                    elseif i in 10:99
                        print(io, "   ")
                    elseif i in 100:999
                        print(io, "  ")
                    end
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
                        # else
                        #     print(io,"not ")
                        #     print(io, "( ")
                        #     for g in 1:length(connections_neg)
                        #         if g-1>0
                        #             print(io, "or ")
                        #         end
                        #         gene_name = string("g",connections_neg[g]," ")
                        #         print(io, gene_name)
                        #     end
                        #     print(io, ")")
                        #     print(io, "\n")
                        end
                    else
                        # if connections_neg==[]
                        #     print(io, "( ")
                        #     for g in 1:length(connections_pos)
                        #         if g-1>0
                        #             print(io, "or ")
                        #         end
                        #         gene_name = string("g",connections_pos[g]," ")
                        #         print(io, gene_name)
                        #     end
                        #     print(io, ")")
                        #     print(io, "\n")
                        # else
                            print(io, "( ")
                            for g in 1:length(connections_pos)
                                if g-1>0
                                    print(io, "or ")
                                end
                                gene_name = string("g",connections_pos[g]," ")
                                print(io, gene_name)
                            end
                            print(io, ")")
                        #     print(io," and not ")
                        #     print(io, "( ")
                        #     for g in 1:length(connections_neg)
                        #         if g-1>0
                        #             print(io, "or ")
                        #         end
                        #         gene_name = string("g",connections_neg[g]," ")
                        #         print(io, gene_name)
                        #     end
                        #     print(io, ")")
                        #     print(io, "\n")
                        end
                    end
                end
            end
    end
end

function BoolODE_textfiles(directory)
    txts=readdir(directory)
    for txt in txts
            # set network name
            network_name_BO = string("BoolODE_",txt)

            # read network
            A_path=string(directory,txt)
            A=readdlm(A_path)
            N=size(A)[1]


            # write txt file for BoolODE
            open(network_name_BO, "w") do io
                g = "Gene"
                r = "Rule"
                writedlm(io, [g r])
                for i in 1:N
                    gene_name = string("g",string(i))

                    connections_pos=[]
                    for j in 1:N
                        if A[j,i]!=0
                            append!(connections_pos,j)
                        end
                    end

                    if connections_pos==[]
                            writedlm(io, [gene_name])
                    else
                        connections = "( "
                        for g in 1:length(connections_pos)

                            if g-1>0
                                connections = string(connections,"or ")
                            end
                            connections = string(connections,"g",string(connections_pos[g])," ")
                        end
                        connections = string(connections, ")")
                        writedlm(io, [gene_name connections])
                    end
                end
            end
    end
end

function BoolODE_textfiles_directed(directory, p)
    txts=readdir(directory)
    for txt in txts
        # set network name
        network_name_BO = string("BoolODE_",txt)

        # read network
        A_path=string(directory,txt)
        A=readdlm(A_path)
        N=size(A)[1]


        # write txt file for BoolODE
        open(network_name_BO, "w") do io
            g = "Gene"
            r = "Rule"
            writedlm(io, [g r])
            for i in 1:N
                gene_name = string("g",string(i))

                connections_pos=[]
                connections_neg=[]
                for j in 1:N
                        if A[j,i]!=0
                            pi = rand(1)[1]
                            if pi<p
                                append!(connections_pos,j)
                            else
                                append!(connections_neg,j)
                            end
                        end
                end

                    if connections_pos==[] && connections_neg==[]
                                writedlm(io, [gene_name])
                    elseif connections_neg==[] && connections_pos!=[]
                            connections = "( "
                            for g in 1:length(connections_pos)

                                if g-1>0
                                    connections = string(connections,"or ")
                                end
                                connections = string(connections,"g",string(connections_pos[g])," ")
                            end
                            connections = string(connections, ")")
                            writedlm(io, [gene_name connections])
                    elseif connections_pos==[] && connections_neg!=[]
                            connections = "not ( "
                            for g in 1:length(connections_neg)

                                if g-1>0
                                    connections = string(connections,"or ")
                                end
                                connections = string(connections,"g",string(connections_neg[g])," ")
                            end
                            connections = string(connections, ")")
                            writedlm(io, [gene_name connections])
                    else
                        connections = "( ( "
                        for g in 1:length(connections_pos)

                            if g-1>0
                                connections = string(connections,"or ")
                            end
                            connections = string(connections,"g",string(connections_pos[g])," ")
                        end
                        connections = string(connections, ") and not ( ")
                        for g in 1:length(connections_neg)

                            if g-1>0
                                connections = string(connections,"or ")
                            end
                            connections = string(connections,"g",string(connections_neg[g])," ")
                        end
                        connections = string(connections, ") )")
                        writedlm(io, [gene_name connections])
                    end
            end
        end
    end
end

function BoolODE_textfiles_directed_random(directory)
    txts=readdir(directory)

    pas = []
    pas_direction = []

    for txt in txts
        # set network name
        network_name_BO = string("BoolODE_",txt)

        # read network
        A_path=string(directory,txt)
        A=readdlm(A_path)
        N=size(A)[1]

        # Set the percentage of activating/inhibiting relations
        p = rand(1)[1]

        # Set the percentage of autoregulations (can also be activating/inhibiting)
        pa = rand(1)[1]
        pa_direction = rand(1)[1]
        append!(pas,pa)
        append!(pas_direction,pa_direction)

        # write txt file for BoolODE
        open(network_name_BO, "w") do io
            g = "Gene"
            r = "Rule"
            writedlm(io, [g r])
            for i in 1:N
                # current gene to look at
                gene_name = string("g",string(i))

                # storing act. and inhib. relations
                connections_pos=[]
                connections_neg=[]

                for j in 1:N
                    # Read act. and inh. relations
                    if A[j,i]!=0
                        pi = rand(1)[1]
                        if pi<p
                            append!(connections_pos,j)
                        else
                            append!(connections_neg,j)
                        end
                    end
                end

                # sample autoregulation prob. for each gene
                ipa = rand(1)[1]
                ipa_direction = rand(1)[1]

                # add autoregulation
                if ipa<pa
                    if ipa_direction<pa_direction
                        append!(connections_pos,i)
                    else
                        append!(connections_neg,i)
                    end
                end

                # Write edges into the file
                if connections_pos==[] && connections_neg==[]
                            writedlm(io, [gene_name])
                elseif connections_neg==[] && connections_pos!=[]
                        connections = "( "
                        for g in 1:length(connections_pos)
                            if g-1>0
                                connections = string(connections,"or ")
                            end
                            connections = string(connections,"g",string(connections_pos[g])," ")
                        end
                        connections = string(connections, ")")
                        writedlm(io, [gene_name connections])
                elseif connections_pos==[] && connections_neg!=[]
                        connections = "not ( "
                        for g in 1:length(connections_neg)
                            if g-1>0
                                connections = string(connections,"or ")
                            end
                            connections = string(connections,"g",string(connections_neg[g])," ")
                        end
                        connections = string(connections, ")")
                        writedlm(io, [gene_name connections])
                else
                    connections = "( ( "
                    for g in 1:length(connections_pos)
                        if g-1>0
                            connections = string(connections,"or ")
                        end
                        connections = string(connections,"g",string(connections_pos[g])," ")
                    end
                    connections = string(connections, ") and not ( ")
                    for g in 1:length(connections_neg)
                        if g-1>0
                            connections = string(connections,"or ")
                        end
                        connections = string(connections,"g",string(connections_neg[g])," ")
                    end
                    connections = string(connections, ") )")
                    writedlm(io, [gene_name connections])
                end
            end
        end
    end
    writedlm("Autoregulation_probs.csv",pas)
    writedlm("Autoregulation_direction_probs.csv",pas_direction)
end

BoolODE_textfiles_directed("C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Julia/samples/temp/",.5)
BoolODE_textfiles("C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Julia/samples/temp/")
BoolODE_textfiles_directed_random("C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Julia/samples/Samples_1/small_world/adjacency_matrices/")
BoolODE_textfiles_directed_random("C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Julia/samples/temp/")
