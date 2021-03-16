using LightGraphs
using Plots
using GraphPlot
using Random
using Distributions

function BEELINE_textfiles(directory)
    txts=readdir(directory)
    for txt in txts
            # set network name
            network_name_BE = string("BEELINE_",txt)

            # read network
            A_path=string("samples/",txt)
            A=readdlm(A_path)
            N=size(A)[1]

            #matrix storing the positive connections
            B=zeros(N,N)
            p1=0
            p2=1

            # write txt file for BEELINE
            open(network_name_BE, "w") do jo
                print(jo, "Gene1,Gene2,Type")
                for i in 1:N
                    gene_name_2 = string("g",string(i))
                    for j in 1:N
                        if A[j,i]!=0
                            print(jo, "\n")
                            gene_name_1 = string("g",string(j))
                            print(jo,gene_name_1,",",gene_name_2,",","+")
                            end
                        end
                    end
                end
            end
    end
end

BEELINE_textfiles("samples")
