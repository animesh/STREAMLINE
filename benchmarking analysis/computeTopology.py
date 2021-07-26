import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from itertools import product, permutations, combinations, combinations_with_replacement
from tqdm import tqdm
import networkx as nx

def Topologies(datasetDict, inputSettings):
    '''
    Computes topological network properties

    :param datasetDict:   A dictionary containing the dataset name, path to reference network.
    :type datasetDict: dict

    :param inputSettings: An object of class :class:`BLEval.InputSettings`.
    :type inputSettings: :class:`BLEval.InputSettings`

    :returns:

    '''

    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ datasetDict['name'] +
                                '/' +datasetDict['trueEdges'],
                                sep = ',',
                                header = 0, index_col = None)

    possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                 r = 2))
    EdgeDict = {'|'.join(p):0 for p in possibleEdges}

    refGraph = nx.DiGraph()

    for key in EdgeDict.keys():
        u = key.split('|')[0]
        v = key.split('|')[1]
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == u) &
               (trueEdgesDF['Gene2'] == v)])>0:
                refGraph.add_edge(u,v)

    numEdges = len(refGraph.edges())
    numNodes = refGraph.number_of_nodes()

    refCCd, refGEFFd, refLEFFd, refNCd, refSPd, refASSd, refiDSd, refoDSd, refiCTd, refoCTd, refJd, refDJd = getTopProp_directed(refGraph, refGraph)
    refCCud, refGEFFud, refLEFFud, refNCud, refSPud, refASSud, refiDSud, refoDSud, refiCTud, refoCTud, refJud, refDJud = getTopProp_undirected(refGraph, refGraph)

    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' + datasetDict['name']
    dataDict = {}

    dataDict['CC'] = {}
    dataDict['GEFF'] = {}
    dataDict['LEFF'] = {}
    dataDict['ASS'] = {}
    dataDict['iDS'] = {}
    dataDict['oDS'] = {}
    dataDict['NC'] = {}
    dataDict['SP'] = {}
    dataDict['iCT'] = {}
    dataDict['oCT'] = {}
    dataDict['HJ'] = {}
    dataDict['DJ'] = {}

    dataDict['CC']["ground truth"] = refCCd
    dataDict['GEFF']["ground truth"] = refGEFFd
    dataDict['LEFF']["ground truth"] = refLEFFd
    dataDict['NC']["ground truth"] = refNCd
    dataDict['SP']["ground truth"] = refSPd
    dataDict['ASS']["ground truth"] = refASSd
    dataDict['iDS']["ground truth"] = refiDSd
    dataDict['oDS']["ground truth"] = refoDSd
    dataDict['iCT']["ground truth"] = refiCTd
    dataDict['oCT']["ground truth"] = refoCTd
    dataDict['HJ']["ground truth"] = refJd
    dataDict['DJ']["ground truth"] = refDJd

    dataDict['CC']["ground truth undirected"] = refCCud
    dataDict['GEFF']["ground truth undirected"] = refGEFFud
    dataDict['LEFF']["ground truth undirected"] = refLEFFud
    dataDict['NC']["ground truth undirected"] = refNCud
    dataDict['SP']["ground truth undirected"] = refSPud
    dataDict['ASS']["ground truth undirected"] = refASSud
    dataDict['iDS']["ground truth undirected"] = refiDSud
    dataDict['oDS']["ground truth undirected"] = refoDSud
    dataDict['iCT']["ground truth undirected"] = refiCTud
    dataDict['oCT']["ground truth undirected"] = refoCTud
    dataDict['HJ']["ground truth undirected"] = refJud
    dataDict['DJ']["ground truth undirected"] = refDJud

    for algo in tqdm(inputSettings.algorithms,
                     total = len(inputSettings.algorithms), unit = " Algorithms"):
        # if algo[0] == 'PPCOR' or algo[0] == 'PIDC':
        #     continue
        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
             # Initialize Precsion

            predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)


            predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
            predDF.drop_duplicates(keep = 'first', inplace=True)
            predDF.reset_index(drop = True,  inplace= True)
            # check if ranked edges list is empty
            # if so, it is just set to an empty set

            if not predDF.shape[0] == 0:

                # we want to ensure that we do not include
                # edges without any edge weight
                # so check if the non-zero minimum is
                # greater than the edge weight of the top-kth
                # node, else use the non-zero minimum value.
                predDF.EdgeWeight = predDF.EdgeWeight.round(6)
                predDF.EdgeWeight = predDF.EdgeWeight.abs()

                # Use num True edges or the number of
                # edges in the dataframe, which ever is lower
                maxk = min(predDF.shape[0], numEdges)

                if algo[0] == 'PPCOR':
                    edgeWeightTopk = predDF.iloc[2*maxk-1+numNodes].EdgeWeight
                elif algo[0] == 'PIDC':
                    edgeWeightTopk = predDF.iloc[2*maxk-1].EdgeWeight
                else:
                    edgeWeightTopk = predDF.iloc[maxk-1].EdgeWeight

                nonZeroMin = np.nanmin(predDF.EdgeWeight.replace(0, np.nan).values)
                bestVal = max(nonZeroMin, edgeWeightTopk)

                if algo[0] == 'PPCOR':
                    newDF = predDF.loc[(predDF['EdgeWeight'] >= bestVal) & (predDF['EdgeWeight'] < 1)]
                else:
                    newDF = predDF.loc[(predDF['EdgeWeight'] >= bestVal)]


                predGraph = nx.DiGraph()


                for key in EdgeDict.keys():
                    u = key.split('|')[0]
                    v = key.split('|')[1]
                    if len(newDF.loc[(newDF['Gene1'] == u) &
                           (newDF['Gene2'] == v)])>0:
                            predGraph.add_edge(u,v)

                # dataDict['Conn. Comp'][algo[0]], dataDict['FBL'][algo[0]], dataDict['FFL'][algo[0]], dataDict['Mutual'][algo[0]] = getNetProp(predGraph)

                # Calculation for random predictors

                # rcc, rg_eff, rl_eff, rnc, rsp, rass, rds, rds, rcentr, rcentr, rh_jac, rd_jac = 0,0,0,0,0,0,0,0,0,0,0,0
                # rdcc, rdg_eff, rdl_eff, rdnc, rdsp, rdass, rdds, rdds, rdcentr, rdcentr, rdh_jac, rdd_jac = 0,0,0,0,0,0,0,0,0,0,0,0
                r_vals = np.zeros(12)
                rd_vals = np.zeros(12)

                for j in range(0,29):
                    rGr = nx.gnm_random_graph(numNodes, numEdges, directed=False)
                    rGrd = nx.gnm_random_graph(numNodes, numEdges, directed=True)
                    r_vals = r_vals+1/30*np.array(getTopProp_undirected(rGr,refGraph))
                    rd_vals = rd_vals+1/30*np.array(getTopProp_directed(rGrd,refGraph))
                    # rcc, rg_eff, rl_eff, rnc, rsp, rass, rds, rds, rcentr, rcentr, rh_jac, rd_jac = rcc, rg_eff, rl_eff, rnc, rsp, rass, rds, rds, rcentr, rcentr, rh_jac, rd_jac + getTopProp_undirected(rGr,refGraph)*1/10
                    # rdcc, rdg_eff, rdl_eff, rdnc, rdsp, rdass, rdds, rdds, rdcentr, rdcentr, rdh_jac, rdd_jac = rdcc, rdg_eff, rdl_eff, rdnc, rdsp, rdass, rdds, rdds, rdcentr, rdcentr, rdh_jac, rdd_jac + getTopProp_directed(rGrd,refGraph)*1/10


                dataDict['CC']["random undirected"], dataDict['GEFF']["random undirected"], dataDict['LEFF']["random undirected"], dataDict['NC']["random undirected"], dataDict['SP']["random undirected"], dataDict['ASS']["random undirected"] = r_vals[0], r_vals[1], r_vals[2], r_vals[3], r_vals[4], r_vals[5]
                dataDict['iDS']["random undirected"], dataDict['oDS']["random undirected"], dataDict['iCT']["random undirected"], dataDict['oCT']["random undirected"], dataDict['HJ']["random undirected"], dataDict['DJ']["random undirected"] = r_vals[6], r_vals[7], r_vals[8], r_vals[9], r_vals[10], r_vals[11]
                dataDict['CC']["random directed"], dataDict['GEFF']["random directed"], dataDict['LEFF']["random directed"], dataDict['NC']["random directed"], dataDict['SP']["random directed"], dataDict['ASS']["random directed"] = rd_vals[0], rd_vals[1], rd_vals[2], rd_vals[3], rd_vals[4], rd_vals[5]
                dataDict['iDS']["random directed"], dataDict['oDS']["random directed"], dataDict['iCT']["random directed"], dataDict['oCT']["random directed"], dataDict['HJ']["random directed"], dataDict['DJ']["random directed"] = rd_vals[6], rd_vals[7], rd_vals[8], rd_vals[9], rd_vals[10], rd_vals[11]

                # Calculation for inferred networks
                if algo[0] == 'PPCOR' or algo[0] == 'PIDC':
                    dataDict['CC'][algo[0]], dataDict['GEFF'][algo[0]], dataDict['LEFF'][algo[0]], dataDict['NC'][algo[0]], dataDict['SP'][algo[0]], dataDict['ASS'][algo[0]], dataDict['iDS'][algo[0]], dataDict['oDS'][algo[0]], dataDict['iCT'][algo[0]], dataDict['oCT'][algo[0]], dataDict['HJ'][algo[0]], dataDict['DJ'][algo[0]] = getTopProp_undirected(predGraph,refGraph)
                else:
                    dataDict['CC'][algo[0]], dataDict['GEFF'][algo[0]], dataDict['LEFF'][algo[0]], dataDict['NC'][algo[0]], dataDict['SP'][algo[0]], dataDict['ASS'][algo[0]], dataDict['iDS'][algo[0]], dataDict['oDS'][algo[0]], dataDict['iCT'][algo[0]], dataDict['oCT'][algo[0]], dataDict['HJ'][algo[0]], dataDict['DJ'][algo[0]] = getTopProp_directed(predGraph,refGraph)

            else:
                # no edges are predicted, set to 0!
                dataDict['CC'][algo[0]] = 0
                dataDict['GEFF'][algo[0]] = 0
                dataDict['LEFF'][algo[0]] = 0
                dataDict['NC'][algo[0]] = 0
                dataDict['SP'][algo[0]] = 0
                dataDict['ASS'][algo[0]] = 0
                dataDict['iDS'][algo[0]] = 0
                dataDict['oDS'][algo[0]] = 0
                dataDict['iCT'][algo[0]] = 0
                dataDict['oCT'][algo[0]] = 0
                dataDict['HJ'][algo[0]] = 0
                dataDict['DJ'][algo[0]] = 0
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')

            dataDict['CC'][algo[0]] = 0
            dataDict['GEFF'][algo[0]] = 0
            dataDict['LEFF'][algo[0]] = 0
            dataDict['NC'][algo[0]] = 0
            dataDict['SP'][algo[0]] = 0
            dataDict['ASS'][algo[0]] = 0
            dataDict['iDS'][algo[0]] = 0
            dataDict['oDS'][algo[0]] = 0
            dataDict['iCT'][algo[0]] = 0
            dataDict['oCT'][algo[0]] = 0
            dataDict['HJ'][algo[0]] = 0
            dataDict['DJ'][algo[0]] = 0

    dataDF = pd.DataFrame(dataDict)
    print("success")
    return dataDF['CC'], dataDF['GEFF'], dataDF['LEFF'], dataDF['NC'], dataDF['SP'], dataDF['ASS'], dataDF['iDS'], dataDF['oDS'], dataDF['iCT'], dataDF['oCT'], dataDF['HJ'], dataDF['DJ']


def getTopProp_undirected(inGraph, refGraph):
    # A helper function to compute
    # topological network properties
    inGraphU = nx.DiGraph.to_undirected(inGraph)
    refGraphU = nx.DiGraph.to_undirected(refGraph)
    cc = nx.average_clustering(inGraphU)
    g_eff = nx.global_efficiency(inGraphU)
    l_eff = nx.local_efficiency(inGraphU)
    nc = nx.number_connected_components(inGraphU)
    #ass = nx.degree_assortativity_coefficient(inGraph)
    ass = nx.degree_assortativity_coefficient(inGraphU)

    if nc==1:
        # sp=nx.average_shortest_path_length(inGraph)
        sp=nx.average_shortest_path_length(inGraphU)
    else:
        sp=0

    degree_hist = nx.degree_histogram(inGraphU)
    degree_hist = np.array(degree_hist, dtype=float)
    degree_prob = degree_hist/inGraphU.order()

    N=inGraphU.order()
    max_d = max(degree_hist)
    centr = float((N*max_d - sum(degree_hist)))/(N-1)**2

    degree_hist_ref = nx.degree_histogram(refGraphU)
    degree_hist_ref = np.array(degree_hist_ref, dtype=float)
    degree_prob_ref = degree_hist_ref/refGraphU.order()

    ds = kl_divergence(degree_prob,degree_prob_ref)

    ranks_1 = nx.pagerank(inGraphU)
    ranks_1_s = sorted(ranks_1.items(), key=lambda kv: kv[1],reverse=True)
    ranks_2 = nx.pagerank(refGraphU)
    ranks_2_s = sorted(ranks_2.items(), key=lambda kv: kv[1],reverse=True)


    hubs1 = [a_tuple1[0] for a_tuple1 in ranks_1_s[0:9]]
    hubs2 = [a_tuple2[0] for a_tuple2 in ranks_2_s[0:9]]
    h_jac = jaccard(hubs1, hubs2)

    deg_centr1 = nx.degree_centrality(inGraphU)
    deg_centr1_s = sorted(deg_centr1.items(), key=lambda kv: kv[1],reverse=True)
    deg_centr2 = nx.degree_centrality(refGraphU)
    deg_centr2_s = sorted(deg_centr2.items(), key=lambda kv: kv[1],reverse=True)

    dhubs1 = [d_tuple1[0] for d_tuple1 in list(deg_centr1_s)[0:4]]
    dhubs2 = [d_tuple2[0] for d_tuple2 in list(deg_centr2_s)[0:4]]

    d_jac = jaccard(dhubs1, dhubs2)


    return cc, g_eff, l_eff, nc, sp, ass, ds, ds, centr, centr, h_jac, d_jac

def getTopProp_directed(inGraph, refGraph):
    # A helper function to compute
    # topological network properties
    inGraphU = nx.DiGraph.to_undirected(inGraph)
    refGraphU = nx.DiGraph.to_undirected(refGraph)
    cc = nx.average_clustering(inGraph)
    g_eff = nx.global_efficiency(inGraphU)
    l_eff = nx.local_efficiency(inGraphU)
    nc = nx.number_connected_components(inGraphU)
    #ass = nx.degree_assortativity_coefficient(inGraph)
    ass = nx.degree_assortativity_coefficient(inGraph)

    if nc==1:
        # sp=nx.average_shortest_path_length(inGraph)
        sp=nx.average_shortest_path_length(inGraph)
    else:
        sp=0

    N_in=inGraph.order()
    indegrees = [x[1] for x in list(inGraph.in_degree)]
    max_in = max(indegrees)
    centr_in = float((N_in*max_in - sum(indegrees)))/(N_in-1)**2

    N_out=inGraph.order()
    outdegrees = [x[1] for x in list(inGraph.out_degree)]
    max_out = max(outdegrees)
    centr_out = float((N_out*max_out - sum(outdegrees)))/(N_out-1)**2

    N_in_ref = refGraph.order()
    indegrees_ref = [x[1] for x in list(refGraph.in_degree)]
    max_in_ref = max(indegrees_ref)

    N_out_ref = refGraph.order()
    outdegrees_ref = [x[1] for x in list(refGraph.out_degree)]
    max_out_ref = max(outdegrees_ref)

    l_in = max(max_in,max_in_ref)
    indegs1=np.zeros(l_in)
    indegs2=np.zeros(l_in)
    for i in range(1,l_in):
        indegs1[i]=indegrees.count(i)
        indegs2[i]=indegrees_ref.count(i)

    ds_in = kl_divergence(indegs1,indegs2)

    l_out = max(max_out,max_out_ref)
    outdegs1=np.zeros(l_out)
    outdegs2=np.zeros(l_out)
    for i in range(1,l_out):
        outdegs1[i]=outdegrees.count(i)
        outdegs2[i]=outdegrees_ref.count(i)

    ds_out = kl_divergence(outdegs1,outdegs2)


    ranks_1 = nx.pagerank(inGraph)
    ranks_1_s = sorted(ranks_1.items(), key=lambda kv: kv[1],reverse=True)
    ranks_2 = nx.pagerank(refGraph)
    ranks_2_s = sorted(ranks_2.items(), key=lambda kv: kv[1],reverse=True)


    hubs1 = [a_tuple1[0] for a_tuple1 in ranks_1_s[0:9]]
    hubs2 = [a_tuple2[0] for a_tuple2 in ranks_2_s[0:9]]
    h_jac = jaccard(hubs1, hubs2)

    deg_centr1 = nx.degree_centrality(inGraph)
    deg_centr1_s = sorted(deg_centr1.items(), key=lambda kv: kv[1],reverse=True)
    deg_centr2 = nx.degree_centrality(refGraph)
    deg_centr2_s = sorted(deg_centr2.items(), key=lambda kv: kv[1],reverse=True)

    dhubs1 = [d_tuple1[0] for d_tuple1 in list(deg_centr1_s)[0:4]]
    dhubs2 = [d_tuple2[0] for d_tuple2 in list(deg_centr2_s)[0:4]]

    d_jac = jaccard(dhubs1, dhubs2)


    return cc, g_eff, l_eff, nc, sp, ass, ds_in, ds_out, centr_in, centr_out, h_jac, d_jac



def kl_divergence(d1, d2):
    max_l=max(d1.shape[0],d2.shape[0])

    p = np.zeros(max_l)
    q = np.zeros(max_l)

    p[:d1.shape[0]] = d1
    q[:d2.shape[0]] = d2

    return np.sum(np.where(np.bitwise_and(p != 0,q != 0), p * np.log2(p / q), 0))

def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union
