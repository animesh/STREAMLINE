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

    refCC, refGEFF, refLEFF, refNC, refSP, refASS, refDS, refCT = getTopProp(refGraph, refGraph)


    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' + datasetDict['name']
    dataDict = {}
    dataDict['CC'] = {}
    dataDict['GEFF'] = {}
    dataDict['LEFF'] = {}
    dataDict['ASS'] = {}
    dataDict['DS'] = {}
    dataDict['NC'] = {}
    dataDict['SP'] = {}
    dataDict['CT'] = {}
    dataDict['CC']["ground truth"] = refCC
    dataDict['GEFF']["ground truth"] = refGEFF
    dataDict['LEFF']["ground truth"] = refLEFF
    dataDict['NC']["ground truth"] = refNC
    dataDict['SP']["ground truth"] = refSP
    dataDict['ASS']["ground truth"] = refASS
    dataDict['DS']["ground truth"] = refDS
    dataDict['CT']["ground truth"] = refCT

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
                    edgeWeightTopk = predDF.iloc[2*maxk-1+25].EdgeWeight
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
                dataDict['CC'][algo[0]], dataDict['GEFF'][algo[0]], dataDict['LEFF'][algo[0]], dataDict['NC'][algo[0]], dataDict['SP'][algo[0]], dataDict['ASS'][algo[0]], dataDict['DS'][algo[0]], dataDict['CT'][algo[0]] = getTopProp(predGraph,refGraph)

            else:
                # no edges are predicted, set to 0!
                dataDict['CC'][algo[0]] = 0
                dataDict['GEFF'][algo[0]] = 0
                dataDict['LEFF'][algo[0]] = 0
                dataDict['NC'][algo[0]] = 0
                dataDict['SP'][algo[0]] = 0
                dataDict['ASS'][algo[0]] = 0
                dataDict['DS'][algo[0]] = 0
                dataDict['CT'][algo[0]] = 0
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')

            dataDict['CC'][algo[0]] = 0
            dataDict['GEFF'][algo[0]] = 0
            dataDict['LEFF'][algo[0]] = 0
            dataDict['NC'][algo[0]] = 0
            dataDict['SP'][algo[0]] = 0
            dataDict['ASS'][algo[0]] = 0
            dataDict['DS'][algo[0]] = 0
            dataDict['CT'][algo[0]] = 0

    dataDF = pd.DataFrame(dataDict)
    print("success")
    return dataDF['CC'], dataDF['GEFF'], dataDF['LEFF'], dataDF['NC'], dataDF['SP'], dataDF['ASS'], dataDF['DS'], dataDF['CT']


def getTopProp(inGraph, refGraph):
    #
    # A helper function to compute
    # topological network properties
    #
    #
    # :param inGraph: An graph object of class :class:`networkx.DiGraph`.
    # :type inGraph: :obj:networkx.DiGraph
    #
    # :returns:
    #
    #
    #
    #cc = nx.average_clustering(inGraph)
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


    return cc, g_eff, l_eff, nc, sp, ass, ds, centr



def kl_divergence(d1, d2):
    max_l=max(d1.shape[0],d2.shape[0])

    p = np.zeros(max_l)
    q = np.zeros(max_l)

    p[:d1.shape[0]] = d1
    q[:d2.shape[0]] = d2

    return np.sum(np.where(np.bitwise_and(p != 0,q != 0), p * np.log2(p / q), 0))









'''
                # dataDict['Conn. Comp'][algo[0]], dataDict['FBL'][algo[0]], dataDict['FFL'][algo[0]], dataDict['Mutual'][algo[0]] = getNetProp(predGraph)
                dataDict['FBL'][algo[0]], dataDict['FFL'][algo[0]], dataDict['Mutual'][algo[0]] = getNetProp(predGraph)


                dataDict['FBL'][algo[0]] = dataDict['FBL'][algo[0]]/refFB
                dataDict['FFL'][algo[0]] = dataDict['FFL'][algo[0]]/refFF
                dataDict['Mutual'][algo[0]] = dataDict['Mutual'][algo[0]]/refMI

            else:
                # no edges are predicted, set to 0!
                dataDict['FBL'][algo[0]] = 0
                dataDict['FFL'][algo[0]] = 0
                dataDict['Mutual'][algo[0]] = 0
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')

            dataDict['FBL'][algo[0]] = 0
            dataDict['FFL'][algo[0]] = 0
            dataDict['Mutual'][algo[0]] = 0

    dataDF = pd.DataFrame(dataDict)

    return dataDF['FBL'], dataDF['FFL'], dataDF['Mutual']



def getNetProp(inGraph):
    #
    # A helper function to compute
    # counts of various network motifs.
    #
    #
    # :param inGraph: An graph object of class :class:`networkx.DiGraph`.
    # :type inGraph: :obj:networkx.DiGraph
    #
    # :returns:
    #     - A value corresponding to the number of three-node feedback loops
    #     - A value corresponding to the number of three-node feedforward loops
    #     - A value corresponding to the number of two-node mutual interaction
    #
    #
    #

    # number of weakly connected components in
    # reference network
    # numCC = len(list(nx.weakly_connected_components(inGraph)))

    # number of feedback loop
    # in reference network
    allCyc = nx.simple_cycles(inGraph)
    cycSet = set()
    for cyc in allCyc:
        if len(cyc) == 3:
            cycSet.add(frozenset(cyc))

    numFB = len(cycSet)

    # number of feedfwd loops
    # in reference network
    allPaths = []
    allPathsSet = set()
    for u,v in inGraph.edges():
        allPaths = nx.all_simple_paths(inGraph, u, v, cutoff=2)
        for p in allPaths:
            if len(p) > 2:
                allPathsSet.add(frozenset(p))

    numFF= len(allPathsSet)


    # number of mutual interactions
    numMI = 0.0
    for u,v in inGraph.edges():
        if (v,u) in inGraph.edges():
            numMI += 0.5

    return numFB, numFF, numMI
'''
