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

def Hubs(datasetDict, inputSettings):
    '''
    Computes hub similarity

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

    refPR, refIC, refOC, refEIG, refBT, refRD = Hubs_directed(refGraph, refGraph)

    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' + datasetDict['name']
    dataDict = {}

    dataDict['PR'] = {}
    dataDict['IC'] = {}
    dataDict['OC'] = {}
    dataDict['EIG'] = {}
    dataDict['BT'] = {}
    dataDict['RD'] = {}

    dataDict['PR']["ground truth"] = refPR
    dataDict['IC']["ground truth"] = refIC
    dataDict['OC']["ground truth"] = refOC
    dataDict['EIG']["ground truth"] = refEIG
    dataDict['BT']["ground truth"] = refBT
    dataDict['RD']["ground truth"] = refRD

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

                r_vals = np.zeros(6)
                rd_vals = np.zeros(6)
                mapping = {0: "g1", 1: "g2", 2: "g3", 3: "g4", 4: "g5", 5: "g6", 6: "g7", 7: "g8", 8: "g9", 9: "g10", 10: "g1", 11: "g12", 12: "g13", 13: "g4", 14: "g15"}

		# This is the statistical approach to calculating the expected jaccard coefficient of a random predictor
		# Here we use the output statistic with N=300 samples 
		# For the results in the paper we used an exact formula for the expected value, this approach here is to check whether the formula holds true
                for j in range(0,299):
                    rGr = nx.gnm_random_graph(numNodes, numEdges, directed=False)
                    rGrd = nx.gnm_random_graph(numNodes, numEdges, directed=True)
                    rGr = nx.relabel_nodes(rGr, mapping)
                    rGrd = nx.relabel_nodes(rGrd, mapping)
                    r_vals = r_vals+1/300*np.array(Hubs_undirected(rGr,refGraph))
                    rd_vals = rd_vals+1/300*np.array(Hubs_directed(rGrd,refGraph))
                    # rcc, rg_eff, rl_eff, rnc, rsp, rass, rds, rds, rcentr, rcentr, rh_jac, rd_jac = rcc, rg_eff, rl_eff, rnc, rsp, rass, rds, rds, rcentr, rcentr, rh_jac, rd_jac + getTopProp_undirected(rGr,refGraph)*1/10
                    # rdcc, rdg_eff, rdl_eff, rdnc, rdsp, rdass, rdds, rdds, rdcentr, rdcentr, rdh_jac, rdd_jac = rdcc, rdg_eff, rdl_eff, rdnc, rdsp, rdass, rdds, rdds, rdcentr, rdcentr, rdh_jac, rdd_jac + getTopProp_directed(rGrd,refGraph)*1/10


                dataDict['PR']["random undirected"], dataDict['IC']["random undirected"], dataDict['OC']["random undirected"], dataDict['EIG']["random undirected"], dataDict['BT']["random undirected"], dataDict['RD']["random undirected"] = r_vals[0], r_vals[1], r_vals[2], r_vals[3], r_vals[4], r_vals[5]
                dataDict['PR']["random directed"], dataDict['IC']["random directed"], dataDict['OC']["random directed"], dataDict['EIG']["random directed"], dataDict['BT']["random directed"], dataDict['RD']["random directed"] = rd_vals[0], rd_vals[1], rd_vals[2], rd_vals[3], rd_vals[4], rd_vals[5]


                if algo[0] == 'PPCOR' or algo[0] == 'PIDC':
                    dataDict['PR'][algo[0]], dataDict['IC'][algo[0]], dataDict['OC'][algo[0]], dataDict['EIG'][algo[0]], dataDict['BT'][algo[0]], dataDict['RD'][algo[0]] = Hubs_undirected(predGraph,refGraph)
                else:
                    dataDict['PR'][algo[0]], dataDict['IC'][algo[0]], dataDict['OC'][algo[0]], dataDict['EIG'][algo[0]], dataDict['BT'][algo[0]], dataDict['RD'][algo[0]] = Hubs_directed(predGraph,refGraph)
            else:
                # no edges are predicted, set to 0!
                dataDict['PR'][algo[0]] = 0
                dataDict['IC'][algo[0]] = 0
                dataDict['OC'][algo[0]] = 0
                dataDict['EIG'][algo[0]] = 0
                dataDict['BT'][algo[0]] = 0
                dataDict['RD'][algo[0]] = 0
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')

            dataDict['PR'][algo[0]] = 0
            dataDict['IC'][algo[0]] = 0
            dataDict['OC'][algo[0]] = 0
            dataDict['EIG'][algo[0]] = 0
            dataDict['BT'][algo[0]] = 0
            dataDict['RD'][algo[0]] = 0

    dataDF = pd.DataFrame(dataDict)
    print("success")
    return dataDF['PR'], dataDF['IC'], dataDF['OC'], dataDF['EIG'], dataDF['BT'], dataDF['RD']


def Hubs_undirected(inGraph, refGraph):
    # A helper function to compute
    # hub similarity for undirected output
    inGraphU = nx.DiGraph.to_undirected(inGraph)
    refGraphU = nx.DiGraph.to_undirected(refGraph)
    pranks_1 = nx.pagerank(inGraphU, max_iter=500)
    pranks_1_s = sorted(pranks_1.items(), key=lambda kv: kv[1],reverse=True)
    pranks_2 = nx.pagerank(refGraphU, max_iter=500)
    pranks_2_s = sorted(pranks_2.items(), key=lambda kv: kv[1],reverse=True)
    prhubs1 = [a_tuple1[0] for a_tuple1 in pranks_1_s[0:3]]
    prhubs2 = [a_tuple2[0] for a_tuple2 in pranks_2_s[0:3]]
    pr = jaccard(prhubs1, prhubs2)

    deg_centr1 = nx.degree_centrality(inGraphU)
    deg_centr1_s = sorted(deg_centr1.items(), key=lambda kv: kv[1],reverse=True)
    deg_centr2 = nx.degree_centrality(refGraphU)
    deg_centr2_s = sorted(deg_centr2.items(), key=lambda kv: kv[1],reverse=True)
    dhubs1 = [d_tuple1[0] for d_tuple1 in list(deg_centr1_s)[0:3]]
    dhubs2 = [d_tuple2[0] for d_tuple2 in list(deg_centr2_s)[0:3]]
    centr = jaccard(dhubs1, dhubs2)

    eig_centr1 = nx.eigenvector_centrality(inGraphU, max_iter=500)
    eig_centr1_s = sorted(eig_centr1.items(), key=lambda kv: kv[1],reverse=True)
    eig_centr2 = nx.eigenvector_centrality(refGraphU, max_iter=500)
    eig_centr2_s = sorted(eig_centr2.items(), key=lambda kv: kv[1],reverse=True)
    eighubs1 = [d_tuple1[0] for d_tuple1 in list(eig_centr1_s)[0:3]]
    eighubs2 = [d_tuple2[0] for d_tuple2 in list(eig_centr2_s)[0:3]]
    eig = jaccard(eighubs1, eighubs2)

    bt_centr1 = nx.betweenness_centrality(inGraphU)
    bt_centr1_s = sorted(bt_centr1.items(), key=lambda kv: kv[1],reverse=True)
    bt_centr2 = nx.betweenness_centrality(refGraphU)
    bt_centr2_s = sorted(bt_centr2.items(), key=lambda kv: kv[1],reverse=True)
    bthubs1 = [d_tuple1[0] for d_tuple1 in list(bt_centr1_s)[0:3]]
    bthubs2 = [d_tuple2[0] for d_tuple2 in list(bt_centr2_s)[0:3]]
    bt = jaccard(bthubs1, bthubs2)

    rad_centr1 = radiality(inGraph)
    rad_centr1_s = sorted(rad_centr1.items(), key=lambda kv: kv[1],reverse=True)
    rad_centr2 = radiality(refGraph)
    rad_centr2_s = sorted(rad_centr2.items(), key=lambda kv: kv[1],reverse=True)
    radhubs1 = [d_tuple1[0] for d_tuple1 in list(rad_centr1_s)[0:3]]
    radhubs2 = [d_tuple2[0] for d_tuple2 in list(rad_centr2_s)[0:3]]
    rad = jaccard(radhubs1, radhubs2)
    #rad = 0

    return pr, centr, centr, eig, bt, rad

def Hubs_directed(inGraph, refGraph):
    # A helper function to compute
    # hub similarity for directed output
    pranks_1 = nx.pagerank(inGraph, max_iter=500)
    pranks_1_s = sorted(pranks_1.items(), key=lambda kv: kv[1],reverse=True)
    pranks_2 = nx.pagerank(refGraph, max_iter=500)
    pranks_2_s = sorted(pranks_2.items(), key=lambda kv: kv[1],reverse=True)
    prhubs1 = [a_tuple1[0] for a_tuple1 in pranks_1_s[0:3]]
    prhubs2 = [a_tuple2[0] for a_tuple2 in pranks_2_s[0:3]]
    pr = jaccard(prhubs1, prhubs2)

    ideg_centr1 = nx.in_degree_centrality(inGraph)
    ideg_centr1_s = sorted(ideg_centr1.items(), key=lambda kv: kv[1],reverse=True)
    ideg_centr2 = nx.in_degree_centrality(refGraph)
    ideg_centr2_s = sorted(ideg_centr2.items(), key=lambda kv: kv[1],reverse=True)
    idhubs1 = [d_tuple1[0] for d_tuple1 in list(ideg_centr1_s)[0:3]]
    idhubs2 = [d_tuple2[0] for d_tuple2 in list(ideg_centr2_s)[0:3]]
    icentr = jaccard(idhubs1, idhubs2)

    odeg_centr1 = nx.out_degree_centrality(inGraph)
    odeg_centr1_s = sorted(odeg_centr1.items(), key=lambda kv: kv[1],reverse=True)
    odeg_centr2 = nx.out_degree_centrality(refGraph)
    odeg_centr2_s = sorted(odeg_centr2.items(), key=lambda kv: kv[1],reverse=True)
    odhubs1 = [d_tuple1[0] for d_tuple1 in list(odeg_centr1_s)[0:3]]
    odhubs2 = [d_tuple2[0] for d_tuple2 in list(odeg_centr2_s)[0:3]]
    ocentr = jaccard(odhubs1, odhubs2)

    # eig_centr1 = nx.eigenvector_centrality(inGraph)
    # eig_centr1_s = sorted(eig_centr1.items(), key=lambda kv: kv[1],reverse=True)
    # eig_centr2 = nx.eigenvector_centrality(refGraph)
    # eig_centr2_s = sorted(eig_centr2.items(), key=lambda kv: kv[1],reverse=True)
    # eighubs1 = [d_tuple1[0] for d_tuple1 in list(eig_centr1_s)[0:4]]
    # eighubs2 = [d_tuple2[0] for d_tuple2 in list(eig_centr2_s)[0:4]]
    # eig = jaccard(eighubs1, eighubs2)
    eig = 0

    bt_centr1 = nx.betweenness_centrality(inGraph)
    bt_centr1_s = sorted(bt_centr1.items(), key=lambda kv: kv[1],reverse=True)
    bt_centr2 = nx.betweenness_centrality(refGraph)
    bt_centr2_s = sorted(bt_centr2.items(), key=lambda kv: kv[1],reverse=True)
    bthubs1 = [d_tuple1[0] for d_tuple1 in list(bt_centr1_s)[0:3]]
    bthubs2 = [d_tuple2[0] for d_tuple2 in list(bt_centr2_s)[0:3]]
    bt = jaccard(bthubs1, bthubs2)

    rad_centr1 = radiality(inGraph)
    rad_centr1_s = sorted(rad_centr1.items(), key=lambda kv: kv[1],reverse=True)
    rad_centr2 = radiality(refGraph)
    rad_centr2_s = sorted(rad_centr2.items(), key=lambda kv: kv[1],reverse=True)
    radhubs1 = [d_tuple1[0] for d_tuple1 in list(rad_centr1_s)[0:3]]
    radhubs2 = [d_tuple2[0] for d_tuple2 in list(rad_centr2_s)[0:3]]
    rad = jaccard(radhubs1, radhubs2)
    #rad = 0

    return pr, icentr, ocentr, eig, bt, rad


### Helper functions

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

def radiality(Graph):
	fw = nx.floyd_warshall(Graph)
	results = {a:dict(b) for a,b in fw.items()}
	N = Graph.order()
	DistanceMatrix = np.zeros((N,N))
	for i in range(0,N):
		for j in range(0,N):
			key1 = list(results.keys())[i]
			key2 = list(results.keys())[j]
			# DistanceMatrix[i,j]=results[i][j]
			DistanceMatrix[i,j]=results[key1][key2]
	if nx.is_connected(Graph.to_undirected()):
		RD=nx.diameter(Graph.to_undirected())*np.ones((N,N))+1-DistanceMatrix
	else:
		S = max(nx.connected_component_subgraphs(Graph.to_undirected()), key=len)
		RD=nx.diameter(S)*np.ones((N,N))+1-DistanceMatrix
	radialities=(RD.sum(axis=1)-RD.diagonal())/(len(RD.diagonal())-1)
	rad_dict = dict((k,radialities[k]) for k in range(0,N))
	return rad_dict
