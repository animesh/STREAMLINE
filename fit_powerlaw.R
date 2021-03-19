require(igraph)
library(doParallel)# for parallel processing!
library(foreach)
require(poweRlaw)
require(tseries)

powerlawfit<-function(degree) {
  poweRlawFit = displ$new(degree)
  estPL = estimate_xmin(poweRlawFit) #estimate parameters! Recommended than "estimate_pars" which uses the smallest Xmin
  poweRlawFit$setXmin(estPL) #set xmin
  poweRlawFit$setPars(estimate_pars(poweRlawFit)) #Set parameters
  return(poweRlawFit)
}


### Example power law fit

A_mESC <- read.matrix("C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Julia/mESC_ChIP_adjacency.txt")
mESC_network=graph_from_adjacency_matrix(A_mESC,"directed")
gridPowerLawIgraph(mESC_network,"out")
degree<-degree(mESC_network, mode = 'in') #get degree#
degree<-degree[degree>0]


#Execute function for fitting a power law using poweRlaw:
poweRlawFit<-powerlawfit(degree)
poweRlawFit$pars