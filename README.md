# MixGGM
Inferring multiple microbial interaction networks from compositional taxa_count matrix using the Mixture of Graphical Guassian Model
"MixGGM" is a package written in R, which receives a sample-taxa count matrix and extracts k(number of components) different interaction networks between taxa.  

## 1-Prerequisites
Following packages must be installed and loaded in R environment before using "MixGGM":
```
library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(cvTools)
library(huge)
library(edgeR)
library(factoextra)
library(igraph)
library(ROCR)
```
## 2-Installation
After installing the previous prerequisities, please run the following lines to install the "MixGGM"
```
library("devtools")
install_github("sahatava/MixGGM")
library("MixGGM")
```
 
## 3-Extracting K interaction networks
Please save your data as matrix which its rows represents samples and its columns represents the taxa names. If you are interested in visualising the final graphs, you need to save the taxa names as colnames of this input matrix. 
K( number of components ) is another input for MixGGM which should be specified by the user.   
For initialize the optimization using the Kmeans clustering please specified init="Kmeans" and rep=1. Otherwise, init="Random" and rep can be any integer number. As higher value for rep results in higher running time, we suggest to set it as 3 to trade off between the time and accuracy.  
##### Usage
```
out_MixGGM = MixGGM(M , K , penalty , init , rep )
```
##### Arguments
```
M--------> name of the csv file which includes the sample-taxa matrix
K-----------> number of components 
penalty------> method to select tuning parameter: "no_sparsity" , "CV" (cross validation) , "StARS" , "fixed" , "iterative" 
init---------> initialization can be "Kmeans" or "Random"
rep----------> number of repeats when initialization is "Random"
out---------> the type of output for each componenet: precision matrix, partial correlation, lambda, mixing coefficiant, clustering membership 
threshold---> a value between 0 and 1. This threshold is used to generate adjacency matrix from partial correlation
```

##### Values
```
out_MixGGM$precision[[i]]-----> i_th precision matrix(dXd)
out_MixGGM$partial[[i]]-------> i_th partial correlation matrix(dXd)
out_MixGGM$lambda[[i]]--------> i_th lambda matrix(NXd)
out_MixGGM$pi-----------------> mixing coefficiant
out_MixGGM$cluster------------> clustering membership
```

## 4-Auxiliary functions
 
### Graphical display
```
visualize(partial , threshold )
```
##### Arguments
```
partial--------> partial correlation matrix corresponding to one component
threshold------> Threshold to get the adjacency matrix of the graph
```
##### Values
```
The graphical plot for the taxa interaction network. 
```

### Converting the partial corelation to the adjacency matrix
```
out_adj = partial2adj(partial, threshold)
```
##### Arguments
```
partial-------------> partial correlation matrix as a result of running the MixMPLN
threshold-----------> a value between 0 and 1
```
##### Values
```
out_adj-------------> 1 and -1 indicate edge existance, 0 indicates no edge 
```

 
## 5-Example on real data
 ` 
```
M=as.matrix(read.csv("real_data.csv",sep=",",header=TRUE,row.names=1,check.names=FALSE))
out_MixGGM = MixMPLN( M , K=2 , penalty="CV", init = "Random" , rep = 3)
visualize(out_MixGGM$partial[[1]] , threshold=0.3 )
```

## Publication
 
