# MixGGM
Inferring multiple microbial interaction networks from compositional taxa_count matrix using the Mixture of Graphical Guassian Model
"MixGGM" is a package written in R, which receives a sample-taxa count matrix and extracts k(number of components) different interaction networks between taxa.  

## 1-Prerequisites
Following packages must be installed in R environment before using "MixMPLN":
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
After installing the previous prerequisities, please run the following lines to install the "MixMPLN"
```
library("devtools")
install_github("sahatava/MixGGM")
library("MixGGM")
```
## 3-Extracting K interaction networks
Please save your matrix in a cvs file with row names of sample IDs and a header with names of the taxa. After entering your csv file name and other arguments in following command, K different file will be save in the current directory. The type of the outputs can be "precision"(precision matrix) , "pratial"(partial correlation) ,"adj"(adjacency matrix).   
##### Usage
```
out_MixMPLN = MixMPLN(M , K , penalty )
```
##### Arguments
```
file--------> name of the csv file which includes the sample-taxa matrix
K-----------> number of components 
penalty------> method to select tuning parameter: "no_sparsity" , "CV" (cross validation) , "StARS" , "fixed" , "iterative" 
out---------> the type of output for each componenet: "precision"(precision matrix) , "pratial"(partial correlation) ,"adj"(adjacency matrix) 
threshold---> a value between 0 and 1. This threshold is used to generate adjacency matrix from partial correlation
```

##### Values
```
out_MixMPLN$precision[[i]]-----> i_th precision matrix(dXd)
out_MixMPLN$partial[[i]]-------> i_th partial correlation matrix(dXd)
out_MixMPLN$lambda[[i]]--------> i_th lambda matrix(NXd)
out_MixMPLN$pi-----------------> clustering coefficiant 
out_MixMPLN$cluster------------> clustering membership
```

## 5-Auxiliary functions

### silhouette method to find the optimal K
```
fviz_nbclust(M , MixMPLN ,penalty ,k.max, method = c("silhouette"))
```
##### Arguments
```
M-----------> sample-taxa count matrix
K.max-------> k=1 to k=k.max will be checked
penalty------> method to select tuning parameter: "no_sparsity" , "CV" (cross validation) , "StARS" , "fixed" , "iterative" 
```
##### Values
```
Running the mentioned commad will plot a diagram which shows optimum K 
```
 
### Graphical display
```
visualize(partial , threshold , taxa_IDs )
```
##### Arguments
```
partial--------> partial correlation matrix corresponding to one component
threshold------> Threshold to get the adjacency matrix of the graph
taxa_IDs-------> taxa ID names as a list of characters
```
##### Values
```
The graphical plot
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

### Comparing the estimated precision matrix and the real one for synthetic data
```
out = compare(K , input1 , input2)
```
##### Arguments
```
input1 -------> real partial correlation
input2 -------> estimated partial correlation
```
##### Values
```
out$frob------> forbenious norm of the difference
out$fr--------> relative difference
out$rms-------> rms difference
out$sn--------> sensitivity
out$sp--------> specificity
out$ROC-------> area under the ROC 
```


## 6-Example
### Example on synthetic data
```
out_generate = generate(K=2 , N=30 , d=30 , sp=0.8, type="orig")
out_MixMPLN = MixMPLN( out_generate$M , K=2 , penalty="CV")
out = compare(K=2 , out_generate$real_precision  , out_MixMPLN$precision)
``` 
### Example on real data
```
M=as.matrix(read.csv("real_data.csv",sep=",",header=TRUE,row.names=1,check.names=FALSE))
fviz_nbclust(M , MixMPLN ,"CV",k.max, method = c("silhouette"))
out_MixMPLN = MixMPLN( M , K=2 , penalty="CV")
visualize(out_MixMPLN$partial[[1]] , threshold=0.5 , taxa_IDs=colnames(M) )
```
