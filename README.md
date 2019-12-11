# ISODATA_clustering
ISODATA Clustering: Unsupervised classification/clustering algorithm

ISODATA stands for Iterative Self-Organizing Data Analysis Techniques. 

ISODATA_clustering includes an algorithm which automatically creates groups of data (clusters) by means of an iterative process. During the process, similar clusters are merged and those with large standard deviations are split. 
The centroids of the clusters are randomly selected. Then, the euclidean distance among the centroids and the SD of the euclidean distance among each point of the cluster and its centroid are calculated. These two parameters are used to establish if some conditions defined by the user are reached and, therefore, the process of clustering is finished.

A dataset (dataset4isodata.RData) is provided for testing the function.
