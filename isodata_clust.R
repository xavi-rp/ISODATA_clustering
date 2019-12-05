
##############################################################################
###########                                                       ############
###########                 ISODATA Clustering                    ############
###########                                                       ############
###########   Unsupervised classification/clustering algorithm    ############
###########                                                       ############
##############################################################################
#
#
# isodata_clust.R
#
# Created on: Autumn 2019
#
# Created by: Xavier Rotllan-Puig (xavier.rotllan.puig@aster-projects.cat)
#
# Description: ISODATA stands for Iterative Self-Organizing Data Analysis 
#              Techniques. This script includes an algorithm which 
#              automatically creates groups of data (clusters) by means of 
#              an iterative process. During the process, similar
#              clusters are merged and those with large standard deviations 
#              are split. 
#              The centroids of the clusters are randomly selected. Then,
#              the euclidean distance among the centroids and the SD of the 
#              euclidean distance among each point of the cluster and its 
#              centroid are calculated. These two parameters are used to 
#              establish if some conditions defined by the user are reached 
#              and, therefore, the process of clustering is finished.
#
#
# ------------------------------------------

#' @param data_ini Data set
#' @param nclust_ini Number of clusters to start the process (first iteration)
#' @param max_num_clustrs Maximum number of cluster allowed
#' @param max_SD_intraclust Maximum standard deviation of euclidean distances among points in the cluster and its centroid
#' @param min_dist_interclust Minimum euclidean distance among centroids
#' @param max_iter Maximum number of iterations
#' @param path2saveResults Where to save results


## Packages ##
library(stats)
library(parallel)
library(rrcov)
library(dplyr)
##


isodata_clustering <- function(data_ini,
                               nclust_ini = 15,
                               max_num_clustrs = 500,
                               max_SD_intraclust,
                               min_dist_interclust,
                               max_iter = 100,
                               path2saveResults = NULL){
  
  if(is.null(path2saveResults)) path2saveResults <- getwd()
  
  ## step1: selecting random centroids of the clusters 
  iter_num <- 1  
  nvars <- length(data_ini)
  
  print(paste0("iteration # ", iter_num, " (max allowed = ", max_iter, ")"))
  print(paste0("Starting at ", Sys.time()))
  
  clust_centr_ini <- data_ini[sample(nrow(data_ini), nclust_ini), ]  #centroids (random)
  
  data_ini_noCentr <- data_ini[!rownames(data_ini) %in% rownames(clust_centr_ini), ]
  
  ## step2: calculating the closest centroid (assigning to a cluster) 
  # parallelizing
  cl <- makeCluster(3)
  clusterExport(cl, c("clust_centr_ini", "nvars"), envir=environment())
  clsts <- parRapply(cl = cl, 
                     x = data_ini_noCentr[, c(1:nvars)], 
                     FUN = function(y) rownames(clust_centr_ini)[which.min(apply(clust_centr_ini, 1, function(x) dist(rbind(y, x))))]
  )
  stopCluster(cl)
  data_ini_noCentr$closest <- as.vector(clsts)
  
  stuff2save <- c("clust_centr_ini", "data_ini_noCentr", "iter_num")
  save(list = stuff2save, file = paste0(path2saveResults, "/results_isodata_iter", iter_num, ".RData"))
  
  
  
  repeat{
    
    iter_num <- iter_num + 1
    print(paste0("iteration # ", iter_num, " (max allowed = ", max_iter, ")"))
    print(paste0("Starting at ", Sys.time()))
    
    ## Step3: SD distance intra-cluster and Distance inter-cluster 
    # if SD of distances to the centroid of the cluster is > max_SD_intraclust, split it in two clusters. How? Assigning two new random centroids?? 
    # if distance among centroids is smaller than a threshold --min_dist_interclust--, merge clusters
    
    ## Distance to closest (centroid)
    cl <- makeCluster(3)
    clusterExport(cl, c("clust_centr_ini", "nvars"), envir=environment())
    
    clsts <- parRapply(cl = cl, 
                       x = data_ini_noCentr, 
                       FUN = function(y) {
                         y1 <- clust_centr_ini[rownames(clust_centr_ini) == y[5], 1:nvars]
                         y2 <- dist(rbind(y1, y[1:nvars]))
                         return(y2)
                       }
    )
    ## verification: dist(rbind(data_ini_noCentr[1, 1:4], clust_centr_ini[rownames(clust_centr_ini) == 6662183, ]))
    stopCluster(cl)
    print(paste0("End of calculating distance to centroids at ", Sys.time()))
    
    data_ini_noCentr$dist2closest <- as.vector(clsts)
    
    
    ## Distance Inter-Clusters (centroids; if smaller than a threshold --min_dist_interclust--, merge clusters)
    dist_interclust <- dist(clust_centr_ini)
    #quantile(dist_interclust, seq(0, 1, 0.05))
    #percentile <- ecdf(dist_interclust)
    #if((round(percentile(min_dist_interclust), 2) * 100) > 5) stop("Your minimum distance intercluster is bigger than 5th percentile.\nThis probably means a lot of cluster mergings. Are you happy with that???") 
    
    small_centroids <- dist_interclust[dist_interclust < min_dist_interclust] 
    dist_interclust <- as.data.frame(as.matrix(dist_interclust))
    
    clust_mergd <- c()
    if(sum(small_centroids) > 0){
      for (sml_ctrds in sort(small_centroids)){ #sorting them we make sure that are taken into account the closest pair of clusters 
        clusts2merge <- rownames(which(dist_interclust == sml_ctrds, arr.ind = TRUE))
        clust_mergd <- c(clust_mergd, clusts2merge)
        
        if(all(clusts2merge %in% rownames(clust_centr_ini))){
          centroids_back <- clust_centr_ini[rownames(clust_centr_ini) %in% clusts2merge, ]
          centroids_back$closest <- NA
          centroids_back$dist2closest <- NA
          data_ini_noCentr <- rbind(data_ini_noCentr, centroids_back)
          data_ini_noCentr <- data_ini_noCentr[order(as.numeric(rownames(data_ini_noCentr))), ]
          clust_centr_ini <- clust_centr_ini[!rownames(clust_centr_ini) %in% clusts2merge, ]
          
          sml_cluster2merge <- data_ini_noCentr[data_ini_noCentr$closest %in% clusts2merge, ]
          sml_cluster2merge <- sml_cluster2merge[sample(nrow(sml_cluster2merge), 1), ][, 1:nvars]
          clust_centr_ini <- rbind(clust_centr_ini, sml_cluster2merge)
          data_ini_noCentr <- data_ini_noCentr[!rownames(data_ini_noCentr) %in% rownames(sml_cluster2merge), ]
          
        }
      }
    }
    print(paste0("End of distance inter-clusters at ", Sys.time()))
    
    
    ## SD intra-cluster (if bigger than a threshold, split the cluster in two)
    #stats_intraclustr <- as.data.frame(data_ini_noCentr[!is.na(data_ini_noCentr$closest), ] %>% group_by(closest) %>% summarise_at(.vars = "dist2closest", .funs = c("mean_intraclust" = mean, "SD_intraclust" = sd, "max_intraclust" = max, "min_intraclust" = min)))
    stats_intraclustr <- as.data.frame(data_ini_noCentr[!is.na(data_ini_noCentr$closest), ] %>% group_by(closest) %>% summarise_at(.vars = "dist2closest", .funs = c("SD_intraclust" = sd)))
    sd_intraclustr <- stats_intraclustr[, colnames(stats_intraclustr) %in% c("closest", "SD_intraclust")]
    big_clusters <- sd_intraclustr$closest[which(sd_intraclustr$SD_intraclust >= max_SD_intraclust)]
    big_clusters <- big_clusters[big_clusters %in% rownames(clust_centr_ini)]
    
    if(length(big_clusters) > 0){
      centroids_back <- clust_centr_ini[rownames(clust_centr_ini) %in% big_clusters, ]
      centroids_back$closest <- NA
      centroids_back$dist2closest <- NA
      data_ini_noCentr <- rbind(data_ini_noCentr, centroids_back)
      data_ini_noCentr <- data_ini_noCentr[order(as.numeric(rownames(data_ini_noCentr))), ]
      
      clust_centr_ini <- clust_centr_ini[!rownames(clust_centr_ini) %in% big_clusters, ]
      
      for (bg_clstr in big_clusters) {
        bg_cluster2split <- data_ini_noCentr[data_ini_noCentr$closest %in% bg_clstr, ]
        bg_cluster2split <- bg_cluster2split[sample(nrow(bg_cluster2split), 2), ][, 1:nvars]
        clust_centr_ini <- rbind(clust_centr_ini, bg_cluster2split)
        data_ini_noCentr <- data_ini_noCentr[!rownames(data_ini_noCentr) %in% rownames(bg_cluster2split), ]
      }
    }
    print(paste0("End of SD inter-centroids at ", Sys.time()))
    
    
    ## Step4: Calculating the closest centroid (assigning to a cluster) with the new centroids
    # parallelizing
    cl <- makeCluster(3)
    clusterExport(cl, c("clust_centr_ini", "nvars"), envir=environment())
    clsts <- parRapply(cl = cl, 
                       x = data_ini_noCentr[, c(1:nvars)], 
                       FUN = function(y) rownames(clust_centr_ini)[which.min(apply(clust_centr_ini, 1, function(x) dist(rbind(y, x))))]
    )
    stopCluster(cl)
    data_ini_noCentr$closest <- as.vector(clsts)
    
    stuff2save <- c("clust_centr_ini", "data_ini_noCentr", "iter_num", "clust_mergd", "big_clusters")
    save(list = stuff2save, file = paste0(path2saveResults, "/results_isodata_iter", iter_num, ".RData"))
    
    print(paste0("End of iteration #", iter_num, " at ", Sys.time()))
    
    
    ## Conditions to stop the ISODATA clustering 
    if(iter_num == max_iter){
      print("maximum number of iterations reached")
      break  
    } 
    if(is.null(clust_mergd) & length(big_clusters) == 0){
      print("all clusters are within defined parameters (max_SD_intraclust and min_dist_interclust)")
      break  
    }
    if(nrow(clust_centr_ini) >= max_num_clustrs){
      print("maximum number of clusters reached")
      break
    }
  }
  
  print(paste0("Writing results..."))
  
  ## Putting results in a table 
  iters <- paste0("iter", 1:iter_num)
  wilks_clust_df <- as.data.frame(matrix(nrow = 0, ncol = 0))
  for (i in iters) {
    load(paste0(path2saveResults, "/results_isodata_", i, ".RData"), verbose = FALSE)
    frmla <- as.formula(paste0("closest ~ ", paste0(colnames(data_ini_noCentr)[1:nvars], collapse = " + ")))
    wilks_formula <- Wilks.test(frmla, data = data_ini_noCentr)
    conon_corr <- cancor(data_ini_noCentr[, 1], data_ini_noCentr[, 2:nvars])$cor
    wilks_clust_df <- rbind(wilks_clust_df, data.frame(i, nrow(clust_centr_ini), wilks_formula$statistic, conon_corr))
  }
  
  names(wilks_clust_df) <- c("Iteration", "NumberOfClusters", "WilksLambda", "CanonicalCorrel")
  save(list = "wilks_clust_df", file = paste0(path2saveResults, "/results_isodata_WilksDF", ".RData"))
  
  
  ## Returning results
  
  data_ini_noCentr <- data_ini_noCentr[, c(1:(nvars + 1))]
  
  names(data_ini_noCentr)[names(data_ini_noCentr) %in% "closest"] <- "centroid"
  clust_centr_ini$centroid <- rownames(clust_centr_ini)
  
  clust_diff <- setdiff(unique(data_ini_noCentr$centroid), clust_centr_ini$centroid)
  clust_diff <- clust_diff[!is.na(clust_diff)]
  if(length(clust_diff) > 0){
    data_ini_noCentr[rownames(data_ini_noCentr) %in% clust_diff, ]$centroid <- rownames(data_ini_noCentr[rownames(data_ini_noCentr) %in% clust_diff, ])
  }
  
  all_data <- rbind(data_ini_noCentr, clust_centr_ini)
  all_data <- all_data[order(as.numeric(rownames(all_data))), ]
  
  print(paste0("End of ISODATA clustering at ", Sys.time()))
  
  
  return(list(isodta_wilks = wilks_clust_df, 
              centroids_isodta = sort(unique(all_data$centroid)),
              data_ini_centroids = all_data,
              iterations = iter_num,
              number_clusters = length(unique(all_data$centroid))))
  
}


