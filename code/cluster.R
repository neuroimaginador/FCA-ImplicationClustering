# This function finds the medoids given a distance matrix and a 
# vector representing the clustering
find_representatives <- function(D, cl) {
  
  D <- as.matrix(D)
  n_clusters <- max(cl)
  res <- vector(mode = "integer", length = n_clusters)
  
  for (i in seq(n_clusters)) {
    
    idx <- which(cl == i)
    
    if (length(idx) > 1) {
      
      
      tmp <- colMeans(D[idx, idx])
      
      res[i] <- idx[which.min(tmp)[1]]
      
    } else {
      
      res[i] <- idx
      
    }

  }
  
  return(res)
  
}

# This function clusters implications, given the
# dissimilarity matrix and the number k of clusters.
# All methods are equivalent, only difference is
# performance in larger datasets.
cluster_implications <- function(diss, k, 
                                 method = "diana") {
  
  res <- list(clustering = NULL,
              representatives = NULL,
              mean_diss = NULL,
              max_diss = NULL,
              inter = NULL)
  
  # Several methods have been implemented, all of them 
  # equivalent to PAM, with the same idea of partitioning, 
  # but performing slightly better with large datasets.
  if (method == "pam") {
    
    cl <- cluster::pam(diss, k = k)
    res$clustering <- cl$clustering
    res$representatives <- cl$id.med
    
  }
  
  if (method == "fanny") {
    
    if (k > 1) {
      
      cl <- cluster::fanny(x = diss, 
                           k = k, 
                           memb.exp = 1.1,
                           maxit = 1000)
      res$clustering <- cl$clustering

    } else {
      
      res$clustering <- rep(1, attr(diss, "Size"))
      
    }
    
    res$representatives <- find_representatives(diss,
                                                res$clustering)
    
  }
  
  if (method == "agnes") {
    
    cl <- cluster::agnes(diss, method = "ward")
    res$clustering <- cutree(cl, k = k)
    res$representatives <- find_representatives(diss,
                                                res$clustering)
    
  }
  
  if (method == "diana") {
    
    cl <- cluster::diana(diss)
    res$clustering <- cutree(cl, k = k)
    res$representatives <- find_representatives(diss,
                                                res$clustering)
    
  }
  
  D <- as.matrix(diss)
  
  # We return the mean and maximum intracluster dissimilarity 
  for (i in seq(k)) {
    
    ii <- which(res$clustering == i)
    Dtmp <- D[res$representatives[i], ii]
    
    res$mean_diss[i] <- mean(Dtmp)
    res$max_diss[i] <- max(Dtmp)
    
  }
  
  # And, for each pair of clusters, find its minimum, mean 
  # and maximum dissimilarity
  res$inter <- tibble(i = NULL, j = NULL, 
                      min_d = NULL, mean_d = NULL,
                      max_d = NULL)
  
  if (k > 1) {
    
    for (i in seq(k - 1)) {
      
      for (j in seq(i + 1, k)) {
        
        foo <- D[which(res$clustering == i), which(res$clustering == j)]
        
        if (!is.vector(foo)) {
          
          foo <- foo %>% 
            apply(1, min)
          
        }
        
        res$inter <- rbind(res$inter,
                           tibble(i = i, j = j, 
                                  min_d = min(foo), 
                                  mean_d = mean(foo, trim = 0.1),
                                  max_d = max(foo)))
        
        
      }
      
    }
    
  } else {
    
    res$inter <- tibble(i = 1, j = 1, 
                        min_d = 0, mean_d = 0,
                        max_d = 0)
    
  }

  return(res)
  
}
