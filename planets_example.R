#' # planets dataset - Clustering of Implications


#' ## Setup
#' 
#' First, we load the needed libraries
# Suppress excess info
options(dplyr.summarise.inform = FALSE)
options(tidyverse.quiet = TRUE)
library("tidyverse")
library("cluster")
library("fcaR")

#' Load other functions needed to compute implication
#' distances and perform clustering.
source("./distances_registry.R")
source("./implication_distance.R")
source("./cluster.R")
source("./one_hot.R")
source("./count_attributes.R")

#' Auxiliary functions
imps2str <- function(imp) {
  sapply(seq(imp$cardinality()), 
         function(j) 
           fcaR:::.implication_to_string(
             imp[j]$get_LHS_matrix(), 
             imp[j]$get_RHS_matrix(), 
             imp$get_attributes()))
}

#' 
#'  
#' ## The planets dataset
#' 
#' We'll use the planets dataset, which is included
#' in the fcaR package, to show how to cluster 
#' implications.
#' 
#' First, load and print the dataset.
#' 
## --------------------------------
library(fcaR)
fc_planets <- FormalContext$new(planets)

fcaR:::.print_binary(planets) %>% 
  kableExtra::kable(align = "c")

#' Let us find the Duquenne-Guigues basis of
#' implications and make a copy in variable "imps"
#' 
## --------------------------------
fc_planets$find_implications()
imps <- fc_planets$implications$clone()
imps

#' 
#' We can compute the dissimilarity between each pair
#' of implications, as
#' $\text{diss}(P\to Q, R\to T) := d(P, R) + d(P^+, R^+)$
#' where $d$ is the Manhattan distance.
#' 
## --------------------------------------------
Distances <- impl_dist(imps, method = "Manhattan")
Dimps <- Distances$left + Distances$closure
# DD is the dissiilarity matrix
DD <- as.matrix(Dimps)
kableExtra::kable(DD, align = "c", row.names = TRUE)

#' 
#' We can use a clustering algorithm to find $k = 2$
#' clusters
## --------------------------------------------
cluster <- cluster_implications(Dimps, k = 2, "pam")

# The central implications
imps[cluster$representatives]

#' 
#' The implications in the first cluster are
#' 
## --------------------------------------------
imps[cluster$clustering == 1]

#' 
#' and in the second one:
#' 
## --------------------------------------------
imps[cluster$clustering == 2]

#' Note that this cluster is composed by all implications
#' in which all attributes appear. Let us remove them
#' to cluster the other implications properly:
imps <- imps[imps$support() > 0]
print(imps)

#' 
#' Let us compute the dissimilarity matrix
#' 
## --------------------------------------------
imps
Distances <- impl_dist(imps, method = "Manhattan")
Dimps <- Distances$left + Distances$closure
DD <- as.matrix(Dimps)
kableExtra::kable(DD, align = "c", row.names = TRUE)

#' 
#' And cluster into two groups of implications
#' 
## --------------------------------------------
cluster <- cluster_implications(Dimps, k = 2, "pam")

# The central implications are
imps[cluster$representatives]

#' 
#' Now, the implications in the first cluster are:
#' 
## --------------------------------------------
imps[cluster$clustering == 1]


#' 
#' and the ones in the second one are:
#' 
## --------------------------------------------
imps[cluster$clustering == 2]

#' 
#' This way, we can tag the clusters according to the
#' closed sets that arise from the central implications
#' 
## ----results="asis"--------------------------
for (i in seq(k)) {
  
  cat(" - Cluster", i, ": ")
  imp_med <- imps[cluster$representatives[i]]$clone()
  M <- fcaR:::.union(imp_med$get_LHS_matrix(),
                     imp_med$get_RHS_matrix())
  S <- SparseSet$new(attributes = fc_planets$attributes,
                     M = M)
  print(S)
  cat("\n")
  
}

#' 
#' This shows that the dataset has two main clusters:
#' 
#' 1) Planets near the Sun and therefore small,
#' 
#' 2) Planets far from the Sun and therefore with satellites. 
#' 
#' __Graphical representation of the implication space__
#' 
#' We can use multidimensional scaling to plot the 
#' implications and their clusters:
imp_mds <- MASS::isoMDS(d = Dimps, k = 2, trace = FALSE) 
df <- imp_mds$point %>% as.data.frame() %>% 
  as_tibble() %>%  
  set_names(c("x", "y")) %>% 
  mutate(cluster = cluster$clustering,
         imp = imps2str(imps))

library(ggrepel)
ggplot(df, aes(x = x, y = y, color = cluster, label = imp)) + 
  geom_label_repel() +
  geom_point() +
  theme_void() + 
  theme(legend.position = "none",
        axis.line = element_line(size = 0))

#' 
#' __What happens if we consider only the left-hand side
#' of the implications when computing the dissimilarity
#' matrix?__
#' 
## --------------------------------------------
Distances <- impl_dist(imps, method = "Manhattan")
Dimps <- Distances$left
DD <- as.matrix(Dimps)
rownames(DD) <- seq(imps$cardinality())
kableExtra::kable(DD, align = "c", row.names = TRUE)

#' 
#' We can use the PAM clustering algorithm to find 
#' 2 clusters:
#' 
## --------------------------------------------
cluster <- cluster_implications(Dimps, k = 2, "pam")

# The central implications in this case are:
imps[cluster$representatives]

#' 
#' The implications in the first cluster:
#' 
## --------------------------------------------
imps[cluster$clustering == 1]


#' 
#' and in the second one:
#' 
## --------------------------------------------
imps[cluster$clustering == 2]

#' This clustering is almost random, due to the uninformative 
#' dissimilarity matrix used.
