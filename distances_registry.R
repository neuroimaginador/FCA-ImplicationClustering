##%######################################################%##
#                                                          #
####               Registro de distancias               ####
#                                                          #
##%######################################################%##

distancesRegistry <- registry::registry(
  registry_class = "distance_registry",
  entry_class = "distance_rule")

distancesRegistry$set_field("method",
                               type = "character",
                               is_key = TRUE,
                               index_FUN = registry::match_partial_ignorecase)

distancesRegistry$set_field("fun",
                               type = "function",
                               is_key = FALSE)
distancesRegistry$set_field("description",
                               type = "character",
                               is_key = FALSE)

##%######################################################%##
#                                                          #
####                     Distancias                     ####
#                                                          #
##%######################################################%##

# Manhattan
distancesRegistry$set_entry(method = "Manhattan",
                               fun = function(X) 
                                 vegan::vegdist(t(as.matrix(X)), 
                                                method = "manhattan", 
                                                binary = TRUE),
                               description = "Manhattan distance")


# jaccard
distancesRegistry$set_entry(method = "jaccard",
                            fun = function(X) 
                              vegan::vegdist(t(as.matrix(X)), 
                                             method = "jaccard", 
                                             binary = TRUE),
                            description = "Jaccard distance")

# cosine
.cosine <- function(X) {
  
  prod <- Matrix::t(X) %*% X
  normas <- Matrix::colSums(X)
  normas <- normas %*% t(normas)
  return(1 - prod / normas)
  
}

distancesRegistry$set_entry(method = "cosine",
                            fun = .cosine,
                            description = "Cosine distance")

