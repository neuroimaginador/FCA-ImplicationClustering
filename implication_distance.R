impl_dist <- function(imps, 
                      method = distancesRegistry$get_entry_names()) {
  
  method <- match.arg(method)
  distance <- distancesRegistry$get_entry(method)$fun
  
  LHS <- imps$get_LHS_matrix()
  RHS <- imps$get_RHS_matrix()
  
  closures1 <- fcaR:::.union(LHS, RHS)
  
  dl <- distance(LHS) %>% as.dist()
  dr <- distance(RHS) %>% as.dist()
  dlr <- distance(closures1) %>% as.dist()
  
  return(list(left = dl, right = dr, closure = dlr))
  
}

