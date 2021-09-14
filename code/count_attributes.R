# This function counts the number of times each
# attribute is present in the implications
# count_attributes_old <- function(imps) {
#   
#   LHS <- imps$get_LHS_matrix()
#   RHS <- imps$get_RHS_matrix()
#   # We are using the Duquenne-Guigues basis, so
#   # computing the closure of the LHS is easy
#   closure <- fcaR:::.union(LHS, RHS)
#   
#   count_left <- Matrix::rowSums(LHS)
#   count_right <- Matrix::rowSums(RHS)
#   count_closure <- Matrix::rowSums(closure)
#   
#   return(list(left = count_left,
#               right = count_right,
#               closure = count_closure))
#   
# }

count_attributes <- function(imps,
                                 only_left = FALSE,
                                 weight = FALSE,
                                 normalize = FALSE) {
  
  LHS <- imps$get_LHS_matrix()
  
  if (only_left) {
    
    M <- LHS
    
  } else {
    
    RHS <- imps$get_RHS_matrix()
    M <- fcaR:::.union(LHS, RHS)
    
  }
  
  att <- imps$get_attributes()
  
  if (weight) {
    
    size <- imps$size()
    sizes <- rowSums(size)
    
    M <- Matrix::t(Matrix::t(M) / sizes)
    
  }
  
  att_imp <- Matrix::rowSums(M) %>% set_names(att)
  
  if (normalize) {
    
    att_imp <- att_imp / sum(att_imp)
    
  }
  
  return(att_imp)
  
}
