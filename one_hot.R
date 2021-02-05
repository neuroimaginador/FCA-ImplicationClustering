one_hot <- function(cl) {
  
  if (is.factor(cl)) {
    cl_num <- as.numeric(cl)
  } else {
    
    cl_num <- cl
    
  }
  A <- matrix(0, nrow = length(cl), ncol = max(cl_num))
  ii <- seq_along(cl)
  A[(cl_num - 1) * length(cl_num) + ii] <- 1
  
  colnames(A) <- paste0(cur_column(), " = ",
                        as.character(levels(cl)))
  
  return(A)
  
}
