
#-----------#
# Functions #
#-----------#

SetMembers <- function(x) {
  # Returns all unique members of a list and flattens it.
  as.numeric(unique(unlist(x)))
}

FilterEmpty <- function(l) {
  # Removes empty vectors within a list.
  l[lapply(l,length)>0]
} 


pairwise_analysis <- function(comp_set,urn_set) {
  # """ Performs pairwise hypergeometric test on two lists """
  # Generate complete, intersecting set
  i_set <- intersect(SetMembers(comp_set), SetMembers(urn_set))
  
  # Filter each set, retaining only elements in common.
  comp_set <- FilterEmpty(lapply(comp_set, function(c) c[c %in% i_set]))
  urn_set <- FilterEmpty(lapply(urn_set, function(c) c[c %in% i_set]))
  
  # (q) - Pairwise intersection matrix (# intersecting)!
  q <- lapply(comp_set, function(c1) {
    lapply(urn_set, function(c2) length(intersect(c1, c2)))
  })
  
  # (m) - Number of white balls in the urn.
  m <- lapply(comp_set, length)
  
  # (n) - Number of black balls in the urn.
  n <- lapply(m, function(x) length(i_set) - x) 
  
  # (k) - Number of balls drawn from the urn.
  k <- lapply(urn_set, length)
  
  # Perform the hypergeometric test!
  results <- sapply(names(comp_set), USE.NAMES=T, function(c1) {
    sapply(names(urn_set), USE.NAMES=T, function(c2) {
      
      phyper(q = q[[c1]][[c2]],
             m = m[[c1]],
             n = n[[c1]],
             k = k[[c2]],
             lower.tail = F,
             log.p = T
      )
    })
  })
  # Replace -Inf results with a max.
  results[results == - Inf] <- min(results[results != - Inf])
  # Return the resutls
  results
}