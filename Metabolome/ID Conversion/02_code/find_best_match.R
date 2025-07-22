library(stringdist)

find_best_match <- function(query, candidates, candidate_ids, max_dist = 2) {
  dists <- stringdist::stringdist(query, candidates, method = "osa")
  min_dist <- min(dists)
  
  if (min_dist <= max_dist) {
    match_index <- which.min(dists)
    return(candidate_ids[match_index])
  } else {
    return(NA)
  }
}
