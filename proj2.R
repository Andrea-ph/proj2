#proj2
make_households <- function(n, hmax = 5, seed = NULL) {
  # Creates vector h of household IDs (integers)
  # Sizes are uniformly distributed between 1 and hmax.
  # Robust trimming ensures exactly n people and no size-0 households.
  # The final vector is shuffled to randomize assignments.
  
  if (!is.null(seed)) set.seed(seed)                 # Set random seed if provided
  
  sizes <- integer(0)                                # Initialize vector to store sampled household sizes
  while (sum(sizes) < n) {                           # Keep adding sampled sizes until total >= n
    sizes <- c(sizes, sample.int(hmax, size = 100, replace = TRUE))  # Sample 100 household sizes at a time
  }
  cs <- cumsum(sizes)                                # Compute cumulative sum to track total population
  last <- which(cs >= n)[1]                          # Find the last household that reaches or exceeds n
  sizes <- sizes[1:last]                             # Keep only the households needed to reach n
  sizes[last] <- sizes[last] - (cs[last] - n)       # Adjust last household size to make total exactly n
  
  h <- rep.int(seq_along(sizes), times = sizes)      # Repeat household IDs according to household sizes
  h <- sample(h, length(h), replace = FALSE)        # Shuffle to assign individuals randomly
  return(h)                                         # Return the household membership vector
}
