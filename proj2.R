## proj2
make_households <- function(n, hmax = 5, seed = NULL, shuffle = TRUE, return_sizes = FALSE) {
  ## Generate a vector h of household IDs (integers) of length n.
  ## Household sizes are uniformly distributed between 1 and hmax.
  ## Includes safety checks and dynamic batch sampling for efficiency.
  ## Optional: reproducible seed, random shuffle, and returning size summary.
  ## Returns integer vector h where same number indicates same household.
  
  ## If seed is provided, set RNG seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  ## Convert to integer explicitly
  n <- as.integer(n)
  hmax <- as.integer(hmax)
  
  ## --- Initialize ---
  sizes <- integer(0)   ## store sampled household sizes
  total <- 0            ## total number of people generated so far
  
  ## Dynamic batch size: ensures efficiency and prevents over-allocation
  batch <- min(max(50, ceiling(n / 10)), 10000)
  
  ## --- Sampling loop ---
  ## Keep sampling household sizes until we have enough to cover n individuals
  while (total < n) {
    new_sizes <- sample.int(hmax, size = batch, replace = TRUE)  ## uniform sizes 1:hmax
    sizes <- c(sizes, new_sizes)   ## append new sampled sizes
    total <- total + sum(new_sizes)## update total count
  }
  
  ## --- Adjust to exactly n ---
  cs <- cumsum(sizes)                 ## cumulative total of individuals
  last <- which(cs >= n)[1]           ## find last needed household
  sizes <- sizes[1:last]              ## keep only needed households
  sizes[last] <- sizes[last] - (cs[last] - n)  ## adjust last one to fit exactly
  
  ## If the last household ends up with size 0, remove it
  if (sizes[length(sizes)] == 0) {
    sizes <- sizes[-length(sizes)]
  }
  
  ## --- Construct household ID vector ---
  h <- rep.int(seq_along(sizes), times = sizes)  ## repeat household id per member
  
  ## Sanity check: ensure final length matches target n
  if (length(h) != n) stop("Internal error: generated h has wrong length")
  
  ## --- Optional shuffle ---
  ## Randomize order of individuals if requested
  if (isTRUE(shuffle)) {
    h <- sample(h, length(h), replace = FALSE)
  }
  
  ## --- Return ---
  ## Optionally return both household vector and size summary
  if (isTRUE(return_sizes)) {
    return(list(h = h, sizes = as.integer(table(h))))
  } else {
    return(h)
  }
}


