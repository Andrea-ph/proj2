#proj2
make_households <- function(n, hmax = 5, seed = NULL) {
## create a vector h of household IDs (integers)
## household sizes are uniformly sampled between 1 and hmax
## trimming ensures exactly n people; no size-0 households
## the final vector is shuffled to randomize assignments

  if (!is.null(seed)) set.seed(seed)              ## set RNG seed if provided

  sizes <- integer(0)                             ## initialize vector to store household sizes
  while (sum(sizes) < n) {                        ## keep sampling until total population >= n
    sizes <- c(sizes, sample.int(hmax, size = 100, replace = TRUE))  ## sample 100 household sizes at a time
  }
  cs <- cumsum(sizes)                             ## cumulative sum of household sizes
  last <- which(cs >= n)[1]                       ## index of last household needed to reach n
  sizes <- sizes[1:last]                          ## keep only necessary households
  sizes[last] <- sizes[last] - (cs[last] - n)    ## adjust last household to make total exactly n

  h <- rep.int(seq_along(sizes), times = sizes)   ## repeat household IDs according to household sizes
  h <- sample(h, length(h), replace = FALSE)     ## shuffle to randomize individual assignments
  h                                              ## return the household membership vector
}
