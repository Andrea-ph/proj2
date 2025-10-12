make_households <- function(n, hmax, seed = NULL, shuffle = TRUE, return_sizes = FALSE) {
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
    result <- list(h = h, sizes = as.integer(table(h)))
  } else {
    result <- h
  }
  
  ## Print output for verification
  print(result)
  
  ## Return invisibly (so print happens only once)
  invisible(result)
}

make_households(n = 1000, hmax = 5, seed = 123, return_sizes = TRUE)

# = Step 2: Build non-household contact network =
# ===============================================
# Function: get.net(beta, nc = 15, h = NULL)
# Builds a symmetric network of regular (non-household) contacts.
# If household vector h is provided, excludes within-household links.
# Returns a list 'alink' where alink[[i]] are contacts of person i.

get.net <- function(beta, nc = 15, h = NULL) {

  n    <- length(beta)                      ## total number of individuals
  bbar <- mean(beta)                        ## mean of contact propensity (beta)
  if (n < 2L) return(vector("list", n))     ## no edges possible if fewer than 2 people

  cst  <- nc / (bbar^2 * (n - 1))           ## constant ensuring expected degree = nc
  alink <- vector("list", n)                ## initialize empty list for adjacency info

  HH <- NULL                                ## initialize household membership index
  if (!is.null(h)) {                        ## only build household index if h is provided
    H_ids <- unique(h)                      ## get all unique household IDs
    HH <- vector("list", length(H_ids))     ## create list to hold members of each household
    names(HH) <- as.character(H_ids)        ## name each list element by its household ID
    for (hid in H_ids) {                    ## loop over each household
      HH[[as.character(hid)]] <- which(h == hid)  ## store member indices of this household
    }
  }

  for (i in 1:(n - 1)) {                    ## loop over all individuals except last
    js <- (i + 1):n                         ## potential partners are those after i (avoid duplicates)
    if (!is.null(HH)) {                     ## if household info exists
      hid <- h[i]                           ## get household ID of person i
      hh_members <- HH[[as.character(hid)]] ## get all members in i’s household
      if (length(hh_members))               ## if household not empty
        js <- setdiff(js, hh_members)       ## remove household members from potential partners
    }
    if (length(js) == 0) next               ## skip if no eligible partners left

    p <- cst * beta[i] * beta[js]           ## edge probability for each (i,j) pair
    p[p < 0] <- 0; p[p > 1] <- 1            ## ensure probabilities remain in [0,1]

    u <- runif(length(js))                  ## draw uniform random numbers for each potential link
    keep <- which(u < p)                    ## select links where random number < probability

    if (length(keep)) {                     ## if any links are formed
      nbrs <- js[keep]                      ## get indices of connected individuals
      alink[[i]] <- c(alink[[i]], nbrs)     ## add neighbors to i’s adjacency list
      for (j in nbrs) {                     ## for each neighbor j
        alink[[j]] <- c(alink[[j]], i)      ## add i to j’s adjacency list (symmetry)
      }
    }
  }

  for (i in seq_len(n)) {                   ## clean adjacency list for each person
    if (length(alink[[i]]) > 1)             ## only if person has more than one neighbor
      alink[[i]] <- sort(unique(alink[[i]]))## remove duplicates and sort for consistency
  }

  return(alink)                             ## return final list of contact links
}
