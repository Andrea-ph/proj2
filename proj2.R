###########################################################
# Practical 2 — Social Structure in SEIR Models 
###########################################################
####proj2 - Group 33 - Extended Statistical Programming ###
#### Group members as below ################################
#### Shuo Li (s2795688), Zhe Zhu (s2841606), Antrea Filippou (s2766374)
#### Contributions as below ################################
#### Shuo Li: xx (xx%) ###
#### Zhe Zhu: xx (xx%) ###
#### Antrea Filippou: xx (xx%) ###
############################################################


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

# ===============================================
# Function: get.net(beta, nc = 15, h = NULL)
# Builds a symmetric network of regular (non-household) contacts.
# beta is the n vector of βi value for each person
# nc is the average number of contacts per person
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

# 
# ============================================================
# Implements:
# - E->I with daily probability gamma, I->R with daily probability delta.
# - S->E via three independent channels (household, contact network, random mix).
# Independence => multiply "no-infection" probabilities from each channel.
# Random mixing uses the standard small-p exponential approximation:
#   prod_i (1 - c_mix * beta_S * beta_i)  ≈  exp(- c_mix * beta_S * sum(beta_I))
# (Exact version is provided below commented; slower for large n.)

nseir <- function(beta, h, alink,
                  alpha = c(0.1, 0.01, 0.01), # c(alpha_h, alpha_c, alpha_r)
                  delta = 0.2,
                  gamma = 0.4,
                  nc    = 15,
                  nt    = 100,
                  pinf  = 0.005) {
  n     <- length(beta)
  stopifnot(length(h) == n, length(alink) == n)
  
  bbar  <- mean(beta)
  tvec  <- seq_len(nt)
  S_CODE <- 1L; E_CODE <- 2L; I_CODE <- 3L; R_CODE <- 4L
  
  # Initial states
  state <- rep.int(S_CODE, n)
  nI0   <- max(1L, round(pinf * n))
  initI <- sample.int(n, size = nI0, replace = FALSE)
  state[initI] <- I_CODE
  
  # Household membership index lists
  H_ids <- unique(h)
  HH <- vector("list", length(H_ids))
  names(HH) <- as.character(H_ids)
  for (hid in H_ids) HH[[as.character(hid)]] <- which(h == hid)
  
  # Output time series
  S_daily <- integer(nt)
  E_daily <- integer(nt)
  I_daily <- integer(nt)
  R_daily <- integer(nt)
  
  # Random mixing constant (shares nc with network, per brief)
  c_mix <- alpha[3] * nc / (bbar^2 * (n - 1))
  
  for (tt in tvec) {
    indS <- which(state == S_CODE)
    indE <- which(state == E_CODE)
    indI <- which(state == I_CODE)
    indR <- which(state == R_CODE)
    
    # Record totals
    S_daily[tt] <- length(indS)
    E_daily[tt] <- length(indE)
    I_daily[tt] <- length(indI)
    R_daily[tt] <- length(indR)
    
    # E -> I
    if (length(indE)) {
      uEI <- runif(length(indE))
      toI <- indE[uEI < gamma]
      if (length(toI)) state[toI] <- I_CODE
    }
    
    # I -> R
    indI <- which(state == I_CODE)
    if (length(indI)) {
      uIR <- runif(length(indI))
      toR <- indI[uIR < delta]
      if (length(toR)) state[toR] <- R_CODE
    }
    
    # New infections S -> E (if both S and I exist)
    indS <- which(state == S_CODE)
    indI <- which(state == I_CODE)
    if (length(indS) && length(indI)) {
      # Household exposures: for each S, count Infectious in same HH
      I_by_hh <- integer(length(H_ids)); names(I_by_hh) <- as.character(H_ids)
      I_tab <- table(h[indI])
      if (length(I_tab)) I_by_hh[names(I_tab)] <- as.integer(I_tab)
      H_I_per_S <- I_by_hh[as.character(h[indS])]
      Pnh <- (1 - alpha[1])^H_I_per_S
      
      # Contact-network exposures: for each S, count Infectious neighbors
      I_flag <- logical(n); I_flag[indI] <- TRUE
      C_I_per_S <- integer(length(indS))
      for (k in seq_along(indS)) {
        j <- indS[k]
        nbrs <- alink[[j]]
        if (length(nbrs)) C_I_per_S[k] <- sum(I_flag[nbrs])
      }
      Pnc <- (1 - alpha[2])^C_I_per_S
      
      # Random mixing exposures: exponential approximation (fast)
      sum_beta_I <- sum(beta[indI])
      Pnr <- exp(- c_mix * beta[indS] * sum_beta_I)
      
      # --- Exact (slower) version for reference ---
      # Pnr <- exp(colSums(log1p(- c_mix * outer(beta[indS], beta[indI]))))
      
      # Combine independent channels
      P_not_infected <- Pnh * Pnc * Pnr
      P_infected     <- 1 - P_not_infected
      
      # Apply infections
      uSE <- runif(length(indS))
      toE <- indS[uSE < P_infected]
      if (length(toE)) state[toE] <- E_CODE
    }
  }
  
  out <- list(S = S_daily, E = E_daily, I = I_daily, R = R_daily, t = tvec)
  return(out)
}

#
# =======================================

plot_nseir <- function(sim, main = "SEIR with Households & Contacts") {
  stopifnot(all(c("S","E","I","R","t") %in% names(sim)))
  mat <- cbind(S = sim$S, E = sim$E, I = sim$I, R = sim$R)
  op <- par(mar = c(4.2, 4.5, 3.5, 1.2))
  on.exit(par(op))
  matplot(sim$t, mat, type = "l", lwd = 2,
          xlab = "Day", ylab = "Population count", main = main,
          lty = 1) # make legend lty consistent
  legend("right", inset = 0.01, lwd = 2, col = 1:4, lty = 1,
         legend = c("S","E","I","R"), bg = "white", cex = 0.9)
}

# 
# =========================================================
# Scenarios:
# A) Full model, beta ~ U(0,1)
# B) Random mixing only (alpha_h = alpha_c = 0, alpha_r = 0.04)
# C) Full model, constant beta = mean(beta)
# D) Random mixing + constant beta
# For fair comparison, we use one RNG seed at the start, but do not
# re-seed inside each simulation to avoid accidental identical draws.

run_four_scenarios <- function(n = 1000, nt = 150, hmax = 5, nc = 15,
                               alpha_full = c(0.1, 0.01, 0.01),
                               alpha_random_only = c(0, 0, 0.04),
                               delta = 0.2, gamma = 0.4, pinf = 0.005,
                               seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  
  betaA <- runif(n, 0, 1)
  h     <- make_households(n, hmax = hmax)  # uses current RNG stream
  alink <- get.net(betaA, nc = nc, h = h)
  
  simA <- nseir(betaA, h, alink,
                alpha = alpha_full, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  simB <- nseir(betaA, h, alink,
                alpha = alpha_random_only, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  betaC <- rep(mean(betaA), n)
  alinkC <- get.net(betaC, nc = nc, h = h)
  
  simC <- nseir(betaC, h, alinkC,
                alpha = alpha_full, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  simD <- nseir(betaC, h, alinkC,
                alpha = alpha_random_only, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  # Helper to compose a compact title with peak I and final R
  mk_title <- function(lbl, sim) {
    peakI <- max(sim$I)
    finR  <- tail(sim$R, 1)
    paste0(lbl, "\npeak I = ", peakI, ", final R = ", finR)
  }
  
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3.5, 1.2))
  on.exit(par(op), add = TRUE)
  plot_nseir(simA, main = mk_title("A) Full model, beta ~ U(0,1)", simA))
  plot_nseir(simB, main = mk_title("B) Random mixing only (αh=αc=0, αr=0.04)", simB))
  plot_nseir(simC, main = mk_title("C) Full model, constant beta = mean(beta)", simC))
  plot_nseir(simD, main = mk_title("D) Random mixing + constant beta", simD))
  
  invisible(list(A = simA, B = simB, C = simC, D = simD,
                 betaA = betaA, betaC = betaC, h = h,
                 alinkA = alink, alinkC = alinkC))
}

# ===========================
# Example: run the scenarios
# ===========================
# (Reproducible with a single seed at the start; no reseeding inside.)

set.seed(42)
res <- run_four_scenarios(
  n = 1000,
  nt = 150,
  hmax = 5,
  nc = 15,
  alpha_full = c(0.1, 0.01, 0.01),
  alpha_random_only = c(0, 0, 0.04),
  delta = 0.2,
  gamma = 0.4,
  pinf = 0.005,
  seed = NULL  # already seeded above; leave NULL to avoid resetting inside
)

# ---------------------------
# Brief commentary (Step 5):
# ---------------------------
#  runs, scenarios with structured mixing (A, C) produce
# - slightly later and/or lower peaks than pure random mixing (B, D),
# - smaller final size when heterogeneity in beta is present (A vs C),
# consistent with the notes: variability and clustering suppress spread.
# Random mixing with elevated αr (B, D) tends to raise peak I and final R.
