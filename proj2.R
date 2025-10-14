###########################################################
####proj2 - Group 33 - Extended Statistical Programming ###
#### Group members as below ################################
#### Shuo Li (s2795688), Zhe Zhu (s2841606), Antrea Filippou (s2766374)
#### Contributions as below ################################
#### Shuo Li: xx (xx%) ###
#### Zhe Zhu: xx (xx%) ###
#### Antrea Filippou: xx (xx%) ###
############################################################

## Code to simulate epidemic dynamics using a social-structured SEIR model.
## The challenge is to extend the basic SEIR model by introducing 
## realistic social structures into how infections spread.
## In the basic SEIR model, people are assumed to infect each other completely at random. 
## Here we assume that people live in households and also have their own regular contact
## networks, reflecting relatively fixed social interactions like households, workplaces, and friend groups.
## Note that n people are uniformly distributed into household of sizes between 1 and hmax=5.
## Mathematically, there are n people each in one of the following four states: 
## Susceptible (S), Exposed (E), Infectious (I), or Recovered (R).
## Transitions between these states occur with the following daily probabilities:
## S → E : as a result of infection by someone in state I. Several ways for a person i in State I to
##   infect a person j in state S:
##   1. Within-household infection occurs with a daily probability αh of person i infecting j.
##   2. Regular contact network based infection occurs with a daily probability αc of i infecting j.
##   3. Random background infection occurs with probability αr * nc * (βi * βj) /
##      (mean(β)^2 * (n - 1)), where βi represents 'sociability' (infectiousness) parameter for the ith person
##      and nc is the average number of contacts per person.
## E → I : with daily probability γ (infection incubation rate)
## I → R : with daily probability δ (recovery rate)
## The aim is to investigate how incorporating social structure changes 
## epidemic dynamic compared to purely random mixing.

hmax<-5  ## use household size maximum to be 5 by default 
n<-1000  ## our code work with any population sizes, here we test and develop with n=1,000

h <- rep(  ## repeat household IDs
  seq_along(household_sizes <- sample(1:hmax, ceiling(n/mean(1:hmax)), replace = TRUE)),  
  ## sample from an integer sequence ranging from 1 to hmax =5 to generate a vector 
  ## representing the household sizes. Repeated selection is allowed in sampling 
  ## to ensure the uniform distribution of household sizes.
  household_sizes  ## repeat IDs according to household size
)[1:n]  ## trim vector to length n

get.net <- function(beta, h, nc = 15) {  
## function to generate regular contact network
## please note that people in the same household are excluded from such contacts. 
  n <- length(beta)  ## total number of individuals
  if (n < 2L) return(vector("list", n))  ## return empty list if less than 2 people
  bbar <- mean(beta)  ## mean infectivity
  cst <- nc / (bbar^2 * (n - 1))  ## constant to match expected degree
  alink <- vector("list", n)  ## initialize adjacency list
  
  H_ids <- unique(h)  ## unique household IDs
  HH <- vector("list", length(H_ids))  ## list to hold household members
  names(HH) <- as.character(H_ids)  ## name each list element by household ID
  for (hid in H_ids) {   ## loop over each household
    HH[[as.character(hid)]] <- which(h == hid)  ## store indices of household members
  }
  
  for (i in 1:(n - 1)) {  ## loop over individuals except last
    js <- (i + 1):n   ## potential partners to avoid duplicates
    hid <- h[i]  ## get household ID of i
    hh_members <- HH[[as.character(hid)]]  ## get members of i's household
    js <- setdiff(js, hh_members)  ## remove household members from partners
    if (length(js) == 0) next  ## skip if no partners left
    
    p <- cst * beta[i] * beta[js] ## compute edge probabilities
    p[p < 0] <- 0  ## ensure probabilities >= 0
    p[p > 1] <- 1  ## ensure probabilities <= 1
    
    u <- runif(length(js))  ## draw uniform random numbers
    keep <- which(u < p)  ## keep edges where u < probability
    
    if (length(keep)) {  ## if any edges are kept
      nbrs <- js[keep]   ## get connected partners
      alink[[i]] <- c(alink[[i]], nbrs)  ## add neighbors to i's adjacency list
      for (j in nbrs) alink[[j]] <- c(alink[[j]], i) ## add i to each neighbor's adjacency list
    }
  }
  
  for (i in seq_len(n)) {  ## loop to clean adjacency list
    if (length(alink[[i]]) > 1)  ## only if person has multiple neighbors
      alink[[i]] <- sort(unique(alink[[i]]))  ## remove duplicates and sort
  }
  
  return(alink) ## return final adjacency list
}

#=========================================================
# SEIR simulator with households and contact network
#=========================================================
nseir <- function(beta, h, alink,                   # SEIR simulation function
                  alpha = c(0.1, 0.01, 0.01),      # alpha_h, alpha_c, alpha_r
                  delta = 0.2,                     # daily probability of recovery
                  gamma = 0.4,                     # daily probability of becoming infectious
                  nc = 15,                          # average number of contacts
                  nt = 100,                         # number of simulated days
                  pinf = 0.005) {                   # initial infection proportion
  
  n <- length(beta)                                 # total number of individuals
  stopifnot(length(h) == n, length(alink) == n)     # check vector lengths match
  
  bbar <- mean(beta)                                # mean infectivity
  tvec <- seq_len(nt)                               # vector of time steps
  
  S_CODE <- 1L; E_CODE <- 2L; I_CODE <- 3L; R_CODE <- 4L # numeric codes for states
  
  state <- rep.int(S_CODE, n)                       # initialize all as susceptible
  nI0 <- max(1L, round(pinf * n))                  # compute initial infected count
  initI <- sample.int(n, size = nI0, replace = FALSE) # randomly select initial infected
  state[initI] <- I_CODE                            # assign initial infections
  
  H_ids <- unique(h)                                # unique household IDs
  HH <- vector("list", length(H_ids))              # initialize household membership list
  names(HH) <- as.character(H_ids)                # name elements by household ID
  for (hid in H_ids) HH[[as.character(hid)]] <- which(h == hid) # store household members
  
  S_daily <- integer(nt)                            # initialize daily susceptible counts
  E_daily <- integer(nt)                            # initialize daily exposed counts
  I_daily <- integer(nt)                            # initialize daily infectious counts
  R_daily <- integer(nt)                            # initialize daily recovered counts
  
  c_mix <- alpha[3] * nc / (bbar^2 * (n - 1))      # random mixing constant
  
  for (tt in tvec) {                                # loop over each day
    indS <- which(state == S_CODE)                 # indices of susceptible individuals
    indE <- which(state == E_CODE)                 # indices of exposed individuals
    indI <- which(state == I_CODE)                 # indices of infectious individuals
    indR <- which(state == R_CODE)                 # indices of recovered individuals
    
    S_daily[tt] <- length(indS)                    # record number of susceptibles
    E_daily[tt] <- length(indE)                    # record number of exposed
    I_daily[tt] <- length(indI)                    # record number of infectious
    R_daily[tt] <- length(indR)                    # record number of recovered
    
    if (length(indE)) {                             # E -> I transition
      uEI <- runif(length(indE))                    # random numbers for transition
      toI <- indE[uEI < gamma]                      # exposed becoming infectious
      if (length(toI)) state[toI] <- I_CODE         # update state to infectious
    }
    
    indI <- which(state == I_CODE)                 # re-identify infectious individuals
    if (length(indI)) {                             # I -> R transition
      uIR <- runif(length(indI))                    # random numbers for recovery
      toR <- indI[uIR < delta]                      # identify recovering individuals
      if (length(toR)) state[toR] <- R_CODE         # update state to recovered
    }
    
    indS <- which(state == S_CODE)                 # re-identify susceptible indices
    indI <- which(state == I_CODE)                 # re-identify infectious indices
    if (length(indS) && length(indI)) {            # proceed if both S and I exist
      
      I_hh <- integer(length(H_ids))               # initialize infected counts per household
      names(I_hh) <- as.character(H_ids)           # name by household ID
      I_tab <- table(h[indI])                       # count infectious per household
      if (length(I_tab)) I_hh[names(I_tab)] <- as.integer(I_tab) # fill household infection counts
      
      I_in_S_hh <- I_hh[as.character(h[indS])]    # infectious in each susceptible's household
      Pnh <- (1 - alpha[1])^I_in_S_hh             # probability of avoiding household infection
      
      I_flag <- logical(n)                          # initialize vector for infectious individuals
      I_flag[indI] <- TRUE                          # mark infectious individuals
      count_I_in_S <- integer(length(indS))        # initialize contact counts per susceptible
      
      for (k in seq_along(indS)) {                 # loop over each susceptible
        j <- indS[k]                                # current susceptible index
        nbrs <- alink[[j]]                          # get neighbors from network
        if (length(nbrs)) count_I_in_S[k] <- sum(I_flag[nbrs]) # count infectious neighbors
      }
      
      Pnc <- (1 - alpha[2])^count_I_in_S            # probability avoiding network infection
      sum_beta_I <- sum(beta[indI])                 # sum of infectiousness
      Pnr <- exp(- c_mix * beta[indS] * sum_beta_I) # probability avoiding random infection
      
      P_not_infected <- Pnh * Pnc * Pnr            # overall probability of avoiding infection
      P_infected <- 1 - P_not_infected             # probability of being infected
      
      uSE <- runif(length(indS))                    # random numbers for S->E
      toE <- indS[uSE < P_infected]                # identify susceptible becoming exposed
      if (length(toE)) state[toE] <- E_CODE        # update state to exposed
    }
  }
  
  list(S = S_daily, E = E_daily, I = I_daily, R = R_daily, t = tvec) # return daily SEIR counts
}


# Plotting of SEIR trajectories
# =======================================

plot_nseir <- function(sim, main = "SEIR with Households & Contacts") {
  stopifnot(all(c("S","E","I","R","t") %in% names(sim)))  # check that sim contains all required components
  mat <- cbind(S = sim$S, E = sim$E, I = sim$I, R = sim$R) # combine SEIR vectors into a single matrix for plotting
  op <- par(mar = c(4.2, 4.5, 3.5, 1.2))                  # set plot margins for readability
  on.exit(par(op))                                        # restore previous plot parameters when function exits
  matplot(sim$t, mat, type = "l", lwd = 2,               # plot all four time series as lines
          xlab = "Day", ylab = "Population count", main = main, # set x and y labels and plot title
          lty = 1)                                        # use solid line type for all series
  legend("right", inset = 0.01, lwd = 2, col = 1:4, lty = 1, # add legend to the right side
         legend = c("S","E","I","R"), bg = "white", cex = 0.52) # specify legend labels, background, and font size
}


# ===============================================
# = Compare 4 scenarios & side-by-side plotting =
# ===============================================
se a single RNG seed at start, do not re-seed inside each simulation.

run_four_scenarios <- function(n = 1000, nt = 150, hmax = 5, nc = 15,
                               alpha_full = c(0.1, 0.01, 0.01),
                               alpha_random_only = c(0, 0, 0.04),
                               delta = 0.2, gamma = 0.4, pinf = 0.005,
                               seed = 1) {
  # Scenarios:
# A) Full model, beta ~ U(0,1)
# B) Random mixing only (alpha_h = alpha_c = 0, alpha_r = 0.04)
# C) Full model, constant beta = mean(beta)
# D) Random mixing + constant beta
# For fair comparison, use a single RNG seed at start, do not re-seed inside each simulation.
  if (!is.null(seed)) set.seed(seed)  # set random seed if provided for reproducibility
  
  betaA <- runif(n, 0, 1)             # generate n beta values uniformly between 0 and 1
  household_sizes <- sample(1:hmax, ceiling(n/mean(1:hmax)), replace = TRUE)
  h <- rep(seq_along(household_sizes), household_sizes)[1:n]
  
  alink <- get.net(betaA, nc = nc, h = h)   # generate contact network using betaA and households
  
  simA <- nseir(betaA, h, alink,        # scenario A: full SEIR model
                alpha = alpha_full, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  simB <- nseir(betaA, h, alink,        # scenario B: random mixing only
                alpha = alpha_random_only, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  betaC <- rep(mean(betaA), n)          # scenario C: constant beta equal to mean(betaA)
  alinkC <- get.net(betaC, nc = nc, h = h) # generate new network using constant beta
  
  simC <- nseir(betaC, h, alinkC,      # scenario C: full SEIR with constant beta
                alpha = alpha_full, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  simD <- nseir(betaC, h, alinkC,      # scenario D: random mixing + constant beta
                alpha = alpha_random_only, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  mk_title <- function(lbl, sim) {
    # Helper function to create a compact title showing peak I and final R
    peakI <- max(sim$I)                 # find maximum number of infectious individuals
    finR  <- tail(sim$R, 1)             # find final number of recovered individuals
    paste0(lbl, "\npeak I = ", peakI, ", final R = ", finR) # compose title string
  }
  
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3.5, 1.2)) # set 2x2 plot layout and margins
  on.exit(par(op), add = TRUE)       # restore original plotting parameters on exit
  plot_nseir(simA, main = mk_title("A) Full model, beta ~ U(0,1)", simA)) # plot scenario A
  plot_nseir(simB, main = mk_title("B) Random mixing only (αh=αc=0, αr=0.04)", simB)) # plot scenario B
  plot_nseir(simC, main = mk_title("C) Full model, constant beta = mean(beta)", simC)) # plot scenario C
  plot_nseir(simD, main = mk_title("D) Random mixing + constant beta", simD)) # plot scenario D
  
  invisible(list(A = simA, B = simB, C = simC, D = simD,  # return results invisibly
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
  nt = 100,
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
# Brief commentary:
# ---------------------------
#  runs, scenarios with structured mixing (A, C) produce
# - slightly later and/or lower peaks than pure random mixing (B, D),
# - smaller final size when heterogeneity in beta is present (A vs C),
# consistent with the notes: variability and clustering suppress spread.
