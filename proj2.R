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
##   3. Random background infection occurs with probability αr * nc * (β_i * β_j) /
##      (mean(β)^2 * (n - 1)), where βi represents 'sociability' (infectiousness) parameter for the ith person
##      and nc is the average number of contacts per person.
## E → I : with daily probability γ (infection incubation rate)
## I → R : with daily probability δ (recovery rate)
## The aim is to investigate how incorporating social structure changes 
## epidemic dynamic compared to purely random mixing.

hmax <- 5  ## use household size maximum to be 5 by default 
n <- 1000  ## our code work with any population sizes, here we test and develop with n=1,000

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
## beta is the n vector of β_i value for each person
## nc is the average number of contacts per person
  n <- length(beta)  ## total number of individuals
  if (n < 2L) return(vector("list", n))  ## return empty list if less than 2 people
  beta_bar <- mean(beta)  ## mean sociability (infectivity) parameter
  cst <- nc / (beta_bar^2 * (n - 1))  ## constant factor for exact pairwise probability: Pr_ij = cst * β_i * β_j.Note: (n-1) because each person can link to n-1 others (excluding themselves)
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
    js <- setdiff(js, hh_members)  ## exclude same household members 
    if (length(js) == 0) next  ## skip if no partners left
    
    p <- cst * beta[i] * beta[js] ## compute link probabilities
    p[p < 0] <- 0  ## ensure probabilities >= 0
    p[p > 1] <- 1  ## ensure probabilities <= 1
    
    u <- runif(length(js))  ## draw uniform random numbers
    keep <- which(u < p)  ## keep edges where u < probability
    
    if (length(keep)) {  ## if any edges are kept
      nbrs <- js[keep]   ## get connected neighbours indices for i
      alink[[i]] <- c(alink[[i]], nbrs)  ## add neighbours to i's adjacency list
      for (j in nbrs) alink[[j]] <- c(alink[[j]], i) ## symmetrically add i to each neighbor's adjacency list
    }
  }
  
  for (i in seq_len(n)) {  ## loop to clean adjacency list
    if (length(alink[[i]]) > 1)  ## only if person has multiple neighbors
      alink[[i]] <- sort(unique(alink[[i]]))  ## remove duplicates and sort
  }
  
  return(alink) ## return adjacency list
}

nseir <- function(beta, h, alink, ## infection rates, household memberships, and regular contacts
                 alpha = c(0.1, 0.01, 0.01), ## infection (household, network, random)
                 delta = 0.2, ## daily probability of recovery
                 gamma = 0.4, ## daily probability of becoming infectious, incubation probability
                 nc = 15, ## average number of contacts per person
                 nt = 100, ## number of days
                 pinf = 0.005, ## proportion of the initial population to randomly start in the I state 
                 seed = NULL, ## optional RNG seed for reproducibility
                 exact_random = FALSE ## if TRUE compute exact random-product (slower)
) {
## This function simulates an SEIR epidemic model with household, regular contact network and random mixing. 
## It tracks transitions between S, E, I, and R states over nt days for a population with given beta, h, and alink.
## The model includes infection spread through household, regular network, and random contacts 
## with respective strengths alpha, and accounts for daily infection, exposure, and 
## recovery probabilities (gamma, delta). It returns a list containing daily totals 
## of S, E, I, R, and the time vector t.
  if (!is.null(seed)) set.seed(seed) ## set RNG seed if provided
  
  n <- length(beta) ## total number of individuals            
  if (!is.numeric(beta) || n < 1) stop("beta must be a non-empty numeric vector")
  if (length(h) != n) stop("h must be same length as beta")
  if (!is.list(alink) || length(alink) != n) stop("alink must be a list of length n")
  if (!is.numeric(alpha) || length(alpha) < 3) stop("alpha must be numeric length >= 3")
  
  S_CODE <- 1L; E_CODE <- 2L; I_CODE <- 3L; R_CODE <- 4L ## codes for states
  
  state <- rep.int(S_CODE, n) ## initial state: all susceptible
  nI0 <- max(1L, round(pinf * n)) ## compute initial infectious count
  initI <- sample.int(n, size = nI0, replace = FALSE) 
  state[initI] <- I_CODE ## set initial infections      
  
  H_ids <- unique(h) ## unique household IDs
  HH <- vector("list", length(H_ids)) ## initialize household membership list
  names(HH) <- as.character(H_ids) ## name elements by household ID
  for (hid in H_ids) HH[[as.character(hid)]] <- which(h == hid) ## store household members
  
  S_daily <- integer(nt) ## initialize daily susceptible counts
  E_daily <- integer(nt) ## initialize daily exposed counts
  I_daily <- integer(nt) ## initialize daily infectious counts
  R_daily <- integer(nt) ## initialize daily recovered counts
  tvec <- seq_len(nt) ## vector of time steps
  
  beta_bar <- mean(beta) ## mean sociability
  if (beta_bar <= 0) stop("mean(beta) must be positive")
  constant_mix <- alpha[3] * nc / (beta_bar^2 * (n - 1)) ## constant for random mixing
  
for (tt in tvec) {
  # record current counts (counts at day start)
  S_daily[tt] <- sum(state == S_CODE)  ## count susceptible individuals
  E_daily[tt] <- sum(state == E_CODE)  ## count exposed individuals
  I_daily[tt] <- sum(state == I_CODE)  ## count infectious individuals
  R_daily[tt] <- sum(state == R_CODE)  ## count recovered individuals
  
  
  indS <- which(state == S_CODE) ## indices of susceptible individuals
  indE <- which(state == E_CODE) ## indices of exposed individuals
  indI <- which(state == I_CODE) ## indices of infectious individuals
  
  ## Compute S -> E based on day-start infected individuals
  if (length(indS) > 0 && length(indI) > 0) {  ## only if there are both S and I individuals
    I_tab <- table(h[indI])  ## count number of infectious individuals per household
    I_in_S_hh <- as.integer(I_tab[as.character(h[indS])])  ## map infected count to susceptibles’ households
    I_in_S_hh[is.na(I_in_S_hh)] <- 0  ## replace NA with 0 if no infected in that household
    P_avoid_hh <- (1 - alpha[1]) ^ I_in_S_hh  ## probability of avoiding infection from household infected individuals
    
    ## Regular Contact Network
    inf_flag <- logical(n); inf_flag[indI] <- TRUE  ## logical flag vector marking currently infectious individuals
    
    ## extract each susceptible’s list of network contacts
    network_lists <- alink[indS]  ## alink provides adjacency lists for all individuals
    
    ## count the number of infectious contacts for each susceptible individual
    ## vapply is used instead of sapply to enforce integer output and ensure performance
    if (length(network_lists) > 0L) {
      count_I_network <- vapply(network_lists, function(net) {  ## iterate over each susceptible’s contact list
        if (length(net) == 0L) return(0L)  ## if the person has no network contacts, count = 0
        sum(inf_flag[net])  ## count how many of their contacts are infectious
      }, integer(1))
    } else {
      count_I_network <- integer(0)  ## if there are no susceptible individuals, create empty vector
    }
    
    P_avoid_network <- (1 - alpha[2]) ^ count_I_network  ## each infectious contact contributes multiplicatively to infection risk
    
    sum_beta_I <- sum(beta[indI]) ## total infectivity across all infectious individuals
    if (sum_beta_I == 0) {
      P_avoid_rand <- rep(1, length(indS)) ## no infectiousness => avoid prob = 1
    } else {
      if (!exact_random) {
        # approximation: avoid probability = exp(- constant_mix * beta_j * sum_beta_I)
        P_avoid_rand <- exp(- constant_mix * beta[indS] * sum_beta_I)  ## compute exponential approximation
        P_avoid_rand[P_avoid_rand < 0] <- 0  ## ensure values are not below 0
        P_avoid_rand[P_avoid_rand > 1] <- 1  ## ensure values are not above 1
      } else {
        ## exact calculation: product over all infectious individuals
        P_avoid_rand <- numeric(length(indS)) ## initialize vector
        for (k in seq_along(indS)) {  ## loop through susceptibles
          j <- indS[k]  ## index of the susceptible
          pij <- constant_mix * beta[indI] * beta[j]  ## compute pairwise infection probabilities
          pij[pij > 1] <- 1  ## cap probabilities at 1
          if (any(pij >= 1)) {
            P_avoid_rand[k] <- 0  ## if any p_ij=1, infection is certain
          } else {
            P_avoid_rand[k] <- exp(sum(log1p(-pij)))  ## compute product of (1 - p_ij)
          }
        }
      }
    }
    
    P_avoid_all <- P_avoid_hh * P_avoid_network * P_avoid_rand  ## combine independent avoidance probabilities
    P_infect <- 1 - P_avoid_all  ## overall infection probability for each susceptible
    
    # perform Bernoulli draws for S->E
    draws_SE <- runif(length(indS)) < P_infect ## random draws for S→E transitions
    newE <- indS[draws_SE] ## susceptibles becoming exposed today
  } else {
    newE <- integer(0) ## no new exposures if no S or no I
  }
  
  ## Compute E -> I (progression)
  if (length(indE) > 0) { ## if there are any exposed individuals today
    draws_EI <- runif(length(indE)) < gamma ## draw a random number to determine if they progress to infectious state
    newI_fromE <- indE[draws_EI] ## select indices of exposed individuals who become infectious today
  } else {
    newI_fromE <- integer(0) ## if no exposed individuals exist, no new infectious individuals will appear
  }
  
  ## Compute I -> R (recoveries)
  if (length(indI) > 0) { ## if there are any infectious individuals today
    draws_IR <- runif(length(indI)) < delta ## draw a random number to determine if they recover today
    newR <- indI[draws_IR] ## select indices of infectious individuals who recover today
  } else {
    newR <- integer(0) ## if no infectious individuals exist, no new recoveries will occur
  }
  
  ## Updates (based on day-start indices)
  ## Please note that updates are based on day-start sets so newly created categories do not act immediately
  if (length(newE) > 0) state[newE] <- E_CODE ## update new exposures
  if (length(newI_fromE) > 0) state[newI_fromE] <- I_CODE ## update new infectious individuals
  if (length(newR) > 0) state[newR] <- R_CODE ## update new recoveries
  ## move to next day
} 

return(list(S = S_daily, E = E_daily, I = I_daily, R = R_daily, t = tvec)) ## return daily time series (counts at day starts)
}

plot_nseir <- function(sim, main = "SEIR with Households & Contacts") {
## combines SEIR into a matrix and plot their trajectories over time.
## use different colors to represent each state and add a legend to distinguish the four lines.
  stopifnot(all(c("S","E","I","R","t") %in% names(sim))) ## check that sim contains all required components
  mat <- cbind(S = sim$S, E = sim$E, I = sim$I, R = sim$R) ## combine SEIR vectors into a single matrix for plotting
  op <- par(mar = c(4.2, 4.5, 3.5, 1.2)) ## set plot margins for readability
  on.exit(par(op)) ## restore previous plot parameters when function exits
  matplot(sim$t, mat, type = "l", lwd = 2, ## plot all four time series as lines
          xlab = "Day", ylab = "Population count", main = main, ## set x and y labels and plot title
          lty = 1) ## use solid line type for all series
  legend("right", inset = 0.01, lwd = 2, col = 1:4, lty = 1, ## add legends to the right side
         legend = c("S","E","I","R"), bg = "white", cex = 0.52) ## specify legend labels, background, and font size
}

run_four_scenarios <- function(n = 1000, nt = 150, hmax = 5, nc = 15,
                               alpha_full = c(0.1, 0.01, 0.01),
                               alpha_random_only = c(0, 0, 0.04),
                               delta = 0.2, gamma = 0.4, pinf = 0.005,
                               seed = 1) {
## Scenarios:
## A) Full model, beta ~ U(0,1)
## B) Random mixing only (alpha_h = alpha_c = 0, alpha_r = 0.04)
## C) Full model, constant beta = mean(beta)
## D) Random mixing + constant beta
## For fair comparison, use a single RNG seed at start, do not re-seed inside each simulation.
  if (!is.null(seed)) set.seed(seed)  ## set random seed if provided for reproducibility
  
  betaA <- runif(n, 0, 1)             ## generate n beta values uniformly between 0 and 1
  household_sizes <- sample(1:hmax, ceiling(n/mean(1:hmax)), replace = TRUE)
  h <- rep(seq_along(household_sizes), household_sizes)[1:n]
  
  alink <- get.net(betaA, nc = nc, h = h)   ## generate contact network using betaA and households
  
  simA <- nseir(betaA, h, alink,        ## scenario A: full SEIR model
                alpha = alpha_full, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  simB <- nseir(betaA, h, alink,        ## scenario B: random mixing only
                alpha = alpha_random_only, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  betaC <- rep(mean(betaA), n)          ## scenario C: constant beta equal to mean(betaA)
  alinkC <- get.net(betaC, nc = nc, h = h) # generate new network using constant beta
  
  simC <- nseir(betaC, h, alinkC,      ## scenario C: full SEIR with constant beta
                alpha = alpha_full, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  simD <- nseir(betaC, h, alinkC,      ## scenario D: random mixing + constant beta
                alpha = alpha_random_only, delta = delta, gamma = gamma,
                nc = nc, nt = nt, pinf = pinf)
  
  mk_title <- function(lbl, sim) {
    # Helper function to create a compact title showing peak I and final R
    peakI <- max(sim$I)                 ## find maximum number of infectious individuals
    finR  <- tail(sim$R, 1)             ## find final number of recovered individuals
    paste0(lbl, "\npeak I = ", peakI, ", final R = ", finR) ## compose title string
  }
  
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3.5, 1.2)) ## set 2x2 plot layout and margins
  on.exit(par(op), add = TRUE)       ## restore original plotting parameters on exit
  plot_nseir(simA, main = mk_title("A) Full model, beta ~ U(0,1)", simA)) ## plot scenario A
  plot_nseir(simB, main = mk_title("B) Random mixing only (αh=αc=0, αr=0.04)", simB)) ## plot scenario B
  plot_nseir(simC, main = mk_title("C) Full model, constant beta = mean(beta)", simC)) ## plot scenario C
  plot_nseir(simD, main = mk_title("D) Random mixing + constant beta", simD)) ## plot scenario D
  
  invisible(list(A = simA, B = simB, C = simC, D = simD,  ## return results invisibly
                 betaA = betaA, betaC = betaC, h = h,
                 alinkA = alink, alinkC = alinkC))
}

## Example: run the scenarios
## (Reproducible with a single seed at the start; no reseeding inside.)

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
  seed = NULL  ## already seeded above; leave NULL to avoid resetting inside
)

## Brief commentary:
## By comparing Figure A and Figure B, when the household and regular contact network structure are removed, 
##    random mixing leads to a faster and wider spread of infection, resulting in a higher short-term epidemic peak.
## By comparing Figure A and Figure C, smaller final size when heterogeneity in beta is present,
##    consistent with the notes -- variability and clustering suppress spread.
