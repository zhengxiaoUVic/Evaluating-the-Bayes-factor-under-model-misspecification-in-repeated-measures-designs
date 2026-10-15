## Reviewer 1: The authors point out that the median Monte Carlo error exceeds 48%. 
## I see in the code on the osf that the authors used 10^5 iterations to reduce the Monte Carlo error. 
## I’m wondering if increasing this further to 10^6 could help reduce the Monte Carlo error further?


##############
## FUNCTION ## SAME
############## ~/Downloads/osfstorage-archive/Simulation 1 continued/sim2x2_CS.R

simData <- function(N, n, B1, sdr=1, sd0=0.5, isRIO=F, sd1=0.5, sd2=0.5, sd12=1, seed) {
  #' Input -
  #' N:         number of participants
  #' n:         number of (balanced) trials
  #' B1:        fixed effect size for factor A of interest
  #' sdr:       a positive value specifying that 
  #'            residuals have sd set to 1 and `sdr` for each level of factor B
  #' sd0:       standard deviation for random intercepts
  #' isRIO:     if TRUE, random slopes have sd suppress to 0 (i.e. random intercepts only)
  #' sd1:       standard deviation for random slopes for (contrasted) factor A
  #' sd2:       standard deviation for random slopes for (contrasted) factor B
  #' sd12:      standard deviation for random slopes for the (contrasted) interaction AB
  #' seed:      random seed
  #'
  #' Output - a long format data frame for a 2x2 repeated-measures design
  
  C1 <- C2 <- c(-.5, .5) # contrast vectors for factors A and B, respectively (sum coding)
  B2 <- 0.5 # fixed effect size for factor B
  
  set.seed(seed)
  s <- rnorm(N, 0, sd0) # random intercepts for each participant
  
  # fixed or random slopes for factor A (this is the SECOND factor in our paper)
  b1 <- rnorm(N, B1, ifelse(isRIO, 0, sd1))
  
  # fixed or random slopes for factor B (this is the FIRST factor in our paper)
  b2 <- rnorm(N, B2, ifelse(isRIO, 0, sd2))
  
  # fixed or random slopes for the interaction AB
  b12 <- rnorm(N, 0, ifelse(isRIO, 0, sd12))
  
  # a matrix of condition contrasts for the 2x2 design
  conds <- expand.grid("A"=C1, "B"=C2)
  conds$AB <- conds$A * conds$B
  
  # each row of `M` represents one participant's set of four condition means
  M <- s + cbind(b1, b2, b12) %*% t(conds) # a1b1, a2b1, a1b2, a2b2 (N×4 matrix)
  
  het <- rep(c(1, sdr), each=N*2) # residual standard deviations (homoscedasticity when sdr = 1)
  DV <- rnorm(N*4*n, c(M), het) # adding noise
  
  data.frame("DV" = DV,
             "ID" = rep(paste0("s",1:N), 4*n),
             "A" = rep(c("a1","a2"), each=N, times=2*n),
             "B" = rep(c("b1","b2"), each=N*2, times=n),
             "Trial" = rep(paste0("r",1:n), each=N*4),
             stringsAsFactors = T) # long format data
}


################
## SIMULATION ##
################

N <- 48 # number of participants
R <- 3 # number of (balanced) trials
ES <- 0.5 # (medium) fixed effect size for factor A of interest
NSim <- 200 # number of simulation runs for each setting
sdr <- 1 # whether residuals are homoscedastic, i.e., sdr = 1
sd0 <- 0.5 # random-intercept standard deviation
sd <- 5 # equal random-slope standard deviations, sd = sd1 = sd2 = 0.5 * sd12 (this sd=5 is σ=2.5 in our paper—Table 5)
(rho <- (sd0^2 - .25 * sd^2) / (sd0^2 + .75 * sd^2 + 1)) # correlation under sdr = 1


dat <- simData(N=N, n=R, B1=ES, sdr=sdr, sd0=sd0, sd1=sd, sd2=sd, sd12=2*sd, seed=2)


library(BayesFactor) # v 0.9.12-4.7

{
  iter <- 1e6 # roughly 22 minutes of run time
  
  set.seed(277)
  BF_full <- lmBF(DV ~ A * B + ID + A:ID + B:ID + A:B:ID, dat, whichRandom="ID", iterations=iter, progress=F) # MRE (unaggregated)
  
  set.seed(277)
  BF_add <- lmBF(DV ~ A + B + ID + A:ID + B:ID + A:B:ID, dat, whichRandom="ID", iterations=iter, progress=F)
  
  set.seed(277)
  BF_omitA <- lmBF(DV ~ B + ID + A:ID + B:ID + A:B:ID, dat, whichRandom="ID", iterations=iter, progress=F)
  
  set.seed(277)
  BF_omitB <- lmBF(DV ~ A + ID + A:ID + B:ID + A:B:ID, dat, whichRandom="ID", iterations=iter, progress=F)
}

## Results: the median Monte Carlo error is reduced to ~20%.

