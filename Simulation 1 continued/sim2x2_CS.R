##############
## FUNCTION ##
##############

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
sd0 <- 1 # random-intercept standard deviation
sd <- 0.5 # equal random-slope standard deviations, sd = sd1 = sd2 = 0.5 * sd12
(rho <- (sd0^2 - .25 * sd^2) / (sd0^2 + .75 * sd^2 + 1)) # correlation under sdr = 1

BF01_AB <- BF_errorAB <- BF01_A <- BF_errorA <- BF01_B <- BF_errorB <- 
  legacy01_AB <- legacy_errorAB <- legacy01_A <- legacy_errorA <- legacy01_B <- legacy_errorB <- 
  pValAB <- pValA <- pValB <- FValAB <- FValA <- FValB <- numeric(NSim)

library(BayesFactor) # v 0.9.12-4.7

for (i in 1:NSim) {
  if (i %% 10 == 0) print (paste0("Run ", i, ": ", Sys.time()))
  
  dat <- simData(N=N, n=R, B1=ES, sdr=sdr, sd0=sd0, sd1=sd, sd2=sd, sd12=2*sd, seed=i)
  
  RMANOVA <- aov(DV ~ (A * B) + Error(ID / (A * B)), dat) # repeated-measures analysis of variance
  sum_stat <- unlist(summary(RMANOVA))
  pValAB[i] <- unname(sum_stat["Error: ID:A:B.Pr(>F)1"])
  pValA[i] <- unname(sum_stat["Error: ID:A.Pr(>F)1"])
  pValB[i] <- unname(sum_stat["Error: ID:B.Pr(>F)1"])
  FValAB[i] <- unname(sum_stat["Error: ID:A:B.F value1"])
  FValA[i] <- unname(sum_stat["Error: ID:A.F value1"])
  FValB[i] <- unname(sum_stat["Error: ID:B.F value1"])
  
  set.seed(277)
  legacy_full <- lmBF(DV ~ A * B + ID, whichRandom="ID", data=dat, iterations=1e5, progress=F)
  
  set.seed(277)
  legacy_add <- lmBF(DV ~ A + B + ID, whichRandom="ID", data=dat, iterations=1e5, progress=F)
  
  set.seed(277)
  legacy_omitA <- lmBF(DV ~ B + ID, whichRandom="ID", data=dat, iterations=1e5, progress=F)
  
  set.seed(277)
  legacy_omitB <- lmBF(DV ~ A + ID, whichRandom="ID", data=dat, iterations=1e5, progress=F)
  
  bf <- extractBF(legacy_add / legacy_full)
  legacy01_AB[i] <- bf$bf
  legacy_errorAB[i] <- bf$error
  
  bf <- extractBF(legacy_omitA / legacy_add)
  legacy01_A[i] <- bf$bf
  legacy_errorA[i] <- bf$error
  
  bf <- extractBF(legacy_omitB / legacy_add)
  legacy01_B[i] <- bf$bf
  legacy_errorB[i] <- bf$error
  
  
  if(R > 1 & rho > -.3) {
    set.seed(277)
    BF_full <- lmBF(DV ~ A * B + ID + A:ID + B:ID + A:B:ID, dat, whichRandom="ID", iterations=1e5, progress=F) # MRE (unaggregated)
    
    set.seed(277)
    BF_add <- lmBF(DV ~ A + B + ID + A:ID + B:ID + A:B:ID, dat, whichRandom="ID", iterations=1e5, progress=F)
    
    set.seed(277)
    BF_omitA <- lmBF(DV ~ B + ID + A:ID + B:ID + A:B:ID, dat, whichRandom="ID", iterations=1e5, progress=F)
    
    set.seed(277)
    BF_omitB <- lmBF(DV ~ A + ID + A:ID + B:ID + A:B:ID, dat, whichRandom="ID", iterations=1e5, progress=F)
    
  } else {
    set.seed(277)
    BF_full <- lmBF(DV ~ A * B + ID + A:ID + B:ID, dat, whichRandom="ID", iterations=1e5, progress=F) # MRE (aggregated)
    
    set.seed(277)
    BF_add <- lmBF(DV ~ A + B + ID + A:ID + B:ID, dat, whichRandom="ID", iterations=1e5, progress=F)
    
    set.seed(277)
    BF_omitA <- lmBF(DV ~ B + ID + A:ID + B:ID, dat, whichRandom="ID", iterations=1e5, progress=F)
    
    set.seed(277)
    BF_omitB <- lmBF(DV ~ A + ID + A:ID + B:ID, dat, whichRandom="ID", iterations=1e5, progress=F)
  }
  
  bf <- extractBF(BF_add / BF_full)
  BF01_AB[i] <- bf$bf
  BF_errorAB[i] <- bf$error
  
  bf <- extractBF(BF_omitA / BF_add)
  BF01_A[i] <- bf$bf
  BF_errorA[i] <- bf$error
  
  bf <- extractBF(BF_omitB / BF_add)
  BF01_B[i] <- bf$bf
  BF_errorB[i] <- bf$error
} # roughly 4.5 / 29.3 / 8 hours of run time for N = 48 / 100 and n = 1 / 3

report <- data.frame(BF01_AB, BF_errorAB, BF01_A, BF_errorA, BF01_B, BF_errorB, 
                     legacy01_AB, legacy_errorAB, legacy01_A, legacy_errorA, legacy01_B, legacy_errorB, 
                     pValAB, pValA, pValB, FValAB, FValA, FValB)
write.csv(report, paste0("~/Downloads/osfstorage-archive/Simulation 1 continued/data/",
                         ifelse(sdr == 1, "CS", paste0("Het", round(sdr,2))), "_N",
                         N, ifelse(rho > -.3 | R==1, "_n", "_s"), R, "_ES", ES,
                         ifelse(sdr == 1, paste0("_rho", round(rho,2)), ""), ".csv"), row.names=F)


###################
## VISUALIZATION ##
###################

library(dplyr) # v 1.1.4
library(ggplot2) # v 3.5.2

# sdr	sd0 	 sd 	 rho	R	      csv
# 1	  0.5 	0.5 	0.13	3	      CS_N48_n3_ES0.5_rho0.13
# 1	    1 	0.5 	0.43	3	      CS_N48_n3_ES0.5_rho0.43
# 1	  0.5 	1.5  -0.11	3	      CS_N48_n3_ES0.5_rho-0.11
# 1	  0.5 	  5 	-0.3	3	      CS_N48_s3_ES0.5_rho-0.3
# 					                    CS_N48_n3_ES0.5_rho-0.3
# 					                    CS_N48_s3_ES0_rho-0.3
# 3	  0.5 	0.5 		    1	      Het3_N48_n1_ES0.5
# 2	  0.5 	0.5		      1	      Het2_N48_n1_ES0.5
# 1 	0.5 	  5	  -0.3	1	      CS_N48_n1_ES0.5_rho-0.3
# 					                    CS_N48_n1_ES0_rho-0.3
# 1	  0.5 	1.5	 -0.11	1	      CS_N48_n1_ES0.5_rho-0.11
# 1 	0.5 	0.5	  0.13	1	      CS_N48_n1_ES0.5_rho0.13
# 					                    CS_N48_n1_ES0.25_rho0.13
# 					                    CS_N48_n1_ES0_rho0.13
# 					                    CS_N100_n1_ES0.5_rho0.13
# 					                    CS_N100_n1_ES0_rho0.13

N <- 48 # number of participants
R <- 3 # number of (balanced) trials
ES <- 0.5 # fixed effect size for factor A of interest (this is the SECOND factor in our paper)
sdr <- 1 # whether residuals are homoscedastic, i.e., sdr = 1
sd0 <- 0.5 # random-intercept standard deviation
sd <- 0.5 # equal random-slope standard deviations, sd = sd1 = sd2 = 0.5 * sd12
(rho <- (sd0^2 - .25 * sd^2) / (sd0^2 + .75 * sd^2 + 1)) # correlation under sdr = 1
df <- read.csv(paste0("~/Downloads/osfstorage-archive/Simulation 1 continued/data/",
                      ifelse(sdr == 1, "CS", paste0("Het", round(sdr,2))), "_N",
                      N, ifelse(rho > -.3 | R==1, "_n", "_s"), R, "_ES", ES,
                      ifelse(sdr == 1, paste0("_rho", round(rho,2)), ""), ".csv"))

lapply(df[, grepl("error", names(df))], quantile, na.rm=T) # check the maximum Monte Carlo error rate

if (sdr == 1) {
  title <- bquote(italic("n") * "=" * .(N) * ", " * italic("R") * "=" * .(R) * ", CS Holds*, " * italic("ρ") * "=" * .(round(rho,2)))
} else {
  title <- bquote(italic("n") * "=" * .(N) * ", " * italic("R") * "=" * .(R) * ", CS is Violated*")
}

fig <- ggplot(filter(df, !if_any(contains("error"), ~ . > .05)),
              aes(-log(legacy01_A), -log(BF01_A), # ln BF₁₀
                  color=factor(cut(pValA, breaks=c(-1E-5, .01, .05, .1, 1),
                                   labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
                  shape=factor(cut(pValA, breaks=c(-1E-5, .01, .05, .1, 1),
                                   labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.5) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("RIO "~italic("BF")[10]),
       y=expression("MRE "~italic("BF")[10]),
       color=expression(paste(italic("p-"), "value")),
       shape=expression(paste(italic("p-"), "value")),
       title=title) +
  scale_x_continuous(limits=range(-log(c(df$legacy01_A, df$BF01_A))), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3), log(100)),
                     labels=c("1/3  ", 1, 3, expression(10^2))) + # white space left intentionally
  scale_y_continuous(limits=range(-log(c(df$legacy01_A, df$BF01_A))), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3), log(100)),
                     labels=c("1/3", 1, 3, expression(10^2))) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  coord_fixed(ratio=1) +
  theme_classic() +
  theme(text=element_text(size=16),
        plot.margin=margin(5, 5, 5, 5, "mm"),
        legend.title=element_text(margin=margin(0,0.3,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.3,0,0.2,"cm")),
        legend.position="bottom", legend.direction="horizontal") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))

ggsave(paste0("~/Downloads/osfstorage-archive/Simulation 1 continued/supp/", 
              ifelse(sdr == 1, "CS", paste0("Het", round(sdr,2))), "_N", N, 
              ifelse(rho > -.3 | R==1, "_n", "_s"), R, "_ES", ES, 
              ifelse(sdr == 1, paste0("_rho", round(rho,2)), ""), ".png"), 
       fig, width=2250, height=2500, units="px", device="png", dpi=300)

