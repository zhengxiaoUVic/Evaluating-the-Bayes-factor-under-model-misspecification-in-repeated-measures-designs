##############
## FUNCTION ##
##############

simData <- function(N=48, eqVar=T, eqCov=F, null=F, seed=277) {
  #' Input -
  #' N:     number of participants
  #' eqVar: whether the variances are equal.
  #' eqCov: whether the covariances are equal across the four conditions (a1b1, a2b1, a1b2, a2b2).
  #' null:  whether simulated under the null for factor B of interest.
  #' seed:  random seed
  #'
  #' Dependency - MASS(v 7.3-61)
  #' Output - the long format data frame for a 2x2 repeated-measures design
  
  set.seed(seed)
  if(eqVar & eqCov) { # compound symmetry holds
    if(null) {
      eff <- c(0, 15, 0, 15) # a1b1, a2b1, a1b2, a2b2
    } else {
      eff <- c(0, 15, 4, 19)
    }
    mat <- MASS::mvrnorm(N, rep(100,4)+eff, diag(c(100,100,100,100)))
    mat <- sweep(mat, 1, rnorm(N, sd=12.25), "+") # correlation 12.25^2/(12.25^2+100)≈.6
    
  } else if(eqVar & !eqCov) { # compound symmetry is violated
    x <- rnorm(N, 100, 6)     # set origin scores
    x1 <- rnorm(N, 0, 10)     # error variability
    y1 <- x+x1                # first set of base scores
    x2 <- rnorm(N, 0, 10)     # error variability
    y2 <- x+x2                # second set of base scores
    
    a1 <- rnorm(N, 0, 2.82)   # error variability
    a1b1 <- round(y1+a1)      # first condition
    
    a1 <- rnorm(N, ifelse(null, 0, 1.3), 2.82) # error variability and effect for factor A
    a1b2 <- round(y1+a1)      # third condition
    
    a2 <- rnorm(N, 15, 2.82)   # error variability and effect for factor B
    a2b1 <- round(y2+a2)       # second condition
    
    a2 <- rnorm(N, ifelse(null, 15, 16.3), 2.82) # error variability and effect for factors A and B
    a2b2 <- round(y2+a2)       # fourth condition
    mat <- cbind(a1b1, a2b1, a1b2, a2b2)
    
  } else { # !eqVar & !eqCov (not consider !eqVar & eqCov)
    x <- rnorm(N, 100, 4.24)  # set origin scores
    x1 <- rnorm(N, 0, 1.41)   # error variability
    a1b1 <- round(x+x1)       # first set of a1 scores
    x2 <- rnorm(N, ifelse(null, 0, 1.8), 1.41) # error variability
    a1b2 <- round(x+x2)       # second set of a1 scores
    
    # base scores for a2
    x1 <- rnorm(N, 0, 9.38)   # base error variability
    a2b1 <- x+x1
    a2b2 <- x+x1
    x2 <- rnorm(N, 16, 3.46)    # extra variability and effect
    a2b1 <- round(a2b1+x2)
    x2 <- rnorm(N, ifelse(null, 16, 17.8), 3.46)  # extra variability and effect
    a2b2 <- round(a2b2+x2)
    mat <- cbind(a1b1, a2b1, a1b2, a2b2)
  }
  
  data.frame("ID"=paste0("s", 1:N),
             "DV"=c(mat),
             "A"=rep(c("a1", "a2"), times=2, each=N),
             "B"=rep(c("b1", "b2"), each=N*2),
             stringsAsFactors=T)
}

#################
## SAMPLE DATA ##
#################

#        BF01_B       pValB
# 10  0.3913389 0.008890373
# 14  0.3800049 0.007031997
# 45  0.3773166 0.007365020
# 78  0.3486783 0.005256488
# 80  0.3426748 0.007280270
# 123 0.3837344 0.008734395
# 130 0.3432300 0.005368126
# 194 0.3908563 0.005488747       # report_vTcF_H1, N=48

df <- simData(seed=44661, N=100000)
mat <- matrix(df$DV, ncol=4, dimnames=list(unique(df$ID), c("a1b1", "a2b1", "a1b2", "a2b2")))
round(cov(mat)) # sample covariance matrix
#      a1b1 a2b1 a1b2 a2b2
# a1b1  145   36  137   36
# a2b1   36  144   36  136
# a1b2  137   36  145   36
# a2b2   36  136   36  144

round(cor(mat),2) # sample correlation matrix
#      a1b1 a2b1 a1b2 a2b2
# a1b1 1.00 0.25 0.94 0.25
# a2b1 0.25 1.00 0.25 0.94
# a1b2 0.94 0.25 1.00 0.25
# a2b2 0.25 0.94 0.25 1.00

df <- simData(seed=44661, N=100)
summary(aov(DV ~ (A * B) + Error(ID / (A * B)), df)) # repeated-measures analysis of variance
# Error: ID
#           Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 99  44154     446               
# 
# Error: ID:A
#           Df Sum Sq Mean Sq F value   Pr(>F)    
# A          1  17109   17109   97.03 2.32e-16 ***
# Residuals 99  17455     176                     
# 
# Error: ID:B
#           Df Sum Sq Mean Sq F value   Pr(>F)    
# B          1  204.5  204.49   30.33 2.88e-07 ***
# Residuals 99  667.5    6.74                     
# 
# Error: ID:A:B
#           Df Sum Sq Mean Sq F value Pr(>F)
# A:B        1    6.2   6.250   0.695  0.406
# Residuals 99  889.8   8.987 

library(lme4) # v 1.1-36

fit <- lmer(DV ~ A * B + (1 | ID), data=df)
vc <- VarCorr(fit)
subject_var <- vc$ID[1,1]                  # variance of random intercept (subject-level)
residual_var <- attr(vc, "sc")^2           # residual variance (within-subject)
subject_var / (subject_var + residual_var) # correlation estimate


################
## SIMULATION ##
################

N <- 48 # number of participants
NSim <- 200 # number of simulation runs for each setting

if (T) { # load pre-run data
  report_vTcF_H1 <- read.csv("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N48_eqVarTRUE_eqCovFALSE_nullFALSE.csv")
  report_vTcF_H0 <- read.csv("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N48_eqVarTRUE_eqCovFALSE_nullTRUE.csv")
  report_vTcT_H1 <- read.csv("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N48_eqVarTRUE_eqCovTRUE_nullFALSE.csv")
  report_vTcT_H0 <- read.csv("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N48_eqVarTRUE_eqCovTRUE_nullTRUE.csv")
  report_vFcF_H1 <- read.csv("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N48_eqVarFALSE_eqCovFALSE_nullFALSE.csv")
  report_vFcF_H0 <- read.csv("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N48_eqVarFALSE_eqCovFALSE_nullTRUE.csv")
} else { # or start from scratch
  
  simNull <- F # whether factor B has no effect
  eqVar <- T   # whether variances are equal
  eqCov <- F   # whether covariances are equal
  BF01_AB <- BF_errorAB <- BF01_A <- BF_errorA <- BF01_B <- BF_errorB <- 
    legacy01_AB <- legacy_errorAB <- legacy01_A <- legacy_errorA <- legacy01_B <- legacy_errorB <- 
    pValAB <- pValA <- pValB <- FValAB <- FValA <- FValB <- numeric(NSim)
  
  library(BayesFactor) # v 0.9.12-4.7
  
  for (i in 1:NSim) {
    if (i %% 10 == 0) print (paste0("Run ", i, ": ", Sys.time()))
    
    dat <- simData(N=N, eqVar=eqVar, eqCov=eqCov, null=simNull, seed=i)
    
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
    
    
    set.seed(277)
    BF_full <- lmBF(DV ~ A * B + ID + A:ID + B:ID, whichRandom="ID", data=dat, iterations=1e5, progress=F)
    
    set.seed(277)
    BF_add <- lmBF(DV ~ A + B + ID + A:ID + B:ID, whichRandom="ID", data=dat, iterations=1e5, progress=F)
    
    set.seed(277)
    BF_omitA <- lmBF(DV ~ B + ID + A:ID + B:ID, whichRandom="ID", data=dat, iterations=1e5, progress=F)
    
    set.seed(277)
    BF_omitB <- lmBF(DV ~ A + ID + A:ID + B:ID, whichRandom="ID", data=dat, iterations=1e5, progress=F)
    
    bf <- extractBF(BF_add / BF_full)
    BF01_AB[i] <- bf$bf
    BF_errorAB[i] <- bf$error
    
    bf <- extractBF(BF_omitA / BF_add)
    BF01_A[i] <- bf$bf
    BF_errorA[i] <- bf$error
    
    bf <- extractBF(BF_omitB / BF_add)
    BF01_B[i] <- bf$bf
    BF_errorB[i] <- bf$error
  } # roughly 5 hours of run time
  
  report <- data.frame(BF01_AB, BF_errorAB, BF01_A, BF_errorA, BF01_B, BF_errorB, 
                       legacy01_AB, legacy_errorAB, legacy01_A, legacy_errorA, legacy01_B, legacy_errorB, 
                       pValAB, pValA, pValB, FValAB, FValA, FValB)
  write.csv(report, paste0("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N",
                           N, "_eqVar", eqVar, "_eqCov", eqCov, "_null", simNull, ".csv"), row.names=F)
}


###################
## VISUALIZATION ##
###################

library(dplyr) # v 1.1.4
library(ggplot2) # v 3.5.2

## Equal Variances, Unequal Covariances, Simulated under the Alternative
ggplot(filter(report_vTcF_H1, !if_any(contains("error"), ~ . > .05)),
       aes(-log(legacy01_B), -log(BF01_B),
           color=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
           shape=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.7) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("Random Intercepts Only  ln"~italic("BF")[10]),
       y=expression("Maximal Set of Random Effects  ln"~italic("BF")[10]),
       color=expression(paste(italic("p-"), "value")),
       shape=expression(paste(italic("p-"), "value")),
       title=expression("Log Bayes Factors for Testing the Main Effect"~italic("B")),
       subtitle=paste0(N, " Participants and 1 Trial"),
       caption=expression(paste("Data-generating: equal variances, unequal covariances. Not significant ",italic("AB"),", significant ",italic(A)," and ",italic(B),"."))) +
  scale_x_continuous(limits=c(-2, 10), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3), 5, 10),
                     labels=c("-ln3", 0, "ln3", 5, 10)) +
  scale_y_continuous(limits=c(-2, 10), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3), 5, 10),
                     labels=c("-ln3", 0, "ln3", 5, 10)) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  coord_fixed(ratio=1) +
  theme_classic() +
  theme(text=element_text(size=14),
        legend.title=element_text(margin=margin(0,0.3,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.3,0,0.2,"cm")),
        legend.position="bottom", legend.direction="horizontal") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))

ggsave("~/Downloads/osfstorage-archive/Simulation 1/supp/fig_vTcF_H1.png", width=2250, height=2500, units="px", device="png", dpi=300)


## Equal Variances, Unequal Covariances, Simulated under the Null
ggplot(filter(report_vTcF_H0, !if_any(contains("error"), ~ . > .05)),
       aes(log(legacy01_B), log(BF01_B),
           color=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
           shape=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.7) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("Random Intercepts Only  ln"~italic("BF")["01"]),
       y=expression("Maximal Set of Random Effects  ln"~italic("BF")["01"]),
       color=expression(paste(italic("p-"), "value")),
       shape=expression(paste(italic("p-"), "value")),
       title=expression("Log Bayes Factors for Testing the Main Effect"~italic("B")),
       subtitle=paste0(N, " Participants and 1 Trial"),
       caption=expression(paste("Data-generating: equal variances, unequal covariances. Not significant ",italic("AB")," and ",italic(B),", significant ",italic(A),"."))) +
  scale_x_continuous(limits=c(-2, 2), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  scale_y_continuous(limits=c(-2, 2), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  coord_fixed(ratio=1) +
  theme_classic() +
  theme(text=element_text(size=14),
        legend.title=element_text(margin=margin(0,0.3,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.3,0,0.2,"cm")),
        legend.position="bottom", legend.direction="horizontal") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))

ggsave("~/Downloads/osfstorage-archive/Simulation 1/supp/fig_vTcF_H0.png", width=2250, height=2500, units="px", device="png", dpi=300)


## Equal Variances and Covariances, Simulated under the Alternative
ggplot(filter(report_vTcT_H1, !if_any(contains("error"), ~ . > .05)),
       aes(-log(legacy01_B), -log(BF01_B),
           color=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
           shape=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.7) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("Random Intercepts Only  ln"~italic("BF")[10]),
       y=expression("Maximal Set of Random Effects  ln"~italic("BF")[10]),
       color=expression(paste(italic("p-"), "value")),
       shape=expression(paste(italic("p-"), "value")),
       title=expression("Log Bayes Factors for Testing the Main Effect"~italic("B")),
       subtitle=paste0(N, " Participants and 1 Trial"),
       caption=expression(paste("Data-generating: equal variances and covariances. Not significant ",italic("AB"),", significant ",italic(A)," and ",italic(B),"."))) +
  scale_x_continuous(limits=c(-2, 12), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  scale_y_continuous(limits=c(-2, 12), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  coord_fixed(ratio=1) +
  theme_classic() +
  theme(text=element_text(size=14),
        legend.title=element_text(margin=margin(0,0.3,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.3,0,0.2,"cm")),
        legend.position="bottom", legend.direction="horizontal") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))

ggsave("~/Downloads/osfstorage-archive/Simulation 1/supp/fig_vTcT_H1.png", width=2250, height=2500, units="px", device="png", dpi=300)


## Equal Variances and Covariances, Simulated under the Null
ggplot(filter(report_vTcT_H0, !if_any(contains("error"), ~ . > .05)),
       aes(log(legacy01_B), log(BF01_B),
           color=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
           shape=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.7) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("Random Intercepts Only  ln"~italic("BF")["01"]),
       y=expression("Maximal Set of Random Effects  ln"~italic("BF")["01"]),
       color=expression(paste(italic("p-"), "value")),
       shape=expression(paste(italic("p-"), "value")),
       title=expression("Log Bayes Factors for Testing the Main Effect"~italic("B")),
       subtitle=paste0(N, " Participants and 1 Trial"),
       caption=expression(paste("Data-generating: equal variances and covariances. Not significant ",italic("AB")," and ",italic(B),", significant ",italic(A),"."))) +
  scale_x_continuous(limits=c(-2.2, 2), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  scale_y_continuous(limits=c(-2.2, 2), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  coord_fixed(ratio=1) +
  theme_classic() +
  theme(text=element_text(size=14),
        legend.title=element_text(margin=margin(0,0.3,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.3,0,0.2,"cm")),
        legend.position="bottom", legend.direction="horizontal") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))

ggsave("~/Downloads/osfstorage-archive/Simulation 1/supp/fig_vTcT_H0.png", width=2250, height=2500, units="px", device="png", dpi=300)


## Unequal Variances and Covariances, Simulated under the Alternative
ggplot(filter(report_vFcF_H1, !if_any(contains("error"), ~ . > .05)),
       aes(-log(legacy01_B), -log(BF01_B),
           color=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
           shape=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.7) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("Random Intercepts Only  ln"~italic("BF")[10]),
       y=expression("Maximal Set of Random Effects  ln"~italic("BF")[10]),
       color=expression(paste(italic("p-"), "value")),
       shape=expression(paste(italic("p-"), "value")),
       title=expression("Log Bayes Factors for Testing the Main Effect"~italic("B")),
       subtitle=paste0(N, " Participants and 1 Trial"),
       caption=expression(paste("Data-generating: unequal variances and covariances. Not significant ",italic("AB"),", significant ",italic(A)," and ",italic(B),"."))) +
  scale_x_continuous(limits=c(-2, 20), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3), 10, 20),
                     labels=c("-ln3", 0, "ln3", 10, 20)) +
  scale_y_continuous(limits=c(-2, 20), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3), 10, 20),
                     labels=c("-ln3", 0, "ln3", 10, 20)) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  coord_fixed(ratio=1) +
  theme_classic() +
  theme(text=element_text(size=14),
        legend.title=element_text(margin=margin(0,0.3,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.3,0,0.2,"cm")),
        legend.position="bottom", legend.direction="horizontal") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))

ggsave("~/Downloads/osfstorage-archive/Simulation 1/supp/fig_vFcF_H1.png", width=2250, height=2500, units="px", device="png", dpi=300)


## Unequal Variances and Covariances, Simulated under the Null
ggplot(filter(report_vFcF_H0, !if_any(contains("error"), ~ . > .05)),
       aes(log(legacy01_B), log(BF01_B),
           color=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
           shape=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.7) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("Random Intercepts Only  ln"~italic("BF")["01"]),
       y=expression("Maximal Set of Random Effects  ln"~italic("BF")["01"]),
       color=expression(paste(italic("p-"), "value")),
       shape=expression(paste(italic("p-"), "value")),
       title=expression("Log Bayes Factors for Testing the Main Effect"~italic("B")),
       subtitle=paste0(N, " Participants and 1 Trial"),
       caption=expression(paste("Data-generating: unequal variances and covariances. Not significant ",italic("AB")," and ",italic(B),", significant ",italic(A),"."))) +
  scale_x_continuous(limits=c(-3, 2), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  scale_y_continuous(limits=c(-3, 2), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  coord_fixed(ratio=1) +
  theme_classic() +
  theme(text=element_text(size=14),
        legend.title=element_text(margin=margin(0,0.3,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.3,0,0.2,"cm")),
        legend.position="bottom", legend.direction="horizontal") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))

ggsave("~/Downloads/osfstorage-archive/Simulation 1/supp/fig_vFcF_H0.png", width=2250, height=2500, units="px", device="png", dpi=300)

