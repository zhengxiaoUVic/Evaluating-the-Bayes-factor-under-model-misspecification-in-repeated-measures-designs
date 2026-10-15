##################################################################### SIMULATION
library(EnvStats)    # v 3.0.0
library(fGarch)      # v 4033.92
library(BayesFactor) # v 0.9.12-4.7

## constant simulation settings
NSim <- 200 # number of simulation runs for each setting
nu <- 4 # shape parameter (degrees of freedom) of the skewed-t
mu.diff <- 0.45 # unweighted mean difference of the two modals
p.mix <- 1/3 # weight of the second modal
SD <- 1 # population standard deviation
var_eq <- SD^2-mu.diff^2/(p.mix*(1-p.mix))
if(var_eq < 0) {
  stop("Such a bimodal cannot achieve the same population variance.")
}
params <- data.frame("B1"=rep(c(0, .5), each=6), # effect size
                     "Case"=rep(1:6, 2)) # Case 1: normality holds
                                         #      2: uniform
                                         #      3: symmetric triangular
                                         #      4: symmetric bimodal
                                         #      5: asymmetric bimodal
                                         #      6: skewed-t

## varying simulation settings
Sk <- 20 # skewness parameter, {20, 1/20}, a measure of the asymmetry; Sk > 1 right-skewed; 0 < Sk < 1 left-skewed
N <- 100 # number of participants, {12, 24, 48, 100}
re <- ifelse(F, -1, 1) # run for the replicate cases

for (rindex in 1:12) {
  list2env(params[rindex,], envir=.GlobalEnv) # assign parameters
  BF10t <- BFt_error <- pVal <- numeric(NSim)
  seed <- re * sample(1:1E5, NSim) # random seeds
  
  for (i in 1:NSim) {
    if (i %% 100 == 0) print (paste0("Run ", i, ": ", Sys.time()))
    
    set.seed(seed[i])
    if (Case == 1) {
      y <- rnorm(N, B1, SD) # normal
    } else if (Case == 2) {
      # mean = (a + b) / 2
      # variance = (b - a)^2 / 12
      y <- runif(N, min=B1-sqrt(3)*SD, max=B1+sqrt(3)*SD) # uniform
    } else if (Case == 3) {
      a <- B1 - sqrt(6) * SD
      b <- B1 + sqrt(6) * SD
      # mean = (a + b) / 2
      # variance = (b - a)^2 / 24
      u <- runif(N) # inverse transformation method
      y <- ifelse(u < 0.5,
                  a+sqrt((b-a)^2*u/2),
                  b-sqrt((b-a)^2*(1-u)/2)) # symmetric triangular
    } else if (Case == 4) {
      y <- rnormMix(N, mean1=B1-mu.diff*2, mean2=B1+mu.diff*2, 
                    sd1=sqrt(SD^2-(mu.diff*2)^2), sd2=sqrt(SD^2-(mu.diff*2)^2),
                    p.mix=0.5) # symmetric bimodal
    } else if (Case == 5) {
      y <- rnormMix(N, mean1=B1-mu.diff/(1-p.mix), mean2=B1+mu.diff/p.mix, 
                    sd1=sqrt(var_eq), sd2=sqrt(var_eq), p.mix=p.mix) # asymmetric bimodal
    } else {
      y <- rsstd(N, mean=B1, sd=SD, nu=nu, xi=Sk) # skewed-t
      if (Sk < 1) Case <- 7 # left-skewed case; otherwise, right-skewed
    }
    
    pVal[i] <- t.test(y)$p.value # one-sample t-test p-value
    set.seed(277)
    bf <- ttestBF(y)
    BF10t[i] <- extractBF(bf)$bf
    BFt_error[i] <- extractBF(bf)$error
  }
  
  if (re == -1 & Case == 1) Case <- 0 # replicates of Case 1
  report <- data.frame(pVal, BF10t, BFt_error, seed,
                       rindex, Case, B1)
  
  write.csv(report, paste0("~/Downloads/osfstorage-archive/Simulation 3/data/Sim_1t_N",N,"_ES",B1,"_Case",Case,".csv"),
            row.names=F)
}


############################################################### DATA PREPARATION
dat <- c()
for (nSubj in c(100, 48, 24, 12)) { # number of participants
  for (es in c(0, 0.5)) { # effect sizes
    for (case in 0:7) {
      temp <- read.csv(paste0("~/Downloads/osfstorage-archive/Simulation 3/data/Sim_1t_N",nSubj,"_ES",es,"_Case",case,".csv"))
      temp$N <- nSubj
      dat <- rbind(dat, temp)
    }
  }
}
write.csv(dat, paste0("~/Downloads/osfstorage-archive/Simulation 3/Figure 4/Sim_1t_robust.csv"), row.names=F) # 884 KB


################################################################## VISUALIZATION
library(dplyr) # v 1.1.4
library(ggplot2) # v 3.5.1
library(transport) # v 0.15-4

dat <- read.csv("~/Downloads/osfstorage-archive/Simulation 3/Figure 4/Sim_1t_robust.csv")
dat <- filter(dat, !if_any(contains("error"), ~ . > .05)) # 0 (out of 12,800) have been filtered out
dat$Label <- factor(dat$Case, labels=c("0: Normality holds",
                                       "1: Normality holds",
                                       "2: Uniform",
                                       "3: Symmetric triangular",
                                       "4: Symmetric bimodal",
                                       "5: Asymmetric bimodal",
                                       "6: Right-skewed",
                                       "7: Left-skewed"))
dat$ES <- factor(dat$B1, labels=paste0(c("Null", "Positive Medium"), " ES"))


nSubj <- 100 # number of participants
df <- subset(dat, N==nSubj)
takeLog <- T
if (takeLog) {
  df$BF <- log(df$BF10t) # ln(BF₁₀)
  refLines <- c(-log(3),log(3))
  yLab <- expression("ln"~italic("BF")[10])
} else {
  df$BF <- df$BF10t # BF₁₀
  refLines <- c(1/3,3)
  yLab <- expression(italic("BF")[10])
}

distBF <- distP <- ksBF <- ksP <- prop <- numeric(16); i <- 1
for (es in c(0, 0.5)) {
  temp <- subset(df, B1==es & Case==0)
  for (case in 0:7) {
    sub <- subset(df, B1==es & Case==case)
    prop[i] <- sum(sub$BF10t > 3) / nrow(sub) # proportion of the Bayes factors supporting the alternative
    distBF[i] <- ifelse(case==0, "EMD", sprintf("%.2f", wasserstein1d(temp$BF, sub$BF)))
    distP[i] <- ifelse(case==0, "EMD", round(wasserstein1d(temp$pVal, sub$pVal),10)) # earth mover's distance (EMD)
    ksBF[i] <- ifelse(case==0, "KS", round(ks.test(temp$BF, sub$BF)$statistic,2))
    ksP[i] <- ifelse(case==0, "KS", round(ks.test(temp$pVal, sub$pVal)$statistic,3)) # Kolmogorov–Smirnov (KS) statistic
    i <- i + 1
  }
}
anno <- data.frame("Case"=0:7,
                   "Label"=factor(levels(df$Label)),
                   "ES"=factor(rep(levels(df$ES), each=8)),
                   "prop"=prop,
                   "distBF"=distBF,
                   "distP"=distP,
                   "ksBF"=ksBF,
                   "ksP"=ksP)

fig_BF <- ggplot(df, aes(x=Case, y=log(BF10t), fill=Label)) +
  geom_boxplot(alpha=.4) +
  facet_wrap(.~ES, scales="free_y", nrow=1) +
  labs(x="", y=yLab, fill="Case",
       title=expression("Bayes Factors for Testing Effects in Paired"~italic("t") * "-Tests"),
       subtitle=paste0("Number of Participants: ", nSubj)) +
  scale_fill_manual(values=c("#d83034", "#d83034", "#f9e858", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  scale_x_continuous(breaks=seq(0,7,by=1)) +
  geom_hline(yintercept=refLines, linetype="dashed", color="gray") +
  geom_text(data=anno, aes(x=Case, y=Inf, label=distBF), vjust=51, col="red", fontface="bold") +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=17), legend.margin=margin(t=-15),
        legend.title=element_text(margin=margin(0,1,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm"))) +
  guides(fill=guide_legend(nrow=2, byrow=T))
ggsave(paste0("~/Downloads/osfstorage-archive/Simulation 3/supp/fig_1t_BF_N",nSubj,".png"), fig_BF, width=3800, height=2300, units="px", device="png", dpi=300)


fig_pval <- ggplot(df, aes(x=Case, y=pVal, fill=Label)) +
  geom_boxplot(alpha=.4, outliers=F) +
  facet_wrap(.~ES, scales="free_y", nrow=1) +
  labs(x="", y=expression(~italic("p")~"-Value"), fill="Case",
       title="Null Hypothesis Significance Testing",
       subtitle=paste0("Number of Participants: ", nSubj)) +
  scale_fill_manual(values=c("#d83034", "#d83034", "#f9e858", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  scale_x_continuous(breaks=seq(0,7,by=1)) +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=17), legend.margin=margin(t=-15),
        legend.title=element_text(margin=margin(0,1,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm"))) +
  guides(fill=guide_legend(nrow=2, byrow=T))
ggsave(paste0("~/Downloads/osfstorage-archive/Simulation 3/supp/fig_1t_pVal_N",nSubj,".png"), fig_pval, width=3800, height=2300, units="px", device="png", dpi=300)

