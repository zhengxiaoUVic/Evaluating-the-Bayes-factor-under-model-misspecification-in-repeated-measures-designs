##################################################################### SIMULATION
library(fGarch)      # v 4052.93
library(BayesFactor) # v 0.9.12-4.7

## constant simulation settings
C1 <- C2 <- c(-.5, .5); B2 <- .5; NSim <- 200
sd <- 0.5 # between-subjects standard deviation
h <- 0.5 # the prior scale on the variability of fixed effects for computing Bayes factors
params <- data.frame("B1"=rep(c(0, .25, .5, .75), each=4), # effect size
                     "Case"=rep(1:4, 4)) # Case 1: normality holds
                                         #      2: random effects are skewed
                                         #      3: residuals are skewed
                                         #      4: both are skewed

## varying simulation settings
Sk <- 1/20 # skewness parameter, {20, 1/20}, a measure of the asymmetry; Sk > 1 right-skewed; 0 < Sk < 1 left-skewed
N <- 48 # number of participants, {12, 24, 48}
n <- 10 # number of trials, {1, 4, 10}
re <- ifelse(F, -1, 1) # change random seeds for the replicate cases

for (rindex in 1:16) { #1:4  5:8  9:12  13:16
  list2env(params[rindex,], envir=.GlobalEnv) # assign parameters
  BF01_AB <- BF_errorAB <- BF01_A <- BF_errorA <- BF01_B <- BF_errorB <-
    BF01_0 <- BF_error0 <-
    pValAB <- pValA <- pValB <- numeric(NSim)
  
  for (i in 1:NSim) {
    if (i %% 10 == 0) print (paste0("Run ", i, ": ", Sys.time()))
    
    set.seed(re*(i+NSim*rindex))
    if (Case %in% c(1,3)) {
      s <- rnorm(N, 0, sd)   # random intercepts
      b1 <- rnorm(N, B1, sd) # random slopes of factor A (this is the SECOND factor in the paper)
      b2 <- rnorm(N, B2, sd) # random slopes of factor B (this is the FIRST factor in the paper)
    } else {
      s <- rsstd(N, mean=0, sd=sd, nu=4, xi=Sk)
      b1 <- rsstd(N, mean=B1, sd=sd, nu=4, xi=Sk)
      b2 <- rsstd(N, mean=B2, sd=sd, nu=4, xi=Sk)
    }
    
    M <- s + cbind(b1, b2) %*% t(expand.grid(C1, C2)) # cell means
    if (Case %in% c(1,2)) {
      DV <- rnorm(N*4*n, c(M), 1) # a1b1, a2b1, a1b2, a2b2
    } else {
      DV <- rsstd(N*4*n, mean=c(M), sd=1, nu=4, xi=Sk) # residual standard deviation = 1
    }
    dat <- data.frame("DV" = DV,
                      "ID" = rep(paste0("s",1:N), 4*n),
                      "A" = rep(c("a1","a2"), each=N, times=2*n),
                      "B" = rep(c("b1","b2"), each=N*2, times=n),
                      "Trial" = rep(paste0("r",1:n), each=N*4),
                      stringsAsFactors = T) # long format data
    
    RMANOVA <- aov(DV ~ (A * B) + Error(ID / (A * B)), dat) # repeated-measures analysis of variance
    sum_stat <- unlist(summary(RMANOVA))
    pValAB[i] <- unname(sum_stat["Error: ID:A:B.Pr(>F)1"])
    pValA[i] <- unname(sum_stat["Error: ID:A.Pr(>F)1"])
    pValB[i] <- unname(sum_stat["Error: ID:B.Pr(>F)1"])
    
    if (n > 1) {
      set.seed(277)
      BF_full <- lmBF(DV ~ A * B + ID + A:ID + B:ID + A:B:ID, whichRandom="ID", rscaleFixed=h,
                      data=dat, iterations=1e5, progress=F) # maximal set of random effects (MRE) model
      # the highest-order interaction A:B:ID is only estimable when n > 1
      
      set.seed(277)
      BF_add <- lmBF(DV ~ A + B + ID + A:ID + B:ID + A:B:ID, whichRandom="ID", rscaleFixed=h,
                     data=dat, iterations=1e5, progress=F)
      
      set.seed(277)
      BF_omitA <- lmBF(DV ~ B + ID + A:ID + B:ID + A:B:ID, whichRandom="ID", rscaleFixed=h,
                       data=dat, iterations=1e5, progress=F)
      
      set.seed(277)
      BF_omitB <- lmBF(DV ~ A + ID + A:ID + B:ID + A:B:ID, whichRandom="ID", rscaleFixed=h,
                       data=dat, iterations=1e5, progress=F)
      
      set.seed(277)
      BF_0 <- lmBF(DV ~ ID + A:ID + B:ID + A:B:ID, whichRandom="ID", rscaleFixed=h,
                   data=dat, iterations=1e5, progress=F)
      
    } else { # aggregated data or n==1
      set.seed(277)
      BF_full <- lmBF(DV ~ A * B + ID + A:ID + B:ID, whichRandom="ID", rscaleFixed=h,
                      data=dat, iterations=1e5, progress=F)
      
      set.seed(277)
      BF_add <- lmBF(DV ~ A + B + ID + A:ID + B:ID, whichRandom="ID", rscaleFixed=h,
                     data=dat, iterations=1e5, progress=F)
      
      set.seed(277)
      BF_omitA <- lmBF(DV ~ B + ID + A:ID + B:ID, whichRandom="ID", rscaleFixed=h,
                       data=dat, iterations=1e5, progress=F)
      
      set.seed(277)
      BF_omitB <- lmBF(DV ~ A + ID + A:ID + B:ID, whichRandom="ID", rscaleFixed=h,
                       data=dat, iterations=1e5, progress=F)
      
      set.seed(277)
      BF_0 <- lmBF(DV ~ ID + A:ID + B:ID, whichRandom="ID", rscaleFixed=h,
                   data=dat, iterations=1e5, progress=F)
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
    
    bf <- extractBF(BF_0 / BF_add)
    BF01_0[i] <- bf$bf
    BF_error0[i] <- bf$error
  }
  
  if (re == -1 & Case == 1) Case <- 0 # replicates of Case 1
  if (Sk < 1 & Case > 1) Case <- Case + 3 # left-skewed case; otherwise, right-skewed
  report <- data.frame(pValAB, pValA, pValB,
                       BF01_AB, BF_errorAB,
                       BF01_A, BF_errorA,
                       BF01_B, BF_errorB,
                       rindex, Case, B1, h, Sk,
                       BF01_0, BF_error0)
  
  write.csv(report,
            paste0("~/Downloads/osfstorage-archive/Simulation 2/data/Sim_N",N,"_n",n,"_ES",B1,"_Case",Case,".csv"),
            row.names=F)
} # roughly 12 hours of run time for each `rindex`


############################################################### DATA PREPARATION
dat <- c()
for (nSubj in c(48, 24, 12)) { # number of participants
  for (nTrial in c(1, 4, 10)) { # number of trials
    for (es in c(0, 0.25, 0.5, 0.75)) { # effect sizes
      for (case in 0:7) {
        temp <- read.csv(paste0("~/Downloads/osfstorage-archive/Simulation 2/data/Sim_N",nSubj,"_n",nTrial,"_ES",es,"_Case",case,".csv"))
        temp$N <- nSubj
        temp$n <- nTrial
        dat <- rbind(dat, temp)
      }
    }
  }
}
write.csv(dat, paste0("~/Downloads/osfstorage-archive/Simulation 2/Figure 3/Sim2x2_robust.csv"), row.names=F) # 13.1 MB


################################################################## VISUALIZATION
library(dplyr) # v 1.1.4
library(ggplot2) # v 4.0.1
library(transport) # v 0.15-4

dat <- read.csv("~/Downloads/osfstorage-archive/Simulation 2/Figure 3/Sim2x2_robust.csv")
dat <- filter(dat, !if_any(contains("error"), ~ . > .05)) # 53 (out of 57,600) have been filtered out
dat$Label <- factor(dat$Case, labels=c("0: Normality holds",
                                       "1: Normality holds",
                                       "2: Random effects are right-skewed",
                                       "3: Residuals are right-skewed",
                                       "4: Both are right-skewed",
                                       "5: Random effects are left-skewed",
                                       "6: Residuals are left-skewed",
                                       "7: Both are left-skewed"))
dat$ES <- factor(dat$B1, labels=paste0(c("Null", paste0("Positive ", c("Small", "Medium", "Large"))), " ES"))

for (nSubj in c(48, 24, 12)) { # number of participants
  for (nTrial in c(1, 4, 10)) { # number of trials
    
    df <- subset(dat, N==nSubj & n==nTrial)
    
    for (BMA in c(F, T)) { # whether to use the inclusion Bayes factor
      if (BMA) {
        df$BF01_A <- (df$BF01_A + df$BF01_0) / (1 + df$BF01_B)
      }
      
      takeLog <- T
      if (takeLog) {
        df$BF <- -log(df$BF01_A) # ln(BF₁₀)
        refLines <- c(-log(3),log(3))
        yLab <- bquote(.(ifelse(BMA, "Inclusion  ", "")) * "MRE  ln " * italic("BF")[10])
      } else {
        df$BF <- 1/df$BF01_A # BF₁₀
        refLines <- c(1/3,3)
        yLab <- bquote(.(ifelse(BMA, "Inclusion  ", "")) * "MRE " * italic("BF")[10])
      }
      
      distBF <- distP <- ksBF <- ksP <- prop <- numeric(32); i <- 1
      for (es in c(0, 0.25, 0.5, 0.75)) {
        temp <- subset(df, B1==es & Case==0)
        for (case in 0:7) {
          sub <- subset(df, B1==es & Case==case)
          prop[i] <- sum(sub$BF01_A < 1/3) / nrow(sub) # proportion of the Bayes factors supporting the alternative
          distBF[i] <- ifelse(case==0, "EMD", sprintf("%.2f", wasserstein1d(temp$BF, sub$BF)))
          distP[i] <- ifelse(case==0, "EMD", round(wasserstein1d(temp$pValA, sub$pValA),10)) # earth mover's distance (EMD)
          ksBF[i] <- ifelse(case==0, "KS", round(ks.test(temp$BF, sub$BF)$statistic,2))
          ksP[i] <- ifelse(case==0, "KS", round(ks.test(temp$pValA, sub$pValA)$statistic,3)) # Kolmogorov–Smirnov (KS) statistic
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
      
      fig_BF <- ggplot(df, aes(x=Case, y=-log(BF01_A), fill=Label)) +
        geom_boxplot(alpha=.4) +
        facet_wrap(.~ES, scales="free_y", nrow=2, ncol=2) +
        labs(x="", y=yLab, fill="Case",
             title=paste0(ifelse(takeLog,"Log ",""),"Bayes Factors for Testing the Main Effect in 2×2 Repeated-Measures Designs"),
             subtitle=paste0(nSubj, " Participants and ", nTrial, ifelse(nTrial > 1, " Trials", " Trial"))) +
        scale_fill_manual(values=c("#d83034", "#d83034", "#f9e858", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
        scale_y_continuous(n.breaks=3) +
        scale_x_continuous(breaks=seq(0,7,by=1)) +
        geom_hline(yintercept=refLines, linetype="dashed", color="gray") +
        geom_text(data=anno, aes(x=Case, y=Inf, label=distBF), vjust=36, col="red", fontface="bold") +
        theme_classic() +
        theme(legend.position="bottom", text=element_text(size=17), legend.margin=margin(t=-15),
              legend.title=element_text(margin=margin(0,1,0,0,"cm")),
              legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm"))) +
        guides(fill=guide_legend(nrow=3, byrow=T))
      ggsave(paste0("~/Downloads/osfstorage-archive/Simulation 2/supp/fig_",ifelse(BMA,"i",""),"BF_N",nSubj,"_n",nTrial,".png"), fig_BF, width=3800, height=3200, units="px", device="png", dpi=300)
    }
    
    
    fig_pval <- ggplot(df, aes(x=Case, y=pValA, fill=Label)) +
      geom_boxplot(alpha=.4, outliers=F) +
      facet_wrap(.~ES, scales="free_y", nrow=2, ncol=2) +
      labs(x="", y=expression(~italic("p")~"-Value"), fill="Case",
           title="Null Hypothesis Significance Testing in 2×2 rANOVA",
           subtitle=paste0(nSubj, " Participants and ", nTrial, " Trials")) +
      scale_fill_manual(values=c("#d83034", "#d83034", "#f9e858", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
      scale_x_continuous(breaks=seq(0,7,by=1)) +
      theme_classic() +
      theme(legend.position="bottom", text=element_text(size=17), legend.margin=margin(t=-15),
            legend.title=element_text(margin=margin(0,1,0,0,"cm")),
            legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm"))) +
      guides(fill=guide_legend(nrow=3, byrow=T))
    ggsave(paste0("~/Downloads/osfstorage-archive/Simulation 2/supp/fig_pVal_N",nSubj,"_n",nTrial,".png"), fig_pval, width=3800, height=3200, units="px", device="png", dpi=300)
    
  }
}

