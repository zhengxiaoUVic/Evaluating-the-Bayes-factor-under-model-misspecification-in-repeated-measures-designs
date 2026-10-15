##################################################################### SIMULATION
library(fGarch)      # v 4033.92
library(BayesFactor) # v 0.9.12-4.7

## constant simulation settings
NSim <- 200 # number of simulation runs for each setting
sd <- 0.5 # between-subjects standard deviation
h <- 0.5 # the prior scale on the variability of fixed effects for computing Bayes factors
params <- data.frame("B1"=rep(c(0, .25, .5, .75), each=4), # effect size
                     "Case"=rep(1:4, 4)) # Case 1: normality holds
                                         #      2: random intercepts are skewed
                                         #      3: residuals are skewed
                                         #      4: both are skewed

## varying simulation settings
Sk <- 1/20 # skewness parameter, {20, 1/20}, a measure of the asymmetry; Sk > 1 right-skewed; 0 < Sk < 1 left-skewed
N <- 48 # number of participants, {12, 24, 48, 100}
re <- ifelse(F, -1, 1) # run for the replicate cases

for (rindex in 1:16) { #1:4  5:8  9:12  13:16
  list2env(params[rindex,], envir=.GlobalEnv) # assign parameters
  BF10 <- BF_error <- BF10t <- BFt_error <- pVal <- numeric(NSim)
  seed <- re * sample(1:1E5, NSim) # random seeds
  
  for (i in 1:NSim) {
    if (i %% 10 == 0) print (paste0("Run ", i, ": ", Sys.time()))
    
    set.seed(seed[i])
    if (Case %in% c(1,3)) {
      s <- rnorm(N, 0, sd)   # subject-specific random effects
    } else {
      s <- rsstd(N, mean=0, sd=sd, nu=4, xi=Sk)
    }
    
    if (Case %in% c(1,2)) {
      DV <- rnorm(N*2, B1 * c(-.5, .5), 1) # a1, a2, between-subjects
    } else {
      DV <- rsstd(N*2, mean=B1 * c(-.5, .5), sd=1, nu=4, xi=Sk) # residual standard deviation = 1
    }
    dat <- data.frame("DV" = c(sweep(matrix(DV, ncol=2, byrow=T), 1, s, "+")), # within-subject
                      "ID" = rep(paste0("s",1:N), 2),
                      "A" = rep(c("a1","a2"), each=N),
                      stringsAsFactors = T) # long format data
    
    fit <- aov(DV ~ A + Error(ID/A), dat) # one-way repeated-measures analysis of variance
    pVal[i] <- summary(fit)[[2]][[1]]["A", "Pr(>F)"] # p-value
    # t.test(subset(dat, A=="a1")$DV, subset(dat, A=="a2")$DV, paired=T)$p.value # SAME
    
    set.seed(277)
    bf <- anovaBF(DV ~ A + ID, dat, whichRandom="ID", rscaleFixed=h, progress=F)
    BF10[i] <- extractBF(bf)$bf
    BF_error[i] <- extractBF(bf)$error
    
    bf_t <- ttestBF(dat$DV[1:N], dat$DV[(N+1):(2*N)], paired=T)
    BF10t[i] <- extractBF(bf_t)$bf
    BFt_error[i] <- extractBF(bf_t)$error
  }
  
  if (re == -1 & Case == 1) Case <- 0 # replicates of Case 1
  if (Sk < 1 & Case > 1) Case <- Case + 3 # left-skewed case; otherwise, right-skewed
  report <- data.frame(pVal, BF10, BF_error, BF10t, BFt_error, seed,
                       rindex, Case, B1, h, Sk)
  
  write.csv(report, paste0("~/Downloads/osfstorage-archive/Simulation 3 continued/data/Sim_1way_N",N,"_ES",B1,"_Case",Case,".csv"),
            row.names=F)
} # roughly 40 seconds of run time for each `rindex`


############################################################### DATA PREPARATION
dat <- c()
for (nSubj in c(100, 48, 24, 12)) { # number of participants
  for (es in c(0, 0.25, 0.5, 0.75)) { # effect sizes
    for (case in 0:7) {
      temp <- read.csv(paste0("~/Downloads/osfstorage-archive/Simulation 3 continued/data/Sim_1way_N",nSubj,"_ES",es,"_Case",case,".csv"))
      temp$N <- nSubj
      dat <- rbind(dat, temp)
    }
  }
}
write.csv(dat, paste0("~/Downloads/osfstorage-archive/Simulation 3 continued/Figure 5/Sim_1way_robust.csv"), row.names=F) # 884 KB


################################################################## VISUALIZATION
library(dplyr) # v 1.1.4
library(ggplot2) # v 3.5.1
library(transport) # v 0.15-4

dat <- read.csv("~/Downloads/osfstorage-archive/Simulation 3 continued/Figure 5/Sim_1way_robust.csv")
dat <- filter(dat, !if_any(contains("error"), ~ . > .1)) # 0 (out of 12,000) have been filtered out
dat$Label <- factor(dat$Case, labels=c("0: Normality holds",
                                       "1: Normality holds",
                                       "2: Random intercepts are right-skewed",
                                       "3: Residuals are right-skewed",
                                       "4: Both are right-skewed",
                                       "5: Random intercepts are left-skewed",
                                       "6: Residuals are left-skewed",
                                       "7: Both are left-skewed"))
dat$ES <- factor(dat$B1, labels=paste0(c("Null", paste0("Positive ", c("Small", "Medium", "Large"))), " ES"))


nSubj <- 100 # number of participants
df <- subset(dat, N==nSubj)
takeLog <- T
if (takeLog) {
  df$BF <- log(df$BF10) # ln(BF₁₀)
  refLines <- c(-log(3),log(3))
  yLab <- expression("ln"~italic("BF")[10])
} else {
  df$BF <- df$BF10 # BF₁₀
  refLines <- c(1/3,3)
  yLab <- expression(italic("BF")[10])
}

distBF <- distP <- ksBF <- ksP <- prop <- numeric(32); i <- 1
for (es in c(0, 0.25, 0.5, 0.75)) {
  temp <- subset(df, B1==es & Case==0)
  for (case in 0:7) {
    sub <- subset(df, B1==es & Case==case)
    prop[i] <- sum(sub$BF10 > 3) / nrow(sub) # proportion of the Bayes factors supporting the alternative
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

fig_BF <- ggplot(df, aes(x=Case, y=log(BF10), fill=Label)) +
  geom_boxplot(alpha=.4) +
  facet_wrap(.~ES, scales="free_y", nrow=2, ncol=2) +
  labs(x="", y=yLab, fill="Case",
       title="Bayes Factors for Testing Effects in One-Way Repeated-Measures Designs",
       subtitle=paste0("Number of Participants: ", nSubj)) +
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
ggsave(paste0("~/Downloads/osfstorage-archive/Simulation 3 continued/supp/fig_1way_BF_N",nSubj,".png"), fig_BF, width=3800, height=3200, units="px", device="png", dpi=300)


fig_pval <- ggplot(df, aes(x=Case, y=pVal, fill=Label)) +
  geom_boxplot(alpha=.4, outliers=F) +
  facet_wrap(.~ES, scales="free_y", nrow=2, ncol=2) +
  labs(x="", y=expression(~italic("p")~"-Value"), fill="Case",
       title="Null Hypothesis Significance Testing",
       subtitle=paste0("Number of Participants: ", nSubj)) +
  scale_fill_manual(values=c("#d83034", "#d83034", "#f9e858", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  scale_x_continuous(breaks=seq(0,7,by=1)) +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=17), legend.margin=margin(t=-15),
        legend.title=element_text(margin=margin(0,1,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm"))) +
  guides(fill=guide_legend(nrow=3, byrow=T))
ggsave(paste0("~/Downloads/osfstorage-archive/Simulation 3 continued/supp/fig_1way_pVal_N",nSubj,".png"), fig_pval, width=3800, height=3200, units="px", device="png", dpi=300)


set.seed(277)
ggplot(df[sample(1:nrow(df), nrow(df)),], # shuffle to prevent overlapping
       aes(log(BF10t), log(BF10),
           color=factor(cut(pVal, breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
           shape=factor(cut(pVal, breaks=c(-1E-5, .01, .05, .1, 1),
                            labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.3) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("Paired"~italic(t) * "-Test  ln"~italic("BF")[10]),
       y=expression("One-Way rANOVA  ln"~italic("BF")[10]),
       color=expression(paste(italic("p-"), "value")),
       shape=expression(paste(italic("p-"), "value")),
       title="Bayes Factors for Testing Effects",
       subtitle=paste0("Number of Participants: ", nSubj)) +
  scale_x_continuous(limits=c(-2, 20), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3), 10, 15),
                     labels=c("-ln3", 0, "ln3", 10, 15)) +
  scale_y_continuous(limits=c(-2, 20), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3), 10, 15),
                     labels=c("-ln3", 0, "ln3", 10, 15)) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  coord_fixed(ratio=1) +
  theme_classic() +
  theme(text=element_text(size=14),
        legend.position="bottom",legend.spacing.x=unit(0.5,"cm"),
        legend.direction="horizontal") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))

ggsave(paste0("~/Downloads/osfstorage-archive/Simulation 3 continued/supp/fig_cf_BF_N",nSubj,".png"), width=3800, height=2300, units="px", device="png", dpi=300)

