library(dplyr) # v 1.1.4
library(ggh4x) # v 0.3.1; library(ggplot2) # v 4.0.1
library(transport) # v 0.15-4

dat <- read.csv("~/Downloads/osfstorage-archive/Simulation 2/Figure 3/Sim2x2_robust.csv")
dat <- filter(dat, !if_any(contains("error"), ~ . > .05)) # 53 (out of 36,000) have been filtered out
dat$Label <- factor(dat$Case, labels=c("0: Normality holds",
                                       "1: Normality holds",
                                       "2: Random effects are right-skewed",
                                       "3: Residuals are right-skewed",
                                       "4: Both are right-skewed",
                                       "5: Random effects are left-skewed",
                                       "6: Residuals are left-skewed",
                                       "7: Both are left-skewed"))
dat$ES <- factor(dat$B1, labels=paste0(c("Null", paste0("Positive ", c("Small", "Medium", "Large"))), " ES"))


nSubj <- 48 # number of participants
nTrial <- 10 # number of trials
df <- subset(dat, N==nSubj & n==nTrial)

distBF <- distP <- ksBF <- ksP <- prop <- numeric(32); i <- 1
for (es in c(0, .25, .5, .75)) {
  temp <- subset(df, B1==es & Case==0)
  for (case in 0:7) {
    sub <- subset(df, B1==es & Case==case)
    prop[i] <- sum(sub$BF01_A < 1/3) / nrow(sub) # proportion of the Bayes factors supporting the alternative
    distBF[i] <- ifelse(case==0, "EMD", sprintf("%.2f", wasserstein1d(-log(temp$BF01_A), -log(sub$BF01_A)))) # ln(BF₁₀)
    distP[i] <- ifelse(case==0, "EMD", round(wasserstein1d(temp$pValA, sub$pValA),10)) # earth mover's distance (EMD)
    ksBF[i] <- ifelse(case==0, "KS", round(ks.test(-log(temp$BF01_A), -log(sub$BF01_A))$statistic,2)) # ln(BF₁₀)
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

fig2 <- ggplot(subset(df, ES %in% c("Null ES", "Positive Medium ES")), aes(x=Case, y=-log(BF01_A), fill=Label)) + # ln(BF₁₀)
  geom_boxplot(alpha=.4) +
  facet_wrap(.~ES, scales="free_y", nrow=1) +
  labs(x="", y=expression("MRE "~italic("BF")[10]), fill="Case",
       title=bquote("2×2 Repeated-Measures Designs: " * italic("n") * "=" * .(nSubj) * ", " * italic("R") * "=" * .(nTrial))) +
  scale_fill_manual(values=c("#d83034", "#d83034", "#f9e858", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  facetted_pos_scales(y=list(scale_y_continuous(breaks=c(-log(10), -log(3), 0, log(3), log(10)),
                                                labels=c("1/10", "1/3", 1, 3, 10)),
                             scale_y_continuous(breaks=c(-log(3), 0, log(3), log(1E4), log(1E8)),
                                                labels=c("1/3", 1, 3, expression(10^4), expression(10^8))) )) +
  scale_x_continuous(breaks=seq(0,7,by=1)) +
  geom_hline(yintercept=c(-log(3), log(3)), linetype="dashed", color="gray") +
  geom_text(data=subset(anno, ES %in% c("Null ES", "Positive Medium ES")), 
            aes(x=Case, y=Inf, label=distBF), vjust=49.5, col="red", fontface="bold") +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=20), legend.margin=margin(t=-15),
        legend.title=element_text(margin=margin(0,1,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm"))) +
  guides(fill=guide_legend(nrow=3, byrow=T))
ggsave("~/Downloads/osfstorage-archive/Simulation 2/Figure 3/Figure 3.png", fig2, width=3800, height=2300, units="px", device="png", dpi=300)

