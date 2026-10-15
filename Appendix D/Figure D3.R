library(dplyr) # v 1.1.4
library(ggh4x) # v 0.3.0; library(ggplot2) # v 3.5.2
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


nSubj <- 12 # number of participants
df <- subset(dat, N==nSubj)

df2 <- subset(df, Case==7 & B1==0.5)
sum(df2$BF10t > 3) / nrow(df2) # Under the alternative in Case 7, 43% of the Bayes factors exceed 3.

distBF <- distP <- ksBF <- ksP <- prop <- numeric(16); i <- 1
for (es in c(0, 0.5)) {
  temp <- subset(df, B1==es & Case==0)
  for (case in 0:7) {
    sub <- subset(df, B1==es & Case==case)
    prop[i] <- sum(sub$BF10t > 3) / nrow(sub) # proportion of the Bayes factors supporting the alternative
    distBF[i] <- ifelse(case==0, "EMD", sprintf("%.2f", wasserstein1d(log(temp$BF10t), log(sub$BF10t)))) # ln(BF₁₀)
    distP[i] <- ifelse(case==0, "EMD", round(wasserstein1d(temp$pVal, sub$pVal),10)) # earth mover's distance (EMD)
    ksBF[i] <- ifelse(case==0, "KS", round(ks.test(log(temp$BF10t), log(sub$BF10t))$statistic,2)) # ln(BF₁₀)
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

fig3 <- ggplot(df, aes(x=Case, y=log(BF10t), fill=Label)) + # ln(BF₁₀)
  geom_boxplot(alpha=.4) +
  facet_wrap(.~ES, scales="free_y", nrow=1) +
  labs(x="", y=expression(italic("BF")[10]), fill="Case",
       title=bquote("Paired" ~italic("t") * "-Test: " * italic("n") * "=" * .(nSubj))) +
  scale_fill_manual(values=c("#d83034", "#d83034", "#f9e858", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  facetted_pos_scales(y=list(scale_y_continuous(breaks=c(-log(10), -log(3), 0, log(3)),
                                                labels=c("1/10", "1/3", 1, 3)),
                             scale_y_continuous(breaks=c(-log(3), log(3), log(1E3)),
                                                labels=c("1/3", 3, expression(10^3))) )) +
  scale_x_continuous(breaks=seq(0,7,by=1)) +
  geom_hline(yintercept=c(-log(3), log(3)), linetype="dashed", color="gray") +
  geom_text(data=anno, aes(x=Case, y=Inf, label=distBF), vjust=53, col="red", fontface="bold") +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=20), legend.margin=margin(t=-15),
        legend.title=element_text(margin=margin(0,1,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm"))) +
  guides(fill=guide_legend(nrow=2, byrow=T))
ggsave("~/Downloads/osfstorage-archive/Appendix D/Figure D3.png", fig3, width=3800, height=2300, units="px", device="png", dpi=300)

