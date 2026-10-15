library(dplyr) # v 1.1.4
library(ggh4x) # v 0.3.0; library(ggplot2) # v 3.5.2
library(transport) # v 0.15-4

dat <- read.csv("~/Downloads/osfstorage-archive/Simulation 3 continued/Figure 5/Sim_1way_robust.csv")
dat <- filter(dat, !if_any(contains("error"), ~ . > .05)) # 0 (out of 25,600) have been filtered out
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

distBF <- distP <- ksBF <- ksP <- prop <- numeric(32); i <- 1
for (es in c(0, 0.25, 0.5, 0.75)) {
  temp <- subset(df, B1==es & Case==0)
  for (case in 0:7) {
    sub <- subset(df, B1==es & Case==case)
    prop[i] <- sum(sub$BF10 > 3) / nrow(sub) # proportion of the Bayes factors supporting the alternative
    distBF[i] <- ifelse(case==0, "EMD", sprintf("%.2f", wasserstein1d(log(temp$BF10), log(sub$BF10)))) # ln(BF₁₀)
    distP[i] <- ifelse(case==0, "EMD", round(wasserstein1d(temp$pVal, sub$pVal),10)) # earth mover's distance (EMD)
    ksBF[i] <- ifelse(case==0, "KS", round(ks.test(log(temp$BF10), log(sub$BF10))$statistic,2)) # ln(BF₁₀)
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

fig4 <- ggplot(subset(df, ES %in% c("Null ES", "Positive Medium ES")), aes(x=Case, y=log(BF10), fill=Label)) + # ln(BF₁₀)
  geom_boxplot(alpha=.4) +
  facet_wrap(.~ES, scales="free_y", nrow=1) +
  labs(x="", y=expression(italic("BF")[10]), fill="Case",
       title=bquote("One-Way Repeated-Measures Designs: " * italic("n") * "=" * .(nSubj))) +
  scale_fill_manual(values=c("#d83034", "#d83034", "#f9e858", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  facetted_pos_scales(y=list(scale_y_continuous(breaks=c(-log(3), 0, log(3), log(10)),
                                                labels=c("1/3", 1, 3, 10)),
                             scale_y_continuous(breaks=c(-log(3), log(3), log(1E2), log(1E5)),
                                                labels=c("1/3", 3, expression(10^2), expression(10^5))) )) +
  scale_x_continuous(breaks=seq(0,7,by=1)) +
  geom_hline(yintercept=c(-log(3), log(3)), linetype="dashed", color="gray") +
  geom_text(data=subset(anno, ES %in% c("Null ES", "Positive Medium ES")), 
            aes(x=Case, y=Inf, label=distBF), vjust=49, col="red", fontface="bold") +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=20), legend.margin=margin(t=-15),
        legend.title=element_text(margin=margin(0,1,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm"))) +
  guides(fill=guide_legend(nrow=3, byrow=T))
ggsave("~/Downloads/osfstorage-archive/Simulation 3 continued/Figure 5/Figure 5.png", fig4, width=3800, height=2300, units="px", device="png", dpi=300)

