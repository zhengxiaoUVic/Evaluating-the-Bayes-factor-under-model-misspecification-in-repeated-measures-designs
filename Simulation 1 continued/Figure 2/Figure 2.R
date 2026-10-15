library(dplyr) # v 1.1.4
library(ggpubr) # v 0.6.0; library(ggplot2) # v 3.5.2

filepath <- c("~/Downloads/osfstorage-archive/Simulation 1 continued/data/Het3_N48_n1_ES0.5.csv",
              "~/Downloads/osfstorage-archive/Simulation 1 continued/data/CS_N48_n1_ES0.5_rho0.13.csv",
              "~/Downloads/osfstorage-archive/Simulation 1 continued/data/CS_N48_n1_ES0.5_rho-0.3.csv",
              "~/Downloads/osfstorage-archive/Simulation 1 continued/data/CS_N48_s3_ES0.5_rho-0.3.csv")
fig <- vector(mode="list", length=4)

for (index in 1:4) {
  df <- read.csv(filepath[index])
  
  if (index==2) {
    title <- bquote(italic("n") * "=48, " * italic("R") * "=1, CS Holds*, " * italic("ρ") * "=.13")
  } else if (index==3) {
    title <- bquote(italic("n") * "=48, " * italic("R") * "=1, CS Holds*, " * italic("ρ") * "=-.3")
  } else if (index==4) {
    title <- bquote(italic("n") * "=48, " * italic("R") * "=3, CS Holds*, " * italic("ρ") * "=-.3")
  } else {
    title <- bquote(italic("n") * "=48, " * italic("R") * "=1, CS is Violated*")
  }
  
  fig[[index]] <- ggplot(filter(df, !if_any(contains("errorA"), ~ . > .05)), # Monte Carlo errors
                         aes(-log(legacy01_A), -log(BF01_A), # ln BF₁₀ (`A` is the SECOND factor in our paper)
                             color=factor(cut(pValA, breaks=c(-1E-5, .01, .05, .1, 1),
                                              labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
                             shape=factor(cut(pValA, breaks=c(-1E-5, .01, .05, .1, 1),
                                              labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
    geom_point(size=2, alpha=.6) +
    scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
    labs(x=expression("RIO "~italic("BF")[10]),
         y=expression("MRE "~italic("BF")[10]),
         color=expression(paste(italic("p-"), "value")),
         shape=expression(paste(italic("p-"), "value")),
         title=title) +
    scale_x_continuous(limits=c(-3, 15), expand=c(0, 0),
                       breaks=c(-log(3), 0, log(3), log(100), log(1E5)),
                       labels=c("1/3  ", 1, 3, expression(10^2), expression(10^5))) + # white space left intentionally
    scale_y_continuous(limits=c(-3, 15), expand=c(0, 0),
                       breaks=c(-log(3), 0, log(3), log(100), log(1E5)),
                       labels=c("1/3", 1, 3, expression(10^2), expression(10^5))) +
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
}

fig2 <- ggarrange(fig[[1]], fig[[2]], fig[[3]], fig[[4]],
                  nrow=2, ncol=2, labels=c("a", "b", "c", "d"), align="h",
                  common.legend=T, legend="bottom", font.label=list(size=23))
ggsave("~/Downloads/osfstorage-archive/Simulation 1 continued/Figure 2/Figure 2.png", fig2, width=3000, height=3000, units="px", device="jpg", dpi=300)





df <- read.csv("~/Downloads/osfstorage-archive/Simulation 1 continued/data/CS_N48_s3_ES0_rho-0.3.csv")
# df <- read.csv("~/Downloads/osfstorage-archive/Simulation 1 continued/data/CS_N48_n1_ES0_rho-0.3.csv")
lapply(df[, grepl("error", names(df))], quantile, na.rm=T) # check the maximum Monte Carlo error rates
BF01_A <- subset(df, BF_errorA <= .05)$BF01_A
legacy01_A <- subset(df, legacy_errorA <= .05)$legacy01_A

sum(BF01_A < 1/3) / length(BF01_A)
sum(legacy01_A < 1/3) / length(legacy01_A) # explore the type I error rate

