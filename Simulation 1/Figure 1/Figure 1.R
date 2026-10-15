library(dplyr) # v 1.1.4
library(ggpubr) # v 0.6.0; library(ggplot2) # v 3.5.2

fig <- vector(mode="list", length=4)
index <- 1

for (N in c(48, 100)) { # number of participants
  for (isCS in c(F, T)) { # whether compound symmetry holds
    df_eff <- read.csv(paste0("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N",
                              N,"_eqVarTRUE_eqCov",isCS,"_nullFALSE.csv")) # simulated under an effect
    df_null <- read.csv(paste0("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N",
                               N,"_eqVarTRUE_eqCov",isCS,"_nullTRUE.csv")) # simulated under the null
    df <- rbind(df_eff, df_null) # simulated data combined by the likelihood principle
    
    if (index %in% c(1,3)) {
      title <- bquote(italic("n") * "=" * .(N) * ", " * italic("R") * "=1, CS is Violated")
    } else {
      title <- bquote(italic("n") * "=" * .(N) * ", " * italic("R") * "=1, CS Holds, " * italic("ρ") * "=.6")
    }

    fig[[index]] <- ggplot(filter(df, !if_any(contains("errorB"), ~ . > .05)), # the maximum Monte Carlo error rate is 2.6%
                           aes(-log(legacy01_B), -log(BF01_B), # ln BF₁₀
                               color=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                                                labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
                               shape=factor(cut(pValB,breaks=c(-1E-5, .01, .05, .1, 1),
                                                labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
      geom_point(size=2, alpha=.5) +
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
    
    index <- index + 1
  }
}

fig1 <- ggarrange(fig[[2]], fig[[4]], fig[[1]], fig[[3]],
                  nrow=2, ncol=2, labels=c("a", "b", "c", "d"), align="h",
                  common.legend=T, legend="bottom", font.label=list(size=23))
ggsave("~/Downloads/osfstorage-archive/Simulation 1/Figure 1/Figure 1.png", fig1, width=3000, height=3000, units="px", device="jpg", dpi=300)





################################################################################ colored by the effect size

fig <- vector(mode="list", length=4)
index <- 1

for (N in c(48, 100)) { # number of participants
  for (isCS in c(F, T)) { # whether compound symmetry holds
    df_eff <- read.csv(paste0("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N",
                              N,"_eqVarTRUE_eqCov",isCS,"_nullFALSE.csv")) # simulated under an effect
    df_eff$ES <- "true"
    df_null <- read.csv(paste0("~/Downloads/osfstorage-archive/Simulation 1/data/Para_N",
                               N,"_eqVarTRUE_eqCov",isCS,"_nullTRUE.csv")) # simulated under the null
    df_null$ES <- "false"
    df <- rbind(df_eff, df_null) # simulated data combined by the likelihood principle
    df$ES <- factor(df$ES)
    
    if (index %in% c(1,3)) {
      title <- bquote(italic("n") * "=" * .(N) * ", " * italic("R") * "=1, CS is Violated")
    } else {
      title <- bquote(italic("n") * "=" * .(N) * ", " * italic("R") * "=1, CS Holds, " * italic("ρ") * "=.6")
    }
    
    fig[[index]] <- ggplot(filter(df, !if_any(contains("errorB"), ~ . > .05)), # the maximum Monte Carlo error rate is 2.6%
                           aes(-log(legacy01_B), -log(BF01_B), # ln BF₁₀
                               color=ES, shape=ES)) +
      geom_point(size=2, alpha=.5) +
      scale_color_manual(values=c("black", "#ff7411")) +
      labs(x=expression("RIO "~italic("BF")[10]),
           y=expression("MRE "~italic("BF")[10]),
           color="The alternative hypothesis is ",
           shape="The alternative hypothesis is ",
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
    
    index <- index + 1
  }
}

figH <- ggarrange(fig[[2]], fig[[4]], fig[[1]], fig[[3]],
                  nrow=2, ncol=2, labels=c("a", "b", "c", "d"), align="h",
                  common.legend=T, legend="bottom", font.label=list(size=23))
ggsave("~/Downloads/osfstorage-archive/Simulation 1/Figure 1/Figure 1S.png", figH, width=3000, height=3000, units="px", device="jpg", dpi=300)

