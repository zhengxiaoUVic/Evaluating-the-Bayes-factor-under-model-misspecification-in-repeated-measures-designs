library(dplyr) # v 1.1.4
library(ggh4x) # v 0.3.1; library(ggplot2) # v 4.0.1
library(ggpubr) # v 0.6.2

dat <- read.csv("~/Downloads/osfstorage-archive/Simulation 2/Figure 3/Sim2x2_robust.csv")
dat <- filter(dat, !if_any(contains("error"), ~ . > .05)) # 53 (out of 57,600) have been filtered out
df_med <- dat %>%
  group_by(N, n, Case, B1) %>%
  summarise(value=median(BF01_A, na.rm=T), .groups="drop") %>%
  select(value, N, n, Case, B1)
df_iqr <- dat %>%
  group_by(N, n, Case, B1) %>%
  summarise(value=IQR(BF01_A, na.rm=T), .groups="drop") %>%
  select(value, N, n, Case, B1)
df_med$Label <- df_iqr$Label <- factor(df_med$Case, labels=c("0: Normality holds",
                                                             "1: Normality holds",
                                                             "2: Random effects are right-skewed",
                                                             "3: Residuals are right-skewed",
                                                             "4: Both are right-skewed",
                                                             "5: Random effects are left-skewed",
                                                             "6: Residuals are left-skewed",
                                                             "7: Both are left-skewed"))
df_med$ES <- df_iqr$ES <- factor(df_med$B1, labels=paste0(c("Null", paste0("Positive ", c("Small", "Medium", "Large"))), " ES"))
df_med$n <- df_iqr$n <- factor(df_med$n, labels=c("1 Trial", "4 Trials", "10 Trials"))


fig1 <- ggplot(df_med, aes(x=N, y=-log(value))) +
  geom_line(aes(linetype=Label, color=Label)) +
  facet_grid2(n~ES, scales="free_y", independent="y") +
  labs(x="Number of Participants", y=expression("Median  MRE  ln"~italic("BF")[10]), 
       color="Case", linetype="Case",
       title="2×2 Repeated-Measures Designs") +
  scale_color_manual(values=c("#d83034", "#d83034", "black", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed", "dashed", "dotted", "dotted", "dotted")) +
  facetted_pos_scales(y=list(scale_y_continuous(breaks=c(-1.3, -1.1, -0.9)), # R1C1
                             scale_y_continuous(breaks=c(-0.8, -0.7, -0.6)), # R1C2
                             scale_y_continuous(breaks=c(0, 1, 2)), # R1C3
                             scale_y_continuous(breaks=c(2, 4, 6)), # R1C4
                             scale_y_continuous(breaks=c(-1.8, -1.6, -1.4)), # R2C1
                             scale_y_continuous(breaks=c(-1, -0.3, 0.4)), # R2C2
                             scale_y_continuous(breaks=c(0, 3, 6)), # R2C3
                             scale_y_continuous(breaks=c(4, 8, 12)), # R2C4
                             scale_y_continuous(breaks=c(-2, -1.7, -1.4)), # R3C1
                             scale_y_continuous(breaks=c(-1, 0, 1)), # R3C2
                             scale_y_continuous(breaks=c(2, 6, 10)), # R3C3
                             scale_y_continuous(breaks=c(5, 10, 15)))) +
  scale_x_continuous(breaks=c(12, 24, 48)) +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=20), legend.margin=margin(t=-15),
        legend.title=element_text(margin=margin(0,1,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm")),
        legend.box.margin=margin(20,0,0,0),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.title.x=element_text(margin=margin(t=10)),
        legend.key.width=unit(1,"cm")) +
  guides(color=guide_legend(nrow=3, byrow=T),
         linetype=guide_legend(override.aes=list(linewidth=1)))

fig2 <- ggplot(df_iqr, aes(x=N, y=-log(value))) +
  geom_line(aes(linetype=Label, color=Label)) +
  facet_grid2(n~ES, scales="free_y", independent="y") +
  labs(x="Number of Participants", y=expression("IQR  MRE  ln"~italic("BF")[10]), 
       color="Case", linetype="Case",
       title="") +
  scale_color_manual(values=c("#d83034", "#d83034", "black", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed", "dashed", "dotted", "dotted", "dotted")) +
  facetted_pos_scales(y=list(scale_y_continuous(breaks=c(-0.3, 0, 0.3)), # R1C1
                             scale_y_continuous(breaks=c(-0.8, -0.4, 0)), # R1C2
                             scale_y_continuous(breaks=c(0, 0.5, 1)), # R1C3
                             scale_y_continuous(breaks=c(0, 2, 4)), # R1C4
                             scale_y_continuous(breaks=c(-1, -0.5, 0)), # R2C1
                             scale_y_continuous(breaks=c(-1, -0.8, -0.6)), # R2C2
                             scale_y_continuous(breaks=c(0, 2, 4)), # R2C3
                             scale_y_continuous(breaks=c(4, 8, 12)), # R2C4
                             scale_y_continuous(breaks=c(-1.4, -1, -0.6)), # R3C1
                             scale_y_continuous(breaks=c(-0.9, -0.5, -0.1)), # R3C2
                             scale_y_continuous(breaks=c(2, 4, 6)), # R3C3
                             scale_y_continuous(breaks=c(5, 10, 15)))) +
  scale_x_continuous(breaks=c(12, 24, 48)) +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=20), legend.margin=margin(t=-15),
        legend.title=element_text(margin=margin(0,1,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm")),
        legend.box.margin=margin(20,0,0,0),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.title.x=element_text(margin=margin(t=10)),
        legend.key.width=unit(1,"cm")) +
  guides(color=guide_legend(nrow=3, byrow=T),
         linetype=guide_legend(override.aes=list(linewidth=1)))

ggarrange(fig1, fig2,
          nrow=2, ncol=1, labels=c("a", "b"), align="v",
          common.legend=T, legend="bottom",
          font.label=list(size=24))

ggsave("~/Downloads/osfstorage-archive/Appendix D/Figure D1.png", width=3800, height=4500, units="px", device="jpg", dpi=300)

