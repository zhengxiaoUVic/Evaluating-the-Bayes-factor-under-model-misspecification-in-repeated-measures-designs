library(dplyr) # v 1.1.4
library(ggh4x) # v 0.3.1; library(ggplot2) # v 4.0.1
library(ggpubr) # v 0.6.2

dat <- read.csv("~/Downloads/osfstorage-archive/Simulation 3/Figure 4/Sim_1t_robust.csv")
dat <- filter(dat, !if_any(contains("error"), ~ . > .05)) # 0 (out of 12,800) have been filtered out
df_med <- dat %>%
  group_by(N, Case, B1) %>%
  summarise(value=median(BF10t, na.rm=T), .groups="drop") %>%
  select(value, N, Case, B1)
df_iqr <- dat %>%
  group_by(N, Case, B1) %>%
  summarise(value=IQR(BF10t, na.rm=T), .groups="drop") %>%
  select(value, N, Case, B1)
df_med$Label <- df_iqr$Label <- factor(df_med$Case, labels=c("0: Normality holds",
                                                             "1: Normality holds",
                                                             "2: Uniform",
                                                             "3: Symmetric triangular",
                                                             "4: Symmetric bimodal",
                                                             "5: Asymmetric bimodal",
                                                             "6: Right-skewed",
                                                             "7: Left-skewed"))
df_med$ES <- df_iqr$ES <- factor(df_med$B1, labels=paste0(c("Null", "Positive Medium"), " ES"))


fig1 <- ggplot(df_med, aes(x=N, y=log(value))) +
  geom_line(aes(linetype=Label, color=Label)) +
  facet_grid2(~ES, scales="free_y", independent="y") +
  labs(x="Number of Participants", y=expression("Median  ln"~italic("BF")[10]), 
       color="Case", linetype="Case",
       title=expression("Paired"~italic("t") * "-Tests")) +
  scale_color_manual(values=c("#d83034", "#d83034", "black", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed", "dashed", "dotted", "dotted", "dotted")) +
  facetted_pos_scales(y=list(scale_y_continuous(breaks=c(-1.8, -1.5, -1.2)),
                             scale_y_continuous(breaks=c(2, 5, 8)))) +
  scale_x_continuous(breaks=c(12, 24, 48, 100)) +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=20), legend.margin=margin(t=-10),
        legend.title=element_text(margin=margin(0,0.5,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm")),
        legend.box.margin=margin(20,0,10,0),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.title.x=element_text(margin=margin(t=10)),
        legend.key.width=unit(1,"cm")) +
  guides(color=guide_legend(nrow=2, byrow=T),
         linetype=guide_legend(override.aes=list(linewidth=1)))

fig2 <- ggplot(df_iqr, aes(x=N, y=log(value))) +
  geom_line(aes(linetype=Label, color=Label)) +
  facet_grid2(~ES, scales="free_y", independent="y") +
  labs(x="Number of Participants", y=expression("IQR  ln"~italic("BF")[10]), 
       color="Case", linetype="Case", title="") +
  scale_color_manual(values=c("#d83034", "#d83034", "black", "#008dff", "#ff73b6", "#ff9d3a", "#003a7d", "#c701ff")) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed", "dashed", "dotted", "dotted", "dotted")) +
  facetted_pos_scales(y=list(scale_y_continuous(breaks=c(-2.4, -2, -1.6)),
                             scale_y_continuous(breaks=c(2, 7, 12)))) +
  scale_x_continuous(breaks=c(12, 24, 48, 100)) +
  theme_classic() +
  theme(legend.position="bottom", text=element_text(size=20), legend.margin=margin(t=-10),
        legend.title=element_text(margin=margin(0,0.5,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.5,0,0.5,"cm")),
        legend.box.margin=margin(20,0,10,0),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.title.x=element_text(margin=margin(t=10)),
        legend.key.width=unit(1,"cm")) +
  guides(color=guide_legend(nrow=2, byrow=T),
         linetype=guide_legend(override.aes=list(linewidth=1)))

ggarrange(fig1, fig2,
          nrow=2, ncol=1, labels=c("a", "b"), align="v",
          common.legend=T, legend="bottom",
          font.label=list(size=24))

ggsave("~/Downloads/osfstorage-archive/Appendix D/Figure D2.png", width=3800, height=3000, units="px", device="jpg", dpi=300)

