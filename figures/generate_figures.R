
# Filename: generate_figures.R
# Version: 2024.04.2+764
# Author: vani taluja
# Date Created: 02 july 2024
# Last Modified: 02 july 2024
# Description: tidy script for generating figures. all legends (i.e. color, shape) were added externally.


# load libraries (and install, if needed)
# install.packages("ggplot2")
# install.packages("patchwork")
library(ggplot2)
library(patchwork)

# read in data
# note: in order to read in correctly,
# all data files must be stored in the same folder as this script
outcome_avgs = read.csv("byCluster_asdoutcome.csv")
mixed_avgs = read.csv("byCluster_mixed.csv")
mixed_avgs_bydx = read.csv("byDx_mixed.csv")
nhb_avgs = read.csv("byCluster_nhb.csv")
cv = read.csv("crossvalidation.csv")
nmi = read.csv("nmi.csv")

# set consistent theme and colors for plotting
theme_set(new = theme_classic() + 
            theme(axis.title = element_text(size = 14),
                  axis.text = element_text(size = 12)))
outcome_colors = c("#BFD3E6","#8C96C6","#465892")
mixed_colors = c("#ff6a6a", "#cd5555", "#8b3a3a")
nhb_colors = c("#A1D99B", "#74C476", "#006D2C")


# Figure 1 ----

# Figure 1a
vec = data.frame(v=c(1.5,2.5,3.5))
ggplot(vec) + geom_point(aes(x=v,y=v, color=factor(v)), shape = 15, size = 30) + 
  scale_color_manual(values=outcome_colors[3:1]) +
  geom_text(aes(x=v, y=v, label=c("Outcome-3", "Outcome-2", "Outcome-1")), size=4) +
  geom_text(aes(x=v, y=v+0.7, label=c("n = 16", "n = 30", "n = 35")), size=3.5) +
  xlim(1,4) + ylim(1,4.3) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"))),
        axis.ticks = element_blank(),
        aspect.ratio = 1) +
  labs(x="Temporal Cortex Activation", y="Clinical Scores",
       title = "A. Hypothesized Relationship")

# Figure 1b
lh_el = ggplot(outcome_avgs) + 
  geom_point(aes(x=LHtemporal, y=el, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=el-el_sd, ymax=el+el_sd), width=0.002) +
  xlim(-0.015,0.05) + ylim(20,120) +
  labs(y="Mullen EL Ratio") +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_el = ggplot(outcome_avgs) + 
  geom_point(aes(x=RHtemporal, y=el, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=el-el_sd, ymax=el+el_sd), width=0.002) +
  xlim(-0.01,0.05) + ylim(20, 120) +
  labs(y="Mullen EL Ratio") +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        legend.position ="none")

lh_rl = ggplot(outcome_avgs) + 
  geom_point(aes(x=LHtemporal, y=rl, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=rl-rl_sd, ymax=rl+rl_sd), width=0.002) +
  xlim(-0.015,0.05) + ylim(33,110) +
  labs(y="Mullen RL Ratio") +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_rl = ggplot(outcome_avgs) + 
  geom_point(aes(x=RHtemporal, y=rl, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=rl-rl_sd, ymax=rl+rl_sd), width=0.002) +
  xlim(-0.01,0.05) + ylim(33,110) +
  labs(y="Mullen RL Ratio") +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        legend.position ="none")

lh_adap = ggplot(outcome_avgs) + 
  geom_point(aes(x=LHtemporal, y=adap, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=adap-adap_sd, ymax=adap+adap_sd), width=0.002) +
  xlim(-0.015,0.05) + ylim(50,105) +
  labs(y="Vineland ABC") +
  theme_classic() +
  theme(legend.position ="none")

rh_adap = ggplot(outcome_avgs) + 
  geom_point(aes(x=RHtemporal, y=adap, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=adap-adap_sd, ymax=adap+adap_sd), width=0.002) +
  xlim(-0.01,0.05) + ylim(50,105) +
  labs(y="Vineland ABC") +
  theme_classic() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position ="none")

lh_el + rh_el + lh_rl + rh_el +
  lh_adap + rh_adap + plot_layout(ncol=2) +
  plot_annotation(title = "B. SNF of Clinical and fMRI Measures") & 
  theme(aspect.ratio = 1.3, title = element_text(size = 9))

# Figure 2 ----
# Figure 2a
lh_elc = ggplot(outcome_avgs) + 
  geom_point(aes(x=LHtemporal, y=elc, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=elc-elc_sd, ymax=elc+elc_sd), width=0.002) +
  xlim(-0.015,0.05) + ylim(45,111) +
  labs(y="Mullen ELC") +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_elc = ggplot(outcome_avgs) + 
  geom_point(aes(x=RHtemporal, y=elc, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=elc-elc_sd, ymax=elc+elc_sd), width=0.002) +
  xlim(-0.01,0.05) + ylim(45,111) +
  labs(y="Mullen ELC") +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position ="none")

lh_et = ggplot(outcome_avgs) + 
  geom_point(aes(x=LHtemporal, y=data.perc.soc, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=data.perc.soc-data.perc.soc_sd, 
                    ymax=data.perc.soc+data.perc.soc_sd), width=0.002) +
  xlim(-0.015,0.05) + ylim(15,87) +
  labs(y="% Fixation Social") +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_et = ggplot(outcome_avgs) + 
  geom_point(aes(x=RHtemporal, y=data.perc.soc, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=data.perc.soc-data.perc.soc_sd, ymax=data.perc.soc+data.perc.soc_sd), width=0.002) +
  xlim(-0.01,0.05) + ylim(15,87) +
  labs(y="% Fixation Social") +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        legend.position ="none")

lh_ados = ggplot(outcome_avgs) + 
  geom_point(aes(x=LHtemporal, y=tot, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=tot-tot_sd, ymax=tot+tot_sd), width=0.002) +
  xlim(-0.015,0.05) + ylim(11,27) +
  labs(y="ADOS Total Score") +
  theme_classic() +
  theme(legend.position ="none")

rh_ados = ggplot(outcome_avgs) + 
  geom_point(aes(x=RHtemporal, y=tot, color=factor(Cluster.last)), size=5) + 
  scale_color_manual(values=outcome_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=tot-tot_sd, ymax=tot+tot_sd), width=0.002) +
  xlim(-0.01,0.05) + ylim(11,27) +
  labs(y="ADOS Total Score") +
  theme_classic() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position ="none")

lh_elc + rh_elc + lh_et + rh_et +
  lh_ados + rh_ados + 
  plot_annotation(title = "A. External Validation") +
  plot_layout(ncol=2) & theme(aspect.ratio = 1.3)

# Figure 2b
theme_set(new = theme_classic())

cv_outcome = ggplot(cv, aes(x=1,y=outcome)) +
  geom_violin(fill = "#8C96C6", alpha = 0.6) +
  geom_boxplot(width = 0.25) +
  geom_point(position = position_jitter(width=0.2, height=0.05, seed = 7)) +
  stat_summary(geom = "point",
               fun = "mean",
               size = 4) +
  labs(y = "Accuracy",
       title = "B. Cross Validation") +
  theme_classic() +
  ylim(0.57,1.03) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

# Figure 2c
nmi_outcome = ggplot(nmi) +
  stat_boxplot(aes(x=factor(percentsampleremoved), y=outcome_index), 
               geom = "errorbar", width=0.25) +
  geom_boxplot(aes(x=factor(percentsampleremoved), y=outcome_index, 
                   group=percentsampleremoved), 
               fill = "#8C96C6", alpha = 0.6) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  scale_x_discrete(labels = c("5%", "10%", "20%", "30%", "40%", "50%")) +
  labs(x="Sample Randomly Removed",
       y="NMI Index",
       title="C. NMI") +
  theme(legend.position = "none")


# Figure 2d
annotation = data.frame(x = c(rep(c("1", "2", "3"), each=3)),
                        y = c(rowMeans(outcome_avgs[c("rl", "rl.1")]),
                              rowMeans(outcome_avgs[c("el", "el.1")]),
                              rowMeans(outcome_avgs[c("adap", "adap.1")])),
                        label = c("p<0.001*","p=0.685","p=0.556",
                                  "p=0.062","p=0.005*","p<0.001*",
                                  "p<0.001*","p=0.011*","p=0.495"))

clinchange_outcome = ggplot(outcome_avgs) + 
  geom_point(aes(x="1", y=rl, color=factor(Cluster.last)), size=7, shape=18) +
  geom_point(aes(x="1", y=rl.1, color=factor(Cluster.last)), size=4) + 
  geom_point(aes(x="2", y=el, color=factor(Cluster.last)), size=7, shape=18) +
  geom_point(aes(x="2", y=el.1, color=factor(Cluster.last)), size=4) +  
  geom_point(aes(x="3", y=adap, color=factor(Cluster.last)), size=7, shape=18) +
  geom_point(aes(x="3", y=adap.1, color=factor(Cluster.last)), size=4) +
  scale_color_manual(values=outcome_colors) +
  geom_segment(aes(x="1", xend="1", y=rl.1, yend=rl), 
               arrow = arrow(length = unit(0.15, "inches"))) +
  geom_segment(aes(x="2", xend="2", y=el.1, yend=el), 
               arrow = arrow(length = unit(0.15, "inches"))) +
  geom_segment(aes(x="3", xend="3", y=adap.1, yend=adap), 
               arrow = arrow(length = unit(0.15, "inches"))) +
  labs(y="Standardized Score", 
       title = "D. Longitudinal Clinical Change") +
  scale_x_discrete(labels=c("1" = "Receptive\nLanguage\nRatio Score",
                            "2" = "Expressive\nLanguage\nRatio Score", 
                            "3" = "Adaptive\nBehavior\nComposite")) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') + 
  geom_text(data=annotation, aes(x=x, y=y, label=label),
            angle=90, vjust = -1.8, size = 3)

# print Figure 2b - d
design <- "AAACCCC
           BBBCCCC"

cv_outcome + free(nmi_outcome) + 
  clinchange_outcome +
  plot_layout(design = design)

# Figure 4 ----
# Figure 4a
theme_set(new = theme_classic() + 
            theme(axis.title = element_text(size = 12),
                  axis.text = element_text(size = 10)))

lh_el_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=LHtemporal, y=el, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=el-el_sd, ymax=el+el_sd), width=0.002) +
  xlim(0,0.055) + ylim(30,130) +
  labs(y="Mullen EL Ratio") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_el_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=RHtemporal, y=el, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=el-el_sd, ymax=el+el_sd), width=0.002) +
  xlim(0,0.07) + ylim(30,130) +
  labs(y="Mullen EL Ratio") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position ="none")

lh_rl_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=LHtemporal, y=rl, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=rl-rl_sd, ymax=rl+rl_sd), width=0.002) +
  xlim(0,0.055) + ylim(40,130) +
  labs(y="Mullen RL Ratio") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_rl_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=RHtemporal, y=rl, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=rl-rl_sd, ymax=rl+rl_sd), width=0.002) +
  xlim(0,0.07) + ylim(40,130) +
  labs(y="Mullen RL Ratio") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position ="none")

lh_adap_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=LHtemporal, y=adap, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=adap-adap_sd, ymax=adap+adap_sd), width=0.002) +
  xlim(0,0.055) + ylim(60,125) +
  labs(y="Vineland ABC") +
  theme(legend.position ="none")

rh_adap_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=RHtemporal, y=adap, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=adap-adap_sd, ymax=adap+adap_sd), width=0.002) +
  xlim(0,0.07) + ylim(60,125) +
  labs(y="Vineland ABC") +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position ="none")

# Figure 4b
lh_elc_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=LHtemporal, y=elc, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=elc-elc_sd, ymax=elc+elc_sd), width=0.002) +
  xlim(0,0.055) + ylim(45,135) +
  labs(y="Mullen ELC") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_elc_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=RHtemporal, y=elc, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=elc-elc_sd, ymax=elc+elc_sd), width=0.002) +
  xlim(0,0.07) + ylim(45,135) +
  labs(y="Mullen ELC") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position ="none")

lh_et_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=LHtemporal, y=data.perc.soc, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=data.perc.soc-data.perc.soc_sd, 
                    ymax=data.perc.soc+data.perc.soc_sd), width=0.002) +
  xlim(0,0.055) + ylim(20,95) +
  labs(y="% Fixation Social") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_et_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=RHtemporal, y=data.perc.soc, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=data.perc.soc-data.perc.soc_sd, 
                    ymax=data.perc.soc+data.perc.soc_sd), width=0.002) +
  xlim(0,0.07) + ylim(20,95) +
  labs(y="% Fixation Social") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position ="none")

lh_ados_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=LHtemporal, y=tot, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=tot-tot_sd, ymax=tot+tot_sd), width=0.002) +
  xlim(0,0.055) + ylim(-2,25) +
  labs(y="ADOS Total Score") +
  theme(legend.position ="none")

rh_ados_mix = ggplot(mixed_avgs) + 
  geom_point(aes(x=RHtemporal, y=tot, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=mixed_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=tot-tot_sd, ymax=tot+tot_sd), width=0.002) +
  xlim(0,0.07) + ylim(-2,25) +
  labs(y="ADOS Total Score") +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position ="none")

# Figure 4c
change_asd_mixed = ggplot(mixed_avgs_bydx[which(mixed_avgs_bydx$dx == "ASD"),]) +
  geom_abline(slope = 0, intercept = 0) +
  geom_point(aes(x="1", y=RLT_z, color=factor(Cluster.last)), size=3) +
  geom_point(aes(x="2", y=ELT_z, color=factor(Cluster.last)), size=3) + 
  geom_point(aes(x="3", y=ABC_z, color=factor(Cluster.last)), size=3) +
  scale_color_manual(values=mixed_colors) +
  labs(title = "ASD", y="Standardized Difference") +
  scale_x_discrete(labels=c("1" = "RL", 
                            "2" = "EL", 
                            "3" = "ABC")) +
  theme(axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none")

change_td_mixed = ggplot(mixed_avgs_bydx[which(mixed_avgs_bydx$dx == "TD"),]) + 
  geom_abline(slope = 0, intercept = 0) +
  geom_point(aes(x="1", y=RLT_z, color=factor(Cluster.last)), size=3) +
  geom_point(aes(x="2", y=ELT_z, color=factor(Cluster.last)), size=3) + 
  geom_point(aes(x="3", y=ABC_z, color=factor(Cluster.last)), size=3) +
  scale_color_manual(values=mixed_colors) +
  labs(title = "TD", y="Standardized Difference") +
  scale_x_discrete(labels=c("1" = "RL", 
                            "2" = "EL", 
                            "3" = "ABC")) +
  theme(axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none")

change_delay_mixed = ggplot(mixed_avgs_bydx[which(mixed_avgs_bydx$dx == "Delay"),]) + 
  geom_abline(slope = 0, intercept = 0) +
  geom_point(aes(x="1", y=RLT_z, color=factor(Cluster.last)), size=3) +
  geom_point(aes(x="2", y=ELT_z, color=factor(Cluster.last)), size=3) + 
  geom_point(aes(x="3", y=ABC_z, color=factor(Cluster.last)), size=3) +
  scale_color_manual(values=mixed_colors) +
  labs(title = "Delay", y="Standardized Difference") +
  scale_x_discrete(labels=c("1" = "RL", 
                            "2" = "EL", 
                            "3" = "ABC")) +
  theme(axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none")

# plotting figure 4a - c
fig4a = lh_el_mix + rh_el_mix +
  lh_rl_mix + rh_rl_mix +
  lh_adap_mix + rh_adap_mix +
  plot_layout(ncol=2) +
  plot_annotation(title = "A. SNF of Clinical and fMRI Measures") & 
  theme(aspect.ratio = 1.2, plot.title = element_text(size = 12)) 

fig4b = lh_elc_mix + rh_elc_mix +
  lh_et_mix + rh_et_mix +
  lh_ados_mix + rh_ados_mix +
  plot_layout(ncol=2) +
  plot_annotation(title = "B. External Validation") &
  theme(aspect.ratio = 1.2, plot.title = element_text(size = 12,
                                                      hjust = 0.25))

fig4c = change_asd_mixed + (change_td_mixed / change_delay_mixed) +
  plot_annotation(title = "C. Standardized Clinical Change by Diagnostic Group") &
  theme(aspect.ratio = 2.5, plot.title = element_text(size = 12))

wrap_elements(fig4a) + plot_spacer() + wrap_elements(fig4b) + 
  plot_spacer() + wrap_elements(fig4c) +
  plot_layout(widths = c(2, -0.35, 2, -0.25, 3))

# Figure 5c ----
lh_el_nhb = ggplot(nhb_avgs) + 
  geom_point(aes(x=LHtemporal, y=el, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=nhb_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=el-el_sd, ymax=el+el_sd), width=0.002) +
  xlim(0.02,0.05) + ylim(45,130) +
  labs(y="Mullen EL Ratio") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_el_nhb = ggplot(nhb_avgs) + 
  geom_point(aes(x=RHtemporal, y=el, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=nhb_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=el-el_sd, ymax=el+el_sd), width=0.002) +
  xlim(0.02,0.08) + ylim(45,130) +
  labs(y="Mullen EL Ratio") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position ="none")

lh_rl_nhb = ggplot(nhb_avgs) + 
  geom_point(aes(x=LHtemporal, y=rl, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=nhb_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=rl-rl_sd, ymax=rl+rl_sd), width=0.002) +
  xlim(0.02,0.05) + ylim(50,120) +
  labs(y="Mullen RL Ratio") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position ="none")

rh_rl_nhb = ggplot(nhb_avgs) + 
  geom_point(aes(x=RHtemporal, y=rl, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=nhb_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=rl-rl_sd, ymax=rl+rl_sd), width=0.002) +
  xlim(0.02,0.08) + ylim(50,120) +
  labs(y="Mullen RL Ratio") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position ="none")

lh_adap_nhb = ggplot(nhb_avgs) + 
  geom_point(aes(x=LHtemporal, y=adap, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=nhb_colors) +
  geom_errorbar(aes(x=LHtemporal, ymin=adap-adap_sd, ymax=adap+adap_sd), width=0.002) +
  xlim(0.02,0.05) + ylim(65,115) +
  labs(y="Vineland ABC") +
  theme(legend.position ="none")

rh_adap_nhb = ggplot(nhb_avgs) + 
  geom_point(aes(x=RHtemporal, y=adap, color=factor(Cluster.last)), size=4) + 
  scale_color_manual(values=nhb_colors) +
  geom_errorbar(aes(x=RHtemporal, ymin=adap-adap_sd, ymax=adap+adap_sd), width=0.002) +
  xlim(0.02,0.08) + ylim(65,115) +
  labs(y="Vineland ABC") +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position ="none")

lh_el_nhb + rh_el_nhb +
  lh_rl_nhb + rh_rl_nhb +
  lh_adap_nhb + rh_adap_nhb +
  plot_layout(ncol=2) &
  theme(aspect.ratio = 1)
