#---
#title: "T-cell manuscript - Supporting Material"
#author: "Mercè Garí"
#date: '2021-03-25'
#---

# Load packages

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(colorspace)
library(readxl)
library(forcats)
library(scales)
library(gridExtra)
library(cowplot)
library(Cairo)
library(stringr)
library(ggpubr)
library(rstatix)
library(EnvStats)
library(lme4)
library(effects)
library(broom)
library(broom.mixed)

library(extrafont)
loadfonts(quiet = TRUE)
library(ggthemes)

# Define theme and colors
theme_set(theme_classic())
theme_update(strip.background = element_rect(colour = "white", fill = "white"))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

here_r = function (...) here::here("Statistics", ...)
here_data_tcell = function (...) here::here("Data-t-cell", ...)
here_output = function (...) here::here("Output", ...) 

# Load data of t-cell
load(here_data_tcell("data-raw-tcell-210315.RData"))

data <- data %>%
  rename(NC = "Capsid sub x5",
         M = "Membrane sub x5",
         SNT = "SNT sub x5",
         SCT = "SCT sub x5",
         `Pepmix (green)` = "grün sbu x5",
         `Pepmix (yellow)` = "yellow sub x5")
data.long <- data.long %>%
  mutate(Antigen = case_when(
    Antigen == "Capsid sub x5" ~ "NC",
    Antigen == "Membrane sub x5" ~ "M",
    Antigen == "SNT sub x5" ~ "SNT",
    Antigen == "SCT sub x5" ~ "SCT",
    Antigen == "grün sbu x5" ~ "Pepmix (green)",
    Antigen == "yellow sub x5" ~ "Pepmix (yellow)"))

data.long <- data.long %>%
  mutate(IgG = as.numeric(IgG))

data <- data %>%
  mutate(IgG = as.numeric(IgG))

#------------------------------------------------ Figure S1

# How is IFGn for this group in relation to IgG +/-

wt <- data.long %>%
  mutate(value = ifelse(value == 0, 0.01, value)) %>%
  filter(Group.5 %in% "PCR+ Seronegative") %>%
 # filter(Antigen %in% "NC") %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  filter(!is.na(IgG)) %>%
  select(LabId, Seropositivity.IgG, 
         IgG, Antigen, value) 
head(wt)
#wilcox.test(wt$value ~ wt$Seropositivity.IgG)

# NON BONFERRONI CORRECTION
stat.test <- wt %>%
  group_by(Antigen) %>%
  wilcox_test(value ~ Seropositivity.IgG) %>%
 # adjust_pvalue(method = "bonferroni") %>%
  add_significance()
data.frame(stat.test)

stat.test <- stat.test %>% add_xy_position(x = "Seropositivity.IgG")
stat.test
data.long %>%
  mutate(value = ifelse(value == 0, 0.01, value)) %>%
  filter(Group.5 %in% "PCR+ Seronegative") %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  filter(!is.na(IgG)) %>%
  ggplot(aes(x=Seropositivity.IgG, y=value)) +
  geom_jitter(width=0.1, alpha=0.4, color=cbPalette[8]) +
  facet_wrap(~(Antigen), ncol=4) +
  geom_hline(yintercept = 40, col="grey", lty=2) +
  # stat_summary(geom = "point", fun = "median", col = "black",
  #              size = 7, shape=95) +
  ylab(expression(paste("IFN", gamma, " (mIU/ml)"))) +
  xlab("EI-S1-IgG") +
  scale_y_continuous(trans="log2", breaks=c(0.01, 5, 40, 250, 1000),
                     labels=c(0, 5, 40, 250, 1000)) +
  stat_n_text(size=3, y.pos=-8) +
  stat_pvalue_manual(stat.test, y.position=12)

ggsave(here_output("Figure_S1.pdf"), height=4, width=6)
ggsave(here_output("Figure_S1.png"), height=4, width=6)

# WITH BONFERRONI CORRECTION
stat.test <- wt %>%
  group_by(Antigen) %>%
  wilcox_test(value ~ Seropositivity.IgG) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
data.frame(stat.test)

stat.test <- stat.test %>% add_xy_position(x = "Seropositivity.IgG")
stat.test
data.long %>%
  mutate(value = ifelse(value == 0, 0.01, value)) %>%
  filter(Group.5 %in% "PCR+ Seronegative") %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  filter(!is.na(IgG)) %>%
  ggplot(aes(x=Seropositivity.IgG, y=value)) +
  geom_jitter(width=0.1, alpha=0.4, color=cbPalette[8]) +
  facet_wrap(~(Antigen), ncol=4) +
  geom_hline(yintercept = 40, col="grey", lty=2) +
  # stat_summary(geom = "point", fun = "median", col = "black",
  #              size = 7, shape=95) +
  ylab(expression(paste("IFN", gamma, " (mU/ml)"))) +
  xlab("EI-S1-IgG") +
  scale_y_continuous(trans="log2", breaks=c(0.01, 5, 40, 250, 1000),
                     labels=c(0, 5, 40, 250, 1000)) +
  stat_n_text(size=3, y.pos=-8) +
  stat_pvalue_manual(stat.test, y.position=12)

#ggsave("Fig-SI-4-bonferroni-adjustment.pdf", height=4, width=6)


#------------------------------------------------ SI Figure
# Figure for sampling dates (Kinetics)
data.long %>%
  filter(Antigen %in% c("NC", "M","SCT", "SNT")) %>%
  filter(!is.na(Group.5)) %>%
  filter(Group.3 %in% "PCR+") %>%
  mutate(value = ifelse(value == 0, 0.1, value)) %>%
  mutate(datediff = VisitDate - TestDate) %>%
  mutate(datediff = as.numeric(datediff)) %>%
  filter(datediff > 100) %>%
  ggplot() +
  geom_point(aes(x=datediff, y=value, color=Group.5), alpha=0.5) +
  scale_y_continuous(trans="log2", breaks=c(0.1, 5, 40, 250, 2500),
                     labels=c("0", "5", "40", "250", "2500")) +
  xlab("Days between PCR test and sample measurement") +
  ylab(expression(paste("IFN", gamma, " (mU/ml)"))) +
#  facet_wrap(~Antigen) +
  geom_hline(yintercept = 40, col="grey", lty=2) +
    stat_smooth(aes(x=datediff, y=value), lty=1, se=FALSE,
                formula = y ~ x, method="lm", color="grey50") +
  # stat_smooth(data = .%>%filter(Antigen %in% "Membrane"),
  #             aes(x=datediff, y=value), lty=1, se=FALSE,
  #             color="red",
  #             formula = y ~ x, method="lm") +
  # stat_smooth(data = .%>%filter(Antigen %in% "Capsid"),
  #             aes(x=datediff, y=value), lty=1, se=FALSE,
  #             color="blue",
  #             formula = y ~ x, method="lm") +
  # 
   # stat_smooth(data=.%>%filter(Group.5 %in% "PCR+ Seropositive"),
   #             aes(x=datediff, y=value), lwd=1, se=FALSE,
   #             color=cbPalette[4],
   #             formula = y ~ poly(x,2), method="lm") +
  # stat_smooth(data=.%>%filter(Group.5 %in% "PCR+ Seronegative"),
  #                          aes(x=datediff, y=value), lty=1, se=FALSE,
  #                          color=cbPalette[2],
  #                          formula = y ~ x, method="lm") +
  # stat_smooth(data=.%>%filter(Group.5 %in% "PCR+ Seronegative"),
  #             aes(x=datediff, y=value), lwd=1, se=FALSE,
  #             color=cbPalette[2],
  #             formula = y ~ poly(x,2), method="lm") +
   scale_color_manual(values=cbPalette[c(8,4)], name="Group") +
  scale_fill_manual(values=cbPalette[c(8,4)], name="Group") +
  facet_wrap(~fct_inorder(Antigen))
  
#ggsave("Fig-Time-PCR-Sampling.pdf", height=4, width=6)

```


