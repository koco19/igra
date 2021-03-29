#---
#title: "T-cell manuscript - Blood donor comparison"
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
here_data_bd = function (...) here::here("Data-blood-donors", ...)
here_output = function (...) here::here("Output", ...) 

# Load data of t-cell
load(here_data_tcell("data-raw-tcell-210315.RData"))

data <- data %>%
  rename(Capsid = "Capsid sub x5",
         Membrane = "Membrane sub x5",
         SNT = "SNT sub x5",
         SCT = "SCT sub x5",
         `Pepmix (green)` = "grün sbu x5",
         `Pepmix (yellow)` = "yellow sub x5")
data.long <- data.long %>%
  mutate(Antigen = case_when(
    Antigen == "Capsid sub x5" ~ "Capsid",
    Antigen == "Membrane sub x5" ~ "Membrane",
    Antigen == "SNT sub x5" ~ "SNT",
    Antigen == "SCT sub x5" ~ "SCT",
    Antigen == "grün sbu x5" ~ "Pepmix (green)",
    Antigen == "yellow sub x5" ~ "Pepmix (yellow)"))

data.long <- data.long %>%
  mutate(IgG = as.numeric(IgG))

data <- data %>%
  mutate(IgG = as.numeric(IgG))

# Load blood donors data
bd.raw <- read.csv(here_data_bd("Final_Lab_Single_20200904.csv"))

head(bd.raw)
table(bd.raw$study_per_ID)


# Load function to extract last values (from lab paper)
source("~/co/CovSKMuc/Statistics/R/functions.R")

dim(bd.raw)
bd.raw <- extract_last_values(bd.raw)
dim(bd.raw)
head(bd.raw)
bd <- bd.raw %>%
  filter(study_per_ID %in% c("BloodDonor_Mar20", "BloodDonor_Oct19")) %>%
  select(tln_ID, Roche, Eur_IgG, study_per_ID) %>%
  rename(LabId = tln_ID) %>%
  rename(IgG = Eur_IgG,
         Group.5 = study_per_ID) %>%
  mutate(Group.5 = case_when(
    Group.5 == "BloodDonor_Mar20" ~ "Blood Donors (Mar20)",
    Group.5 == "BloodDonor_Oct19" ~ "Blood Donors (Oct19)"))
  
data <- bind_rows(data, bd)
dim(data)


#------------------------------------------------ Figure S1

# Boxplot with Wilcoxon test

# Set up the comparisons
my.comparisons <- list(c("Blood Donors (Oct19)", "Controls"),
                       c("Blood Donors (Oct19)", "Exposed Seronegative"),
                       c("Blood Donors (Oct19)", "PCR+ Seronegative"),
                       c("Controls", "Exposed Seronegative"),
                       c("Controls", "PCR+ Seronegative"),
                       c("Exposed Seronegative", "PCR+ Seronegative"))

data.wilcox <- data %>%
  filter(Group.5 %in% c("Blood Donors (Oct19)", "Controls",
                        "Exposed Seronegative", "PCR+ Seronegative"))
wilcox.results <- wilcox_test(Roche~Group.5, data=data.wilcox,
            p.adjust.method = "bonferroni")

fig.s1 <- data %>%
  filter(Group.5 %in% c("Blood Donors (Oct19)",
                        "Controls", "Exposed Seronegative",
                        "PCR+ Seronegative")) %>%
  # mutate(Group.5 = ifelse(Group.5 == "Blood Donors (Oct19)",
  #                         "Blood Donors", Group.5)) %>%
  ggplot(aes(x=Group.5, y=Roche)) +
  geom_jitter(alpha=0.4, aes(color=Group.5)) +
  geom_boxplot(width=0.2, col="black") +
  scale_y_log10() +
  scale_color_manual(values=cbPalette[c(1,2,3,8)], name="Group") +
  geom_hline(yintercept = 0.442, color="grey", lty=2) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("") + ylab("Ro-N-Ig") +
  stat_compare_means(label.y = 0.25, show.legend=FALSE) +
  # stat_pvalue_manual(wilcox.results, label = "p.adj.signif", 
  #                    step.increase = 0.06, size=3)
   stat_compare_means(comparisons = my.comparisons, method="wilcox.test",
                      step.increase=0.1, label="p.signif", hide.ns = TRUE) +
  theme(legend.position = "none")
fig.s1

ggsave(here_output("Figure_S1.pdf"), height=5, width=6)
ggsave(here_output("Figure_S1.png"), height=5, width=6)


# Figure S2

data %>%
  filter(Group.5 %in% c("PCR+ Seronegative")) %>%
  filter(!is.na(IgG)) %>%
  mutate(Seropositivity.IgG = case_when(
    Seropositivity.IgG == "Negative" ~ "EI-S1-IgG negative",
    Seropositivity.IgG == "Positive" ~ "EI-S1-IgG positive")) %>%
  ggplot(aes(x=Seropositivity.IgG, y=IgG)) +
   stat_summary(geom = "point", fun = "median", col = "black",
               size = 7, shape=95) +
  geom_jitter(width=0.1, alpha=0.4, color=cbPalette[8]) +
  theme_classic() +
  ylab("EI-S1-IgG") +
  xlab("") +
  geom_hline(yintercept = 1.015, col="grey", lty=2) +
     stat_n_text(size=3, y.pos=-1) 

ggsave(here_output("Figure_S2.pdf"), height=4, width=4)
ggsave(here_output("Figure_S2.png"), height=4, width=4)
