#---
#title: "T-cell Report in Deutsch"
#author: "Mercè Garí"
#date: '2021-04-07'
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

#source("R_rainclouds.R")

######################################################## Load data
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


#------------------------------------------------------- Figure 3A in German

# Figure 3A (percentage of the 3 antigens)
n.ind <- data.long %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("Capsid", "Membrane", "SCT")) %>%
  select(LabId, Group.5, Antigen, value) %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(Capsid) | is.na(Membrane) | is.na(SCT), 
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  group_by(Group.5) %>%
  count() %>%
  rename(n.ind = n)
n.ind

fig3a <- data.long %>%
  select(LabId, Group.5, Antigen, value) %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("Capsid", "SCT", "Membrane")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(Capsid) | is.na(Membrane) | is.na(SCT), 
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5) %>%
  mutate(threshold = ifelse(value >40, 1, 0)) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5) %>%
  summarize(detected = (sum(threshold))) %>%
  group_by(Group.5, detected) %>%
  count() %>%
  group_by(Group.5) %>%
  mutate(prop = n / sum(n)*100) %>%
  mutate(`Nachgewiesene\nEiweiße`= factor(detected)) %>%
  left_join(n.ind) %>%
  mutate(Group.5 = case_when(
     Group.5 == "Controls" ~ "Kontrollen",
    Group.5 == "PCR+ Seropositive" ~ "PCR+\nmit Antikörper",
    Group.5 == "PCR+ Seronegative" ~ "PCR+\nohne Antikörper",
    Group.5 == "Exposed Seropositive" ~ "Haushaltsmitglieder\nmit Antikörper",
    Group.5 == "Exposed Seronegative" ~ "Haushaltsmitglieder\nohne Antikörper")) %>%
  mutate(Group.5 = paste(Group.5, "\n(n=", n.ind, ")")) %>%
  ggplot(aes(x=fct_rev(fct_inorder(Group.5)), y=prop, fill=`Nachgewiesene\nEiweiße`)) +
  geom_bar(position="stack", stat="identity", alpha=0.7) +
  ylab("Anteil (%)") +
  xlab("") +
  scale_fill_manual(values=cbPalette[c(2,4,3,6)])
fig3a

ggsave(here_output("Figure_3A-Deutsch.pdf"), height=5, width=8)
ggsave(here_output("Figure_3A-Deutsch.png"), height=5, width=8)


