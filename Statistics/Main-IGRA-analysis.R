#---
#title: "T-cell manuscript - Data visualization"
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

######################################################## Load data
load(here_data_tcell("data-raw-tcell-210315.RData"))
head(data)

# Minor modifications to the data
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

#---------------------------------------------------- Figure 2

#---------------------------------------------------- Figure 2a
s <- data.long %>%
  filter(Group.5 %in% c("Controls", "PCR+ Seropositive")) %>%
  filter(Antigen %in% c("NC", "M", "SNT", "SCT")) %>%
  mutate(value = ifelse(value == 0, 0.01, value)) %>%
  mutate(Group.5 = case_when(
    Group.5 == "Controls" ~ "Controls",
    Group.5 == "PCR+ Seropositive" ~ "PCR+\nSeropositive")) %>%
  filter(!is.na(value))
s
ss <- s %>%
  group_by(Antigen, Group.5) %>%
  count() %>%
  rename(n.ind = n)
ss

sss <- s %>%
  mutate(zero = ifelse(value == 0.01, TRUE, FALSE)) %>%
  group_by(Antigen, Group.5, zero) %>%
  count() 
sss  

fig2a <- s %>%
  select(LabId, Group.5, Antigen, value) %>%
  ggplot(aes(x=Group.5, y=value, color=Group.5)) +
  geom_jitter(alpha=0.4, width=0.2) +
  facet_wrap(~fct_inorder(Antigen), ncol=4) +
  scale_y_continuous(trans="log2", breaks=c(0.01, 5, 40, 250, 2500),
                     labels=c(0, 5, 40, 250, 2500)) +
  xlab("") + ylab(expression(paste("IFN", gamma, " (mIU/ml)"))) +
  theme(legend.position = "none") +
    scale_color_manual(values=cbPalette[c(2,4)]) +
  geom_hline(yintercept=40, color="grey", lty=2) +
  stat_summary(geom = "point", fun = "median", col = "black",
               size = 7, shape=95) +
  geom_text(data=sss %>% filter(zero), aes(y=0.05,
            label=paste0("n=", n)), size=3) +
  geom_text(data=sss %>% filter(!zero), aes(y=20000,
            label=paste0("n=", n)), size=3) +
  stat_n_text(size=3, y.pos=-8) 
fig2a

#---------------------------------------------------- Figure 2b
n.ind <- data.long %>%
  filter(!is.na(Group.5)) %>%
  filter(Group.5 %in% c("Exposed Seropositive", "PCR+ Seropositive",
                        "PCR+ Seronegative")) %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  filter(!is.na(value)) %>%
  group_by(Antigen) %>%
  count() %>%
  rename(n.ind = n)
n.ind

s <- data.long %>%
  filter(Group.5 %in% c("PCR+ Seropositive", "PCR+ Seronegative",
                        "Exposed Seropositive")) %>%
  filter(Antigen %in% c("NC", "M", "SNT", "SCT")) %>%
  mutate(value = ifelse(value == 0, 0.01, value)) %>%
  mutate(Group.5 = case_when(
    Group.5 == "PCR+ Seropositive" ~ "PCR+\nSeropositive",
    Group.5 == "PCR+ Seronegative" ~ "PCR+\nSeronegative",
    Group.5 == "Exposed Seropositive" ~ "Exposed\nSeropositive")) %>%
  filter(!is.na(value))
s
ss <- s %>%
  group_by(Antigen, Group.5) %>%
  count()
ss

sss <- s %>%
  mutate(zero = ifelse(value == 0.01, TRUE, FALSE)) %>%
  group_by(Antigen, zero) %>%
  count() %>%
  left_join(n.ind) %>%
  mutate(Antigen = paste0(Antigen, "\n(n=", n.ind, ")"))
sss  

mann.whitney <- list(c(sss$Antigen[1], sss$Antigen[3]),
                     c(sss$Antigen[1], sss$Antigen[5]),
                     c(sss$Antigen[1], sss$Antigen[7]),
                     c(sss$Antigen[3], sss$Antigen[5]),
                     c(sss$Antigen[3], sss$Antigen[7]),
                     c(sss$Antigen[5], sss$Antigen[7]))

data_tests <- data.long %>%
  filter(!is.na(Group.5)) %>%
  filter(Group.5 %in% c("Exposed Seropositive", 
                        "PCR+ Seropositive","PCR+ Seronegative")) %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  left_join(n.ind) %>%
  mutate(Antigen.new = paste0(Antigen, "\n(n=", n.ind, ")")) %>%
  mutate(value = ifelse(value == 0, 0.1, value)) %>%
  select(LabId, Group.5, Antigen, Antigen.new, value, n.ind) %>%
 # mutate(value = log10(value))
  mutate(value = log2(value))

wilcox_test <- wilcox_test(formula = value~Antigen.new, 
                           data = data_tests, 
                           comparisons = mann.whitney, 
                           p.adjust.method = "bonferroni",
                           paired=TRUE)

fig2b <- s %>%
  left_join(n.ind) %>%
  mutate(Antigen = paste0(Antigen, "\n(n=", n.ind, ")")) %>%
  ggplot(aes(x=fct_inorder(Antigen), y=value)) +
  geom_jitter(aes(color=Group.5), alpha=0.4, width=0.2) +
  #scale_y_log10() +
   scale_y_continuous(trans="log2", breaks=c(0.01, 5, 40, 250, 2500),
                      labels=c(0, 5, 40, 250, 2500)) +
  xlab("") + ylab(expression(paste("IFN", gamma, " (mIU/ml)"))) +
  scale_color_manual(values=cbPalette[c(6,8,4)], name="Group",
                     breaks=c("Exposed\nSeropositive", "PCR+\nSeronegative",
                              "PCR+\nSeropositive"),
                     labels=c("Exposed Seropositive", "PCR+ Seronegative",
                              "PCR+ Seropositive")) +
  geom_hline(yintercept=40, color="grey", lty=2) +
  stat_summary(geom = "point", fun = "median", col = "black",
               size = 7, shape=95) +
 #  stat_compare_means(comparisons = mann.whitney, method = "wilcox.test",
  #                   size=3) 
  stat_pvalue_manual(wilcox_test[c(1,4,3),], label = "p.adj.signif", paired=TRUE,
                     y.position = 13, step.increase = 0.06, size=3)

  
fig2b

# Additional information for Figure 2b
# Exact p-values
fig2bpv <- s %>%
  left_join(n.ind) %>%
  mutate(Antigen = paste0(Antigen, "\n(n=", n.ind, ")")) %>%
  ggplot(aes(x=Antigen, y=value)) +
  geom_jitter(aes(color=Group.5), alpha=0.4, width=0.2) +
  #scale_y_log10() +
  scale_y_continuous(trans="log2", breaks=c(0.01, 5, 40, 250, 2500),
                     labels=c(0, 5, 40, 250, 2500)) +
  xlab("") + ylab(expression(paste("IFN", gamma, " (pg/ml)"))) +
  scale_color_manual(values=cbPalette[c(6,8,4)], name="Group",
                     breaks=c("Exposed\nSeropositive", "PCR+\nSeronegative",
                              "PCR+\nSeropositive"),
                     labels=c("Exposed Seropositive", "PCR+ Seronegative",
                              "PCR+ Seropositive")) +
  geom_hline(yintercept=40, color="grey", lty=2) +
  stat_summary(geom = "point", fun = "median", col = "black",
               size = 7, shape=95) +
  stat_pvalue_manual(wilcox_test, label = "p.adj", paired=TRUE,
                     y.position = 15, step.increase = 0.1, size=3)
fig2bpv

# Calculate medians of this plot
s %>%
  group_by(Antigen) %>%
  summarize(median = median(value))


#---------------------------------------------------- Figure 2c

s <- data.long %>%
  #filter(Group.5 %in% c("Controls", "PCR+ Seropositive")) %>%
  filter(Antigen %in% c("NC", "M", "SNT", "SCT")) %>%
  mutate(value = ifelse(value == 0, 0.01, value)) %>%
  mutate(Group.5 = case_when(
    Group.5 == "Controls" ~ "Controls",
    Group.5 == "Exposed Seronegative" ~ "Exposed\nSeronegative",
    Group.5 == "Exposed Seropositive" ~ "Exposed\nSeropositive",
    Group.5 == "PCR+ Seronegative" ~ "PCR+\nSeronegative",
    Group.5 == "PCR+ Seropositive" ~ "PCR+\nSeropositive")) %>%
  filter(!is.na(value))
s
ss <- s %>%
  group_by(Antigen, Group.5) %>%
  count() %>%
  rename(n.ind = n)
ss

sss <- s %>%
  mutate(zero = ifelse(value == 0.01, TRUE, FALSE)) %>%
  group_by(Antigen, Group.5, zero) %>%
  count() 
sss  

fig2c <- s %>%
  select(LabId, Group.5, Antigen, value) %>%
  filter(!is.na(Group.5)) %>%
  ggplot(aes(x=Group.5, y=value, color=Group.5)) +
  geom_jitter(alpha=0.4, width=0.2) +
  facet_wrap(~fct_inorder(Antigen), ncol=4) +
  scale_y_continuous(trans="log2", breaks=c(0.01, 5, 40, 250, 2500),
                     labels=c(0, 5, 40, 250, 2500)) +
  xlab("") + ylab(expression(paste("IFN", gamma, " (mIU/ml)"))) +
  theme(legend.position = "none") +
  scale_color_manual(values=cbPalette[c(2,3, 6, 8, 4)]) +
  geom_hline(yintercept=40, color="grey", lty=2) +
  stat_summary(geom = "point", fun = "median", col = "black",
               size = 7, shape=95) +
  geom_text(data=sss %>% filter(zero), aes(y=0.05,
                                           label=paste0("n=", n)), size=3) +
  geom_text(data=sss %>% filter(!zero), aes(y=20000,
                                            label=paste0("n=", n)), size=3) +
  stat_n_text(size=3, y.pos=-8) +
  theme(axis.text.x = element_text(angle = 90))

fig2c

plot_grid(
  plot_grid(fig2a, fig2b, labels=c("A", "B"),
          rel_widths=c(0.5, 0.5)),
  plot_grid(fig2c, labels="C"), nrow=2)

ggsave(here_output("Figure_2.pdf"), height=8, width=12)
ggsave(here_output("Figure_2.png"), height=8, width=12)


#---------------------------------------------------------- Figure 3

#---------------------------------------------------- Figure 3

# Figure 3A (percentage of the 3 antigens)
n.ind <- data.long %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "M", "SCT")) %>%
  select(LabId, Group.5, Antigen, value) %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | is.na(SCT), 
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  group_by(Group.5) %>%
  count() %>%
  rename(n.ind = n)
n.ind

fig3a <- data.long %>%
  select(LabId, Group.5, Antigen, value) %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "SCT", "M")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | is.na(SCT), 
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
  mutate(`Detected\nAntigens`= factor(detected)) %>%
  left_join(n.ind) %>%
  mutate(Group.5 = case_when(
    Group.5 == "Controls" ~ "Controls",
    Group.5 == "PCR+ Seropositive" ~ "PCR+\nSeropositive",
    Group.5 == "PCR+ Seronegative" ~ "PCR+\nSeronegative",
    Group.5 == "Exposed Seropositive" ~ "Exposed\nSeropositive",
    Group.5 == "Exposed Seronegative" ~ "Exposed\nSeronegative")) %>%
  mutate(Group.5 = paste(Group.5, "\n(n=", n.ind, ")")) %>%
  ggplot(aes(x=fct_rev(Group.5), y=prop, fill=`Detected\nAntigens`)) +
  geom_bar(position="stack", stat="identity", alpha=0.7) +
  ylab("Proportion (%)") +
  xlab("") +
  scale_fill_manual(values=cbPalette[c(2,4,3,6)])
fig3a

# Additional information from this figure
data.chisq.test <- data.long %>%
  select(LabId, Group.5, Antigen, value) %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  #gather(Antigen, value, -LabId, -Group.5) %>%
  mutate(threshold = ifelse(value >40, 1, 0)) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5) %>%
  summarize(detected = (sum(threshold))) %>%
  group_by(Group.5, detected) %>%
  count() %>%
  filter(Group.5 %in% c("Controls", "Exposed Seronegative")) %>%
  mutate(detected.bin = case_when(
    detected == 0 | detected == 1 ~ "Less than 2",
    detected == 2 | detected == 3 | detected == 4 ~ "More than 2")) %>%
  group_by(Group.5, detected.bin) %>%
  summarize(n2=sum(n)) %>%
  spread(detected.bin, n2) %>%
  ungroup()
data.chisq.test

M <- rbind(c(81, 77), c(4, 13))
chisq.test(M)

# Calculate percentage of seropositives with reaction to 3 antigens out of 3
table.fig3a <- data.long %>%
  select(LabId, Group.5, Antigen, value) %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "SCT", "M")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | is.na(SCT), 
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
  mutate(prop = n / sum(n)*100) 
table.fig3a
#write.csv(table.fig3a, "Table-Fig3a.csv")  

#---------------------------------------------------- Figure 3b

# Number of individuals for the 4 antigens 
# (subset of the previous one)
n.ind2 <- data.long %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  select(LabId, Group.5, Antigen, value) %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  group_by(Group.5) %>%
  count() %>%
  rename(n.ind = n)
n.ind2

fig3b <- data.long %>%
  select(LabId, Group.5, Antigen, value) %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5) %>%
  mutate(threshold = ifelse(value >40, "Above", "Below")) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5, threshold) %>%
  count() %>%
  spread(threshold, n) %>%
  gather(threshold, n, -c(LabId, Group.5)) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(threshold %in% "Above") %>%
  left_join(n.ind2) %>%
  mutate(Group.5 = case_when(
    Group.5 == "Controls" ~ "Controls",
    Group.5 == "PCR+ Seropositive" ~ "PCR+\nSeropositive",
    Group.5 == "PCR+ Seronegative" ~ "PCR+\nSeronegative",
    Group.5 == "Exposed Seropositive" ~ "Exposed\nSeropositive",
    Group.5 == "Exposed Seronegative" ~ "Exposed\nSeronegative")) %>%
  mutate(Group.5 = paste(Group.5, "\n(n=", n.ind, ")")) %>%
  ggplot(aes(y=n, x=reorder(LabId, n))) +
  geom_point(shape="|") +
  facet_wrap(~fct_rev(Group.5), scales="free_x", ncol=5) +
  ylab("Number of Detected Antigens") + xlab("Number of individuals") +
  theme(axis.text.x = element_blank()) 
fig3b  

# Additional information on Figure 3b
# Calculate percentage of seropositives with reaction to 4 out of 4 antigens
table.fig3b <- data.long %>%
  select(LabId, Group.5, Antigen, value) %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | 
                             is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5) %>%
  mutate(threshold = ifelse(value >40, "Above", "Below")) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5, threshold) %>%
  count() %>%
  spread(threshold, n) %>%
  gather(threshold, n, -c(LabId, Group.5)) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(threshold %in% "Above") %>%
  rename(Detected.antigens = n) %>%
  group_by(Group.5, Detected.antigens) %>%
  count() %>%
  group_by(Group.5) %>%
   mutate(prop = n / sum(n)*100) 
data.frame(table.fig3b)
#write.csv(table.fig3b, "Table-Fig3b.csv")

plot_grid(fig3a, fig3b, labels=c("A", "B"),
          rel_widths=c(0.5, 0.5))
ggsave(here_output("Figure_3.pdf"), height=4, width=12)
ggsave(here_output("Figure_3.png"), height=4, width=12)


# How many with IgG positives?
ids.pcrpos.seroneg <- data.long %>%
  select(LabId, Group.5, Antigen, value) %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | 
                             is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5) %>%
  mutate(threshold = ifelse(value >40, "Above", "Below")) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5, threshold) %>%
  count() %>%
  spread(threshold, n) %>%
  gather(threshold, n, -c(LabId, Group.5)) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(threshold %in% "Above") %>%
  left_join(n.ind2) %>%
  filter(Group.5 %in% "PCR+ Seronegative") %>%
  select(LabId)
data %>% filter(LabId %in% ids.pcrpos.seroneg$LabId) %>%
  group_by(Seropositivity.IgG) %>%
  count()


# Figure 3 with new seronegativity/positivity definitions

# Figure 3A (percentage of the 3 antigens)
data.long.b <- data.long %>%
  mutate(Seropositivity2 = case_when(
    Seropositivity.Roche == "Positive" | Seropositivity.IgG == "Positive" ~ "Positive",
    Seropositivity.Roche == "Negative" & Seropositivity.IgG == "Negative" ~ "Negative")) %>%
  mutate(Group.5b = case_when(
    Group.3 == "PCR+" & Seropositivity2 == "Positive" ~ "PCR+ Seropositive",
    Group.3 == "PCR+" & Seropositivity2 == "Negative" ~ "PCR+ Seronegative",
    Group.3 == "Exposed" & Seropositivity2 == "Positive" ~ "Exposed Seropositive",
    Group.3 == "Exposed" & Seropositivity2 == "Negative" ~ "Exposed Seronegative",
    Group.3 == "Controls" ~ "Controls"))

n.ind <- data.long.b %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "M", "SCT")) %>%
  select(LabId, Group.5b, Antigen, value) %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) |
                             is.na(SCT), 
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  group_by(Group.5b) %>%
  count() %>%
  rename(n.ind = n)
n.ind

fig3a.new <- data.long.b %>%
  select(LabId, Group.5b, Antigen, value) %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "SCT", "M")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | is.na(SCT), 
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5b) %>%
  mutate(threshold = ifelse(value >40, 1, 0)) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5b) %>%
  summarize(detected = (sum(threshold))) %>%
  group_by(Group.5b, detected) %>%
  count() %>%
  group_by(Group.5b) %>%
  mutate(prop = n / sum(n)*100) %>%
  mutate(`Detected\nAntigens`= factor(detected)) %>%
  left_join(n.ind) %>%
  mutate(Group.5b = case_when(
    Group.5b == "Controls" ~ "Controls",
    Group.5b == "PCR+ Seropositive" ~ "PCR+\nSeropositive",
    Group.5b == "PCR+ Seronegative" ~ "PCR+\nSeronegative",
    Group.5b == "Exposed Seropositive" ~ "Exposed\nSeropositive",
    Group.5b == "Exposed Seronegative" ~ "Exposed\nSeronegative")) %>%
  mutate(Group.5b = paste(Group.5b, "\n(n=", n.ind, ")")) %>%
  ggplot(aes(x=fct_rev(Group.5b), y=prop, fill=`Detected\nAntigens`)) +
  geom_bar(position="stack", stat="identity", alpha=0.7) +
  ylab("Proportion (%)") +
  xlab("") +
  scale_fill_manual(values=cbPalette[c(2,4,3,6)])
fig3a.new

data.chisq.test <- data.long.b %>%
  select(LabId, Group.5b, Antigen, value) %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  #gather(Antigen, value, -LabId, -Group.5) %>%
  mutate(threshold = ifelse(value >40, 1, 0)) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5b) %>%
  summarize(detected = (sum(threshold))) %>%
  group_by(Group.5b, detected) %>%
  count() %>%
  filter(Group.5b %in% c("Controls", "Exposed Seronegative")) %>%
  mutate(detected.bin = case_when(
    detected == 0 | detected == 1 ~ "Less than 2",
    detected == 2 | detected == 3 | detected == 4 ~ "More than 2")) %>%
  group_by(Group.5b, detected.bin) %>%
  summarize(n2=sum(n)) %>%
  spread(detected.bin, n2) %>%
  ungroup()
data.chisq.test

M <- rbind(c(81, 72), c(4, 9))
chisq.test(M)


# Calculate percentage of seropositives with reaction to 3 antigens out of 3
table.fig3a.new <- data.long.b %>%
  select(LabId, Group.5b, Antigen, value) %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "SCT", "M")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | is.na(SCT), 
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5b) %>%
  mutate(threshold = ifelse(value >40, 1, 0)) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5b) %>%
  summarize(detected = (sum(threshold))) %>%
  group_by(Group.5b, detected) %>%
  count() %>%
  group_by(Group.5b) %>%
  mutate(prop = n / sum(n)*100) 
table.fig3a.new
#write.csv(table.fig3a.new, "Table-Fig3a-b.csv")  

# Figure 3b new
# Number of individuals for the 4 antigens (subset of the previous one)
n.ind2 <- data.long.b %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  select(LabId, Group.5b, Antigen, value) %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | 
                             is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  group_by(Group.5b) %>%
  count() %>%
  rename(n.ind = n)
n.ind2

fig3b.new <- data.long.b %>%
  select(LabId, Group.5b, Antigen, value) %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | 
                             is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5b) %>%
  mutate(threshold = ifelse(value >40, "Above", "Below")) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5b, threshold) %>%
  count() %>%
  spread(threshold, n) %>%
  gather(threshold, n, -c(LabId, Group.5b)) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(threshold %in% "Above") %>%
  left_join(n.ind2) %>%
  mutate(Group.5b = case_when(
    Group.5b == "Controls" ~ "Controls",
    Group.5b == "PCR+ Seropositive" ~ "PCR+\nSeropositive",
    Group.5b == "PCR+ Seronegative" ~ "PCR+\nSeronegative",
    Group.5b == "Exposed Seropositive" ~ "Exposed\nSeropositive",
    Group.5b == "Exposed Seronegative" ~ "Exposed\nSeronegative")) %>%
  mutate(Group.5b = paste(Group.5b, "\n(n=", n.ind, ")")) %>%
  ggplot(aes(y=n, x=reorder(LabId, n))) +
  geom_point(shape="|") +
  facet_wrap(~fct_rev(Group.5b), scales="free_x", ncol=5) +
  ylab("Number of Detected Antigens") + xlab("Number of individuals") +
  theme(axis.text.x = element_blank()) 
fig3b.new

# Calculate percentage of seropositives with reaction to 4 out of 4 antigens
table.fig3b.new <- data.long.b %>%
  select(LabId, Group.5b, Antigen, value) %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | 
                             is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5b) %>%
  mutate(threshold = ifelse(value >40, "Above", "Below")) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5b, threshold) %>%
  count() %>%
  spread(threshold, n) %>%
  gather(threshold, n, -c(LabId, Group.5b)) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(threshold %in% "Above") %>%
  rename(Detected.antigens = n) %>%
  group_by(Group.5b, Detected.antigens) %>%
  count() %>%
  group_by(Group.5b) %>%
   mutate(prop = n / sum(n)*100) 
data.frame(table.fig3b.new)

#write.csv(table.fig3b.new, "Table-Fig3b-b.csv")

plot_grid(fig3a.new, fig3b.new, labels=c("A", "B"),
          rel_widths=c(0.5, 0.5))
# ggsave(here_output("Figure_3_New-seropositivity-definition.pdf"), 
#        height=4, width=12)
# ggsave(here_output("Figure_3_New-seropositivity-definition.png"), 
#        height=4, width=12)




#------------------------------------------------Figure 4

# Correlation Roche vs. t-cell antigens
fig4a <- data.long %>%
  filter(!is.na(Group.5)) %>%
    filter(Group.5 %in% c("Exposed Seropositive", "PCR+ Seropositive",
                        "PCR+ Seronegative")) %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  mutate(value = ifelse(value == 0, 0.1, value)) %>%
  ggplot() +
  geom_point(aes(y=value, x=Roche, color=Group.5), alpha=0.4) +
  facet_wrap(~fct_inorder(Antigen)) +
  scale_color_manual(values=cbPalette[c(6, 8, 4)], name="Group") +
  geom_hline(yintercept=40, color="grey5", lty=3) +
  geom_vline(xintercept = 0.422, color="grey5", lty=3) +
  xlab("Ro-N-Ig") + ylab(expression(paste("IFN", gamma, " (mIU/ml)"))) +
  stat_smooth(aes(x=Roche, y=value), lty=1, se=FALSE,
               color="grey50", formula = y ~ x, method="lm") +
  ggpubr::stat_cor(aes(y=value, x=Roche, hjust=0, vjust=0),
                       p.accuracy=0.001, r.accuracy=0.01,
                   label.x.npc=0.5, label.y.npc = 0.1,
                   show.legend=FALSE, size=3) +
    scale_y_continuous(trans="log2", breaks=c(0.1, 40, 2500), 
                      labels=c(0, 40, 2500)) +
  scale_x_continuous(trans="log2", breaks=c(0.5, 5, 100),
                     labels=c(0.5, 5, 100)) +
  theme(legend.position="none")
fig4a


# Figure 4b
# Number of individuals for the 4 antigens (subset of the previous one)
n.ind2 <- data.long %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  select(LabId, Group.5, Antigen, value) %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  group_by(Group.5) %>%
  count() %>%
  rename(n.ind = n)
n.ind2

r <- data %>% select(LabId, Roche)

fig4b <- data.long %>%
  select(LabId, Group.5, Antigen, value, Roche) %>%
  filter(!is.na(Group.5)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | 
                             is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5, -Roche) %>%
  mutate(threshold = ifelse(value >40, "Above", "Below")) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5, threshold) %>%
  count() %>%
  spread(threshold, n) %>%
  gather(threshold, n, -c(LabId, Group.5)) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(threshold %in% "Above") %>%
  left_join(n.ind2) %>%
  left_join(r) %>%
  filter(Group.5 %in% c("PCR+ Seropositive", "PCR+ Seronegative",
                        "Exposed Seropositive")) %>%
#   mutate(Group.5 = case_when(
# #    Group.5 == "Controls" ~ "Controls",
#     Group.5 == "PCR+ Seropositive" ~ "PCR+\nSeropositive",
#     Group.5 == "PCR+ Seronegative" ~ "PCR+\nSeronegative",
#     Group.5 == "Exposed Seropositive" ~ "Exposed\nSeropositive")) %>%#,
    #Group.5 == "Exposed Seronegative" ~ "Exposed\nSeronegative")) %>%
 # mutate(Group.5 = paste(Group.5, "(n=", n.ind, ")")) %>%
  ggplot(aes(y=Roche, x=as.numeric(n))) +
  geom_jitter(aes(color=Group.5), width=0.3, alpha=0.6) +
  # scale_y_continuous(trans="log2") +
   scale_y_continuous(trans="log2", breaks=c(0.1, 1, 10, 100), 
                       labels=c(0.1, 1, 10, 100)) +
  xlab("Detected Antigens") + ylab("Ro-N-Ig") +
  geom_hline(yintercept = 0.422, col="grey", lty=2) +
  scale_color_manual(values=cbPalette[c(6, 8, 4)], name="Group") +
  ggpubr::stat_cor(aes(color = Group.5, y=Roche, 
                       x=as.numeric(n), hjust=1, vjust=1),
                   method="spearman", cor.coef.name = "rho",
                   p.accuracy=0.001, r.accuracy=0.01,
                   label.x.npc="left", label.y.npc = 0.2,
                   show.legend=FALSE, size=3) +
  ggpubr::stat_cor(aes(y=Roche, x=as.numeric(n), hjust=0, vjust=0),
                   method="spearman", cor.coef.name = "rho",
                   p.accuracy=0.001, r.accuracy=0.01,
                   label.x.npc="left", label.y.npc = 0.25,
                   show.legend=FALSE, size=3)

fig4b

plot_grid(fig4a, fig4b, labels=c("A", "B"),
          rel_widths=c(0.45, 0.55))
ggsave(here_output("Figure_4.pdf"), height=5, width=12)
ggsave(here_output("Figure_4.png"), height=5, width=12)


# Additional Figures on the correlation of 
# IgG vs. t-cell antigens

figSI1 <- data.long %>%
  filter(!is.na(Group.5)) %>%
    filter(Group.5 %in% c("Exposed Seropositive", "PCR+ Seropositive",
                        "PCR+ Seronegative")) %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  mutate(value = ifelse(value == 0, 0.1, value)) %>%
  ggplot() +
  geom_point(aes(y=value, x=IgG, color=Group.5), alpha=0.4) +
  facet_wrap(~fct_inorder(Antigen)) +
  scale_color_manual(values=cbPalette[c(6, 8, 4)], name="Group") +
  geom_hline(yintercept=40, color="grey5", lty=3) +
  geom_vline(xintercept = 1.015, color="grey5", lty=3) +
  xlab("EI-S1-IgG") + ylab(expression(paste("IFN", gamma, " (mIU/ml)"))) +
  stat_smooth(aes(x=IgG, y=value), lty=1, se=FALSE,
               color="grey50", formula = y ~ x, method="lm") +
  ggpubr::stat_cor(aes(y=value, x=IgG, hjust=0, vjust=0),
                       p.accuracy=0.001, r.accuracy=0.01,
                   label.x.npc="left", label.y.npc = 0.25,
                   show.legend=FALSE, size=3) +
   scale_y_continuous(trans="log2", breaks=c(0.1, 40, 2500), 
                      labels=c(0, 40, 2500)) +
  scale_x_continuous(trans="log2", breaks=c(0.25, 1, 5),
                     labels=c(0.25, 1, 5)) 
figSI1

ggsave(here_output("Figure_S1.pdf"), height=5, width=8)
ggsave(here_output("Figure_S1.png"), height=5, width=8)


# Number of individuals for the 4 antigens 
# (subset of the previous one)
n.ind2 <- data.long.b %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "M", "SCT", "SNT")) %>%
  select(LabId, Group.5b, Antigen, value) %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | 
                             is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  group_by(Group.5b) %>%
  count() %>%
  rename(n.ind = n)
n.ind2

r <- data %>% select(LabId, Roche)

fig4b <- data.long.b %>%
  select(LabId, Group.5b, Antigen, value, Roche) %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | 
                             is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5b, -Roche) %>%
  mutate(threshold = ifelse(value >40, "Above", "Below")) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5b, threshold) %>%
  count() %>%
  spread(threshold, n) %>%
  gather(threshold, n, -c(LabId, Group.5b)) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(threshold %in% "Above") %>%
  left_join(n.ind2) %>%
  left_join(r) %>%
  filter(Group.5b %in% c("PCR+ Seropositive", "PCR+ Seronegative",
                        "Exposed Seropositive", "Exposed Seronegative")) %>%
#   mutate(Group.5 = case_when(
# #    Group.5 == "Controls" ~ "Controls",
#     Group.5 == "PCR+ Seropositive" ~ "PCR+\nSeropositive",
#     Group.5 == "PCR+ Seronegative" ~ "PCR+\nSeronegative",
#     Group.5 == "Exposed Seropositive" ~ "Exposed\nSeropositive")) %>%#,
    #Group.5 == "Exposed Seronegative" ~ "Exposed\nSeronegative")) %>%
 # mutate(Group.5 = paste(Group.5, "(n=", n.ind, ")")) %>%
  ggplot(aes(y=Roche, x=as.numeric(n))) +
  geom_jitter(aes(color=Group.5b), width=0.3, alpha=0.6) +
  # scale_y_continuous(trans="log2") +
   scale_y_continuous(trans="log2", breaks=c(0.1, 1, 10, 100), 
                       labels=c(0.1, 1, 10, 100)) +
  xlab("Detected Antigens") + ylab("Ro-N-Ig") +
  geom_hline(yintercept = 0.422, col="grey", lty=2) +
  scale_color_manual(values=cbPalette[c(1, 6, 8, 4)], name="Group") +
  ggpubr::stat_cor(aes(color = Group.5b, y=Roche, 
                       x=as.numeric(n), hjust=1, vjust=1),
                   method="spearman", cor.coef.name = "rho",
                   p.accuracy=0.001, r.accuracy=0.01,
                   label.x.npc="left", label.y.npc = 0.2,
                   show.legend=FALSE, size=3) +
  ggpubr::stat_cor(aes(y=Roche, x=as.numeric(n), hjust=0, vjust=0),
                   method="spearman", cor.coef.name = "rho",
                   p.accuracy=0.001, r.accuracy=0.01,
                   label.x.npc="left", label.y.npc = 0.25,
                   show.legend=FALSE, size=3)

fig4b

r <- data %>% select(LabId, IgG) %>%
  mutate(IgG = as.numeric(IgG))

fig4c <- data.long.b %>%
  select(LabId, Group.5b, Antigen, value, IgG) %>%
  filter(!is.na(Group.5b)) %>%
  filter(Antigen %in% c("NC", "SCT", "M", "SNT")) %>%
  droplevels() %>%
  spread(Antigen, value) %>%
  mutate(complete = ifelse(is.na(NC) | is.na(M) | 
                             is.na(SCT) | is.na(SNT),
                           "remove", "keep")) %>%
  filter(complete %in% "keep") %>%
  droplevels() %>%
  select(-complete) %>%
  gather(Antigen, value, -LabId, -Group.5b, -IgG) %>%
  mutate(threshold = ifelse(value >40, "Above", "Below")) %>%
  filter(!is.na(threshold)) %>%
  group_by(LabId, Group.5b, threshold) %>%
  count() %>%
  spread(threshold, n) %>%
  gather(threshold, n, -c(LabId, Group.5b)) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(threshold %in% "Above") %>%
  left_join(n.ind2) %>%
  left_join(r) %>%
  filter(Group.5b %in% c("PCR+ Seropositive", "PCR+ Seronegative",
                        "Exposed Seropositive", "Exposed Seronegative")) %>%
  ggplot(aes(y=IgG, x=as.numeric(n))) +
  geom_jitter(aes(color=Group.5b), width=0.3, alpha=0.6) +
  #scale_y_continuous(trans="log2") +
  scale_y_continuous(trans="log2", breaks=c(0.1, 1, 10, 100), 
                       labels=c(0.1, 1, 10, 100)) +
  xlab("Detected Antigens") + ylab("EI-S1-IgG") +
  geom_hline(yintercept = 1.015, col="grey", lty=2) +
  scale_color_manual(values=cbPalette[c(1, 6, 8, 4)], name="Group") +
  ggpubr::stat_cor(aes(color = Group.5b, y=IgG, 
                       x=as.numeric(n), hjust=1, vjust=1),
                   method="spearman", cor.coef.name = "rho",
                   p.accuracy=0.001, r.accuracy=0.01,
                   label.x.npc="left", label.y.npc = 0.2,
                   show.legend=FALSE, size=3) +
  ggpubr::stat_cor(aes(y=IgG, x=as.numeric(n), hjust=0, vjust=0),
                   method="spearman", cor.coef.name = "rho",
                   p.accuracy=0.001, r.accuracy=0.01,
                   label.x.npc="left", label.y.npc = 0.25,
                   show.legend=FALSE, size=3)

fig4c

plot_grid(fig4b, fig4c, labels=c("A", "B"),
          rel_widths=c(0.5, 0.5))
