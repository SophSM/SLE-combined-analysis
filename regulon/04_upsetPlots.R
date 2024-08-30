# Upset plot of regulons and associated functions
# Sofia Salazar
# August 2024
# -----------

library(tidyverse)
library(UpSetR)
library(ggplot2)

regulon_dir <- "/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/regulon"

downTerms <- read.csv(file = glue::glue("{regulon_dir}/results/parentTerms_down.csv"),
                      header = T, row.names = "X") 
downTerms$direction <- "Downregulated"
upTerms <- read.csv(file = glue::glue("{regulon_dir}/results/parentTerms_up.csv"),
                      header = T, row.names = "X") 
upTerms$direction <- "Upregulated"

data_up <- upTerms %>% mutate(value = 1) %>%
  select(regulon, value, Small) %>%
  pivot_wider(names_from = Small, values_from = value, values_fill = 0) %>%
  as.data.frame()


sets_up <- colnames(data_up)[-1]

png(filename = glue::glue("{regulon_dir}/results/upsetTerms_up.png"), 
    height = 10, width = 15, res = 300, units = "cm")
upset(data_up, sets = sets_up, sets.bar.color = "firebrick", main.bar.color = "firebrick")
dev.off()

data_up_T <- upTerms %>% mutate(value = 1) %>%
  select(regulon, value, Small) %>%
  pivot_wider(names_from = regulon, values_from = value, values_fill = 0) %>%
  as.data.frame()

sets_up_T <- colnames(data_up_T)[-1]
png(filename = glue::glue("{regulon_dir}/results/upsetRegs_up.png"), 
    height = 30, width = 40, res = 300, units = "cm")
upset(data_up_T, sets = sets_up_T, sets.bar.color = "firebrick", main.bar.color = "firebrick",)

dev.off()
# -----

data_down <- downTerms %>% mutate(value = 1) %>%
  select(regulon, value, Small) %>%
  pivot_wider(names_from = Small, values_from = value, values_fill = 0) %>%
  as.data.frame()

sets_down <- colnames(data_down)[-1]
png(filename = glue::glue("{regulon_dir}/results/upsetTerms_down.png"), 
    height = 10, width = 15, res = 300, units = "cm")
upset(data_down, sets = sets_down, sets.bar.color = "dodgerblue", main.bar.color = "dodgerblue",)

dev.off()

data_down_T <- downTerms %>% mutate(value = 1) %>%
  select(regulon, value, Small) %>%
  pivot_wider(names_from = regulon, values_from = value, values_fill = 0) %>%
  as.data.frame()

sets_down_T <- colnames(data_down_T)[-1]
png(filename = glue::glue("{regulon_dir}/results/upsetRegs_down.png"), 
    height = 50, width = 40, res = 300, units = "cm")
upset(data_down_T, sets = sets_down_T, sets.bar.color = "dodgerblue", main.bar.color = "dodgerblue",)

dev.off()

# ---------

# Upset of functions and genes

regulons_targets <- read.csv(glue::glue("{regulon_dir}/results/DEregulon_targets.csv"),
                             header = T, row.names = "X")

regulonsUp_targets <- regulons_targets %>% filter(direction == "Upregulated") %>%
  select(-c(direction)) %>% full_join(upTerms, by = c("regulon" = "regulon"))


targets_up_data <- regulonsUp_targets %>% 
  distinct(geneName, Small) %>%
  mutate(value = 1) %>%
  select(geneName, value, Small) %>%
  pivot_wider(names_from = Small, values_from = value, values_fill = 0) %>%
  as.data.frame()

setsU <- colnames(targets_up_data)[-1]
png(filename = glue::glue("{regulon_dir}/results/upsetTargetRegs_up.png"), 
    height = 20, width = 20, res = 300, units = "cm")
upset(targets_up_data, sets = setsU, sets.bar.color = "firebrick",
      main.bar.color = "firebrick",order.by = "freq")

dev.off()
# ---

regulonsDown_targets <- regulons_targets %>% filter(direction == "Downregulated") %>%
  select(-c(direction)) %>% full_join(downTerms, by = c("regulon" = "regulon"))


targets_down_data <- regulonsDown_targets %>% 
  distinct(geneName, Small) %>%
  mutate(value = 1) %>%
  select(geneName, value, Small) %>%
  pivot_wider(names_from = Small, values_from = value, values_fill = 0) %>%
  as.data.frame()

setsD <- colnames(targets_down_data)[-1]
png(filename = glue::glue("{regulon_dir}/results/upsetTargetRegs_down.png"), 
    height = 20, width = 20, res = 300, units = "cm")
upset(targets_down_data, sets = setsD, sets.bar.color = "dodgerblue",
      main.bar.color = "dodgerblue",order.by = "freq")

dev.off()