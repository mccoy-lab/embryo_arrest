library(data.table)
library(tidyverse)
library(readxl)
library(ggsankey)

dt <- read_excel(path = "~/Dropbox/papers/2022_handyside/obo_deid.xlsx", sheet = 2) %>%
  as.data.table() %>%
  setnames(., make.names(colnames(.)))

dt[, "Copy number results" := factor(Aneuploidy.category)]
dt[`Copy number results` == 0, `Copy number results` := "Euploid"]
dt[`Copy number results` == 1, `Copy number results` := "Full only"]
dt[`Copy number results` == 2, `Copy number results` := "Full plus"]
dt[`Copy number results` == 3, `Copy number results` := "Intermed. only"]
dt[`Copy number results` == 4, `Copy number results` := "Intermed. plus"]

dt$`Copy number results` <- factor(dt$`Copy number results`,
                                   levels = c("Full only", "Full plus", "Intermed. plus", "Intermed. only", "Euploid"))

dev_levels_dt <- dt[!duplicated(Developmental.stage.code), c("Grade.at.biopsy", "Stage.cell.number.at.arrest", "Developmental.stage.code")] %>%
  setorder(., Developmental.stage.code)
dev_levels <- c(dev_levels_dt[1:3]$Stage.cell.number.at.arrest, dev_levels_dt[4:13]$Grade.at.biopsy)

dt[, dev_level := unname(dev_levels[dt$Developmental.stage.code])]
dt$dev_level <- factor(dt$dev_level, unname(dev_levels))
dt$dev_level <- factor(dt$dev_level, levels = c(levels(dt$dev_level), " "))
dt$dev_level <- factor(dt$dev_level, levels = c(levels(dt$dev_level)[1:3], " ", levels(dt$dev_level)[4:13]))

###

dt_complete <- dt[!is.na(dev_level)]

dt_cycle <- read_excel(path = "~/Dropbox/papers/2022_handyside/OBO data (RMcC) 2-11-21.xlsm", sheet = 1, skip = 1) %>%
  as.data.table() %>%
  setnames(., make.names(colnames(.))) %>%
  .[1:165,] %>%
  setnames(., c("patient_id", "maternal_age", "date", "num_2pn", "num_blast", "num_arrested", "num_tested_pgta", "num_tested_blast", "num_tested_arrested"))
dt_cycle[, num_2pn := as.integer(num_2pn)]
dt_cycle[, num_blast := as.integer(num_blast)]
dt_cycle[, num_arrested := as.integer(num_arrested)]
dt_cycle[, num_tested_pgta := as.integer(num_tested_pgta)]
dt_cycle[, num_tested_blast := as.integer(num_tested_blast)]
dt_cycle[, num_tested_arrested := as.integer(num_tested_arrested)]

# one embryo missing, so actual total is 1180 (+ 52 = 1232)
dt_cycle[num_blast + num_arrested != num_2pn]

# add (1180 + 52) - 909 = 323 additional embryos that were not tested
dt_untested <- rbindlist(lapply(1:323, function(x) dt_complete[NA]))
# 558 - 297 + 52 = 313 are arrested embryos
dt_untested[1:313, dev_level := "Untested arrested"]
# 622 - 612 = 10 are blastocysts
dt_untested[314:323, dev_level := "Untested blastocyst"]

dt_complete <- rbind(dt_complete, dt_untested)

dt_complete[, all_embryos := TRUE]
dt_complete[, is_untested_arrest := dev_level == "Untested arrested"]
dt_complete[, is_early_arrest := dev_level == "Early"]
dt_complete[is_untested_arrest == TRUE, is_early_arrest := NA]
dt_complete[, is_mid_arrest := dev_level == "Mid"]
dt_complete[(is_untested_arrest == TRUE | is_early_arrest == TRUE), is_mid_arrest := NA]
dt_complete[, is_late_arrest := dev_level == "Late"]
dt_complete[(is_untested_arrest == TRUE | is_early_arrest == TRUE | is_mid_arrest == TRUE), is_late_arrest := NA]
dt_complete[(is_untested_arrest == TRUE | is_early_arrest == TRUE | is_mid_arrest == TRUE | is_late_arrest == TRUE), is_blast := as.character(NA)]
dt_complete[dev_level %in% c("AA", "AB", "BA", "BB", "BC", "CB", "CC", "DC", "CD", "DD"), is_blast := dev_level]
dt_complete[dev_level %in% c("AA", "AB", "BA", "BB", "CB"), is_blast := "High grade: AA, BA, AB, BB, CB"]
dt_complete[dev_level %in% c("BC", "CC", "DC", "CD", "DD"), is_blast := "Low grade: BC, CC, DC, CD, DD"]
dt_complete[dev_level == "Untested blastocyst", is_blast := "Untested"]

# TE grade C or D = low-grade
# TE grade A or B = high-grade

dt_long <- make_long(dt_complete[, c("all_embryos", "is_untested_arrest", "is_early_arrest", "is_mid_arrest", "is_late_arrest", "is_blast")], 
                     all_embryos, is_untested_arrest, is_early_arrest, is_mid_arrest, is_late_arrest, is_blast)

set.seed(4)
my_palette <- sample(paste0("#", c("68A2B7", "7D97B1", "928BAB", "BA7898", 
                                   "CC7384", "DD6F70", "EF6A5C", "F0744E", 
                                   "F17F41", "F28933", "DE9135", "B79746",
                                   "8F9C57", "68A268")), 13)

# Figure 1

ggplot(dt_long, aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                    fill = factor(node))) +
  geom_sankey(flow.alpha = .6, space = 80) +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# Figure 1 inset barplots

ggplot(data = dt_complete[!is.na(`Copy number results`) & !is.na(is_blast)], 
       aes(x = is_blast, y = ..count.., fill = `Copy number results`)) +
  geom_bar(position = "fill") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Proportion of embryos") +
  xlab("\n\n") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values =c("#e78ac3", "#fc8d62", "#66c2a5", "#a6d854", "#8da0cb")) +
  coord_cartesian(ylim = c(0, 1), clip = "off")

ggplot(data = dt_complete[!is.na(`Copy number results`) & dev_level %in% (c("Early", "Mid", "Late"))], 
       aes(x = dev_level, y = ..count.., fill = `Copy number results`)) +
  geom_bar(position = "fill") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Proportion of embryos") +
  xlab("\n\n") +
  scale_x_discrete(drop = TRUE) +
  scale_fill_manual(values =c("#e78ac3", "#fc8d62", "#66c2a5", "#a6d854", "#8da0cb")) +
  coord_cartesian(ylim = c(0, 1), clip = "off")
