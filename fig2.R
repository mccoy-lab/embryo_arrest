library(readxl)
library(tidyverse)
library(data.table)
library(cowplot)
library(cli)
library(lme4)
library(margins)

# load data and rename columns
dt <- read_excel(path = "~/Dropbox/papers/2022_handyside/obo_deid.xlsx", sheet = 2) %>%
  as.data.table() %>%
  setnames(., make.names(colnames(.)))

# give embryos a unique identifier
dt[, embryo_id := paste(Patient.ID, .I, sep = "_")]

# tested embryos that are not blastocysts are arrested embryos
dt[, is_arrested := is.na(Grade.at.biopsy)]

# rename aneuploidy categories and set the order of levels
dt[, "Copy number results" := factor(Aneuploidy.category)]
dt[`Copy number results` == 0, `Copy number results` := "Euploid"]
dt[`Copy number results` == 1, `Copy number results` := "Full only"]
dt[`Copy number results` == 2, `Copy number results` := "Full plus"]
dt[`Copy number results` == 3, `Copy number results` := "Intermed. only"]
dt[`Copy number results` == 4, `Copy number results` := "Intermed. plus"]
dt$`Copy number results` <- factor(dt$`Copy number results`,
                                   levels = c("Full only", "Full plus", 
                                              "Intermed. plus", "Intermed. only", "Euploid"))

dev_levels_dt <- dt[!duplicated(Developmental.stage.code), 
                    c("Grade.at.biopsy", "Stage.cell.number.at.arrest", "Developmental.stage.code")] %>%
  setorder(., Developmental.stage.code)
dev_levels <- c(dev_levels_dt[1:3]$Stage.cell.number.at.arrest, dev_levels_dt[4:13]$Grade.at.biopsy)

dt[, dev_level := unname(dev_levels[dt$Developmental.stage.code])]
dt$dev_level <- factor(dt$dev_level, unname(dev_levels))
dt$dev_level <- factor(dt$dev_level, levels = c(levels(dt$dev_level), " "))
dt$dev_level <- factor(dt$dev_level, 
                       levels = c(levels(dt$dev_level)[1:3], " ", levels(dt$dev_level)[4:13]))

# parental age distribution
nrow(dt[!duplicated(Patient.ID)])
mean(dt[!duplicated(Patient.ID)]$Age.at.egg.collection, na.rm = TRUE)
range(dt[!duplicated(Patient.ID)]$Age.at.egg.collection, na.rm = TRUE)

# add missing embryos
dt_complete <- dt[!is.na(dev_level)]

dt_cycle <- read_excel(path = "~/Dropbox/papers/2022_handyside/obo_deid.xlsx", sheet = 1, skip = 1) %>%
  as.data.table() %>%
  setnames(., make.names(colnames(.))) %>%
  .[1:165,] %>%
  setnames(., c("patient_id", "maternal_age", "cycle_id", "num_2pn", "num_blast", "num_arrested", 
                "num_tested_pgta", "num_tested_blast", "num_tested_arrested"))
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
dt_untested[1:313, is_arrested := TRUE]
# 622 - 612 = 10 are blastocysts
dt_untested[314:323, dev_level := "Untested blastocyst"]
dt_untested[314:323, is_arrested := FALSE]

dt_complete <- rbind(dt_complete, dt_untested)

dt_cycle[, maternal_age := as.numeric(maternal_age)]
summary(glm(data = dt_cycle, formula = cbind(num_blast, num_arrested) ~ maternal_age, 
            family = "quasibinomial"))

table(dt_complete$is_arrested) # counts of arrested and unarrested
mean(dt_cycle$num_blast) # mean number of blastocysts
mean(dt_cycle$num_2pn) # mean total number of 2pn embryos

# count euploid and aneuploid embryos
nrow(dt_complete[Aneuploidy.category == 0])
nrow(dt_complete[Aneuploidy.category != 0])

### aneuploidy by chromosome

setnames(dt, "Copy number results", "Copy.number.results")
dt_subset <- dt[, c(1:2, 6:9, 86, 13:58, 84, 85, 87)]

setnames(dt_subset, c("patient_id", "maternal_age", "cell_num_after_1st_division", "type_of_1st_division", 
                      "cell_num_after_2nd_division", "type_of_2nd_division", "Copy.number.results", 
                      paste(rep(c("copy_number", "copy_number_level"), 23), rep(c(1:22, "XY"), each = 2), sep = "_"),
                      "embryo_id", "is_arrested", "dev_level"))

dt_subset <- dt_subset %>% 
  mutate(across(starts_with("copy_number_level_"),
                ~ suppressWarnings(as.numeric(as.character(.))))) %>%
  as.data.table()

# convert to long format

dt_melted_chrom <- melt(dt_subset, measure.vars = paste(rep(c("copy_number_"), 23), c(1:22, "XY"), sep = ""))
setnames(dt_melted_chrom, "variable", "chrom")
setnames(dt_melted_chrom, "value", "ploidy")
dt_melted_chrom[, chrom := gsub("copy_number_", "chr", chrom)]

dt_melted_level <- melt(dt_subset, measure.vars = paste(rep(c("copy_number_level"), 23), c(1:22, "XY"), sep = "_"))
dt_melted_level[, chrom := gsub("copy_number_level_", "chr", variable)]
setnames(dt_melted_level, "value", "level")

dt_melted <- merge(dt_melted_chrom[, c("patient_id", "embryo_id", "maternal_age", "chrom", "ploidy", "is_arrested")], 
                   dt_melted_level[, c("embryo_id", "chrom", "level")], 
                   by = c("embryo_id", "chrom"))

dt_melted[, chrom := gsub("chr", "", chrom)]
dt_melted$chrom <- factor(dt_melted$chrom, levels = c(1:22, "XY"))

dt_melted[, is_arrested_char := as.character(NA)]
dt_melted[is_arrested == TRUE, is_arrested_char := "Arrested embryos"]
dt_melted[is_arrested == FALSE, is_arrested_char := "Blastocysts"]

# classify putative meiotic and mitotic aneuploidies

dt_melted[, is_full_char := as.character(NA)]
dt_melted[level > 0.7, is_full_char := "Full (putative meiotic)"]
dt_melted[level <= 0.7, is_full_char := "Intermediate (putative mitotic)"]

dt_melted[ploidy == "G", ploidy := "Gain"]
dt_melted[ploidy == "L", ploidy := "Loss"]

# Figure 2A
ggplot(data = dt_melted[!is.na(ploidy) & !is.na(level) & 
                          !grepl("p", ploidy) & !grepl("q", ploidy)], 
       aes(x = chrom, y = after_stat(count), fill = ploidy)) +
  geom_bar(position = "stack") +
  facet_grid(is_full_char ~ is_arrested_char) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top") +
  xlab("Chromosome") +
  ylab("Number of embryos") +
  scale_fill_discrete(name = "Whole-chromosome abnormality:")

# count whole-chromosome aneuploidies
table(dt_melted[chrom != "XY" & ploidy %in% c("Gain", "Loss")]$is_full_char)
nrow(dt_melted[chrom == "XY" & ploidy %in% c("Gain", "Loss")])

# count segmental aneuploidies
sum(table(dt_melted[chrom != "XY" & !(ploidy %in% c("Gain", "Loss"))]$ploidy))
sum(table(dt_melted[chrom == "XY" & !(ploidy %in% c("Gain", "Loss"))]$ploidy))

# test age associations
aneuploidy_summary <- dt_melted %>%
  group_by(., embryo_id) %>%
  summarize(., is_arrested = unique(is_arrested_char),
            patient_id = unique(patient_id),
            mat_age = unique(maternal_age),
            n_wchr_meiotic = sum(ploidy %in% c("Gain", "Loss") & is_full_char == "Full (putative meiotic)"), 
            n_wchr_mitotic = sum(ploidy %in% c("Gain", "Loss") & is_full_char == "Intermediate (putative mitotic)")) %>%
  as.data.table()

aneuploidy_by_age <- aneuploidy_summary %>%
  group_by(., patient_id) %>%
  summarize(., mat_age = unique(mat_age), 
            n_meiotic_embryos = sum(n_wchr_meiotic > 0), 
            n_mitotic_embryos = sum(n_wchr_mitotic > 0), 
            tot_embryos = n()) %>%
  as.data.table()

summary(glm(data = aneuploidy_by_age, 
            formula = cbind(n_meiotic_embryos, tot_embryos - n_meiotic_embryos) ~ mat_age, 
            family = "quasibinomial"))

summary(glm(data = aneuploidy_by_age, 
            formula = cbind(n_mitotic_embryos, tot_embryos - n_mitotic_embryos) ~ mat_age, 
            family = "quasibinomial"))


# test relationship between meiotic and mitotic aneuploidy

all_embryos <- unique(dt_melted$embryo_id)
all_meiotic <- unique(dt_melted[ploidy %in% c("Gain", "Loss") & is_full_char == "Full (putative meiotic)"]$embryo_id)
all_mitotic <- unique(dt_melted[ploidy %in% c("Gain", "Loss") & is_full_char == "Intermediate (putative mitotic)"]$embryo_id)
no_whole_chrom_aneuploidy <- all_embryos[!(all_embryos %in% c(all_meiotic, all_mitotic))]
meiotic_only <- all_meiotic[!(all_meiotic %in% all_mitotic)]
mitotic_only <- all_mitotic[!(all_mitotic %in% all_meiotic)]
meiotic_and_mitotic <- all_meiotic[all_meiotic %in% all_mitotic]

fisher.test(rbind(cbind(length(no_whole_chrom_aneuploidy), length(meiotic_only)),
                  cbind(length(mitotic_only), length(meiotic_and_mitotic))))

# restrict to normally fertilized blastocysts

blast_only <- dt_complete[is_arrested == FALSE]$embryo_id

all_embryos <- unique(dt_melted[embryo_id %in% blast_only]$embryo_id)
all_meiotic <- unique(dt_melted[embryo_id %in% blast_only & ploidy %in% c("Gain", "Loss") & is_full_char == "Full (putative meiotic)"]$embryo_id)
all_mitotic <- unique(dt_melted[embryo_id %in% blast_only & ploidy %in% c("Gain", "Loss") & is_full_char == "Intermediate (putative mitotic)"]$embryo_id)
no_whole_chrom_aneuploidy <- all_embryos[!(all_embryos %in% c(all_meiotic, all_mitotic))]
meiotic_only <- all_meiotic[!(all_meiotic %in% all_mitotic)]
mitotic_only <- all_mitotic[!(all_mitotic %in% all_meiotic)]
meiotic_and_mitotic <- all_meiotic[all_meiotic %in% all_mitotic]

fisher.test(rbind(cbind(length(no_whole_chrom_aneuploidy), length(meiotic_only)),
                  cbind(length(mitotic_only), length(meiotic_and_mitotic))))

# embryo grade correlation; Figure S2

embryo_grades <- data.table(str_split_fixed(dt[is_arrested == FALSE]$dev_level, n = 2, pattern = "")) %>%
  setnames(., c("ICM Grade", "TE Grade"))

embryo_grades$`ICM Grade` <- factor(embryo_grades$`ICM Grade`, levels = rev(c("A", "B", "C", "D")))
embryo_grades$`TE Grade` <- factor(embryo_grades$`TE Grade`, levels = rev(c("A", "B", "C", "D")))

ggplot(data = embryo_grades, aes(x = `ICM Grade`, y = `TE Grade`)) +
  geom_bin2d() +
  scale_fill_viridis_c(name = "Number of\nblastocysts") +
  theme_bw() +
  theme(panel.grid = element_blank())

# test relationship between meiotic and mitotic aneuploidy and embryo arrest

# mitotic and embryo arrest
all_embryos <- unique(dt_melted$embryo_id)
all_mitotic <- unique(dt_melted[ploidy %in% c("Gain", "Loss") & is_full_char == "Intermediate (putative mitotic)"]$embryo_id)
arrested_mitotic <- unique(dt_melted[is_arrested == TRUE & (embryo_id %in% all_mitotic)]$embryo_id)
arrested_no_mitotic <- unique(dt_melted[is_arrested == TRUE & !(embryo_id %in% all_mitotic)]$embryo_id)
unarrested_mitotic <- unique(dt_melted[is_arrested == FALSE & (embryo_id %in% all_mitotic)]$embryo_id)
unarrested_no_mitotic <- unique(dt_melted[is_arrested == FALSE & !(embryo_id %in% all_mitotic)]$embryo_id)

fisher.test(rbind(cbind(length(unarrested_no_mitotic), length(unarrested_mitotic)),
                  cbind(length(arrested_no_mitotic), length(arrested_mitotic))))

# meiotic and embryo arrest
all_embryos <- unique(dt_melted$embryo_id)
all_meiotic <- unique(dt_melted[ploidy %in% c("Gain", "Loss") & is_full_char == "Full (putative meiotic)"]$embryo_id)
arrested_meiotic <- unique(dt_melted[is_arrested == TRUE & (embryo_id %in% all_meiotic)]$embryo_id)
arrested_no_meiotic <- unique(dt_melted[is_arrested == TRUE & !(embryo_id %in% all_meiotic)]$embryo_id)
unarrested_meiotic <- unique(dt_melted[is_arrested == FALSE & (embryo_id %in% all_meiotic)]$embryo_id)
unarrested_no_meiotic <- unique(dt_melted[is_arrested == FALSE & !(embryo_id %in% all_meiotic)]$embryo_id)

fisher.test(rbind(cbind(length(unarrested_no_meiotic), length(unarrested_meiotic)),
                  cbind(length(arrested_no_meiotic), length(arrested_meiotic))))

# test relationshiop between meiotic and mitotic aneuploidy and initial cell divisions

# mitotic and abnormal division

all_tl_embryos <- unique(dt_subset[!is.na(cell_num_after_1st_division)]$embryo_id)
all_2n4n <- unique(dt_subset[cell_num_after_1st_division == 2 & type_of_1st_division == "N" &
                             cell_num_after_2nd_division == 4 & type_of_2nd_division == "N"]$embryo_id)
all_abnormal <- all_tl_embryos[!(all_tl_embryos %in% all_2n4n)]
all_mitotic <- unique(dt_melted[ploidy %in% c("Gain", "Loss") & is_full_char == "Intermediate (putative mitotic)"]$embryo_id)
normal_mitotic <- all_2n4n[all_2n4n %in% all_mitotic]
normal_no_mitotic <- all_2n4n[!(all_2n4n %in% all_mitotic)]
abnormal_mitotic <- all_abnormal[all_abnormal %in% all_mitotic]
abnormal_no_mitotic <- all_abnormal[!(all_abnormal %in% all_mitotic)]

fisher.test(rbind(cbind(length(normal_no_mitotic), length(abnormal_no_mitotic)),
                  cbind(length(normal_mitotic), length(abnormal_mitotic))))

# meiotic and abnormal division

all_tl_embryos <- unique(dt_subset[!is.na(cell_num_after_1st_division)]$embryo_id)
all_2n4n <- unique(dt_subset[cell_num_after_1st_division == 2 & type_of_1st_division == "N" &
                               cell_num_after_2nd_division == 4 & type_of_2nd_division == "N"]$embryo_id)
all_abnormal <- all_tl_embryos[!(all_tl_embryos %in% all_2n4n)]
all_meiotic <- unique(dt_melted[ploidy %in% c("Gain", "Loss") & is_full_char == "Full (putative meiotic)"]$embryo_id)
normal_meiotic <- all_2n4n[all_2n4n %in% all_meiotic]
normal_no_meiotic <- all_2n4n[!(all_2n4n %in% all_meiotic)]
abnormal_meiotic <- all_abnormal[all_abnormal %in% all_meiotic]
abnormal_no_meiotic <- all_abnormal[!(all_abnormal %in% all_meiotic)]

fisher.test(rbind(cbind(length(normal_no_meiotic), length(abnormal_no_meiotic)),
                  cbind(length(normal_meiotic), length(abnormal_meiotic))))

# test relationship between cell division abnormalities and embryo arrest

all_tl_embryos <- unique(dt_subset[!is.na(cell_num_after_1st_division)]$embryo_id)
all_2n4n <- unique(dt_subset[cell_num_after_1st_division == 2 & type_of_1st_division == "N" &
                             cell_num_after_2nd_division == 4 & type_of_2nd_division == "N"]$embryo_id)
all_abnormal <- all_tl_embryos[!(all_tl_embryos %in% all_2n4n)]
all_arrested <- dt_subset[is_arrested == TRUE]$embryo_id
all_unarrested <- dt_subset[is_arrested == FALSE]$embryo_id
normal_unarrested <- all_2n4n[(all_2n4n %in% all_unarrested)]
normal_arrested <- all_2n4n[(all_2n4n %in% all_arrested)]
abnormal_unarrested <- all_abnormal[all_abnormal %in% all_unarrested]
abnormal_arrested <- all_abnormal[all_abnormal %in% all_arrested]

fisher.test(rbind(cbind(length(normal_unarrested), length(abnormal_unarrested)),
                  cbind(length(normal_arrested), length(abnormal_arrested))))

# same as above, but conditioning on euploid blastocysts

all_euploid <- unique(dt_complete[Aneuploidy.category == 0]$embryo_id)
all_tl_embryos <- unique(dt_subset[!is.na(cell_num_after_1st_division) & embryo_id %in% all_euploid]$embryo_id)
all_2n4n <- unique(dt_subset[cell_num_after_1st_division == 2 & type_of_1st_division == "N" &
                             cell_num_after_2nd_division == 4 & type_of_2nd_division == "N" & 
                             embryo_id %in% all_euploid]$embryo_id)
all_abnormal <- all_tl_embryos[!(all_tl_embryos %in% all_2n4n)]
all_arrested <- dt_subset[is_arrested == TRUE & embryo_id %in% all_euploid]$embryo_id
all_unarrested <- dt_subset[is_arrested == FALSE & embryo_id %in% all_euploid]$embryo_id
normal_unarrested <- all_2n4n[(all_2n4n %in% all_unarrested)]
normal_arrested <- all_2n4n[(all_2n4n %in% all_arrested)]
abnormal_unarrested <- all_abnormal[all_abnormal %in% all_unarrested]
abnormal_arrested <- all_abnormal[all_abnormal %in% all_arrested]

fisher.test(rbind(cbind(length(normal_unarrested), length(abnormal_unarrested)),
                  cbind(length(normal_arrested), length(abnormal_arrested))))

# test relationship between cell meiotic/mitotic aneuploidies and blastocyst morphology

embryo_grade_tabs <- group_by(dt_subset[is_arrested == FALSE], Copy.number.results == "Euploid", dev_level) %>%
  summarize(n=n()) %>%
  pivot_wider(names_from = dev_level, values_from = n, values_fill = 0)

chisq.test(embryo_grade_tabs[, 2:ncol(embryo_grade_tabs)])

dt_subset[, is_meiotic := embryo_id %in% all_meiotic]
dt_subset[, is_mitotic := embryo_id %in% all_mitotic]

embryo_grade_meiotic_tabs <- group_by(dt_subset[is_arrested == FALSE], is_meiotic, dev_level) %>%
  summarize(n=n()) %>%
  pivot_wider(names_from = dev_level, values_from = n, values_fill = 0)

chisq.test(embryo_grade_meiotic_tabs[, 2:ncol(embryo_grade_tabs)])

embryo_grade_mitotic_tabs <- group_by(dt_subset[is_arrested == FALSE], is_mitotic, dev_level) %>%
  summarize(n=n()) %>%
  pivot_wider(names_from = dev_level, values_from = n, values_fill = 0)

chisq.test(embryo_grade_mitotic_tabs[, 2:ncol(embryo_grade_tabs)])

group_by(dt_subset[is_arrested == FALSE], Copy.number.results, dev_level) %>%
  summarize(n=n()) %>%
  ggplot(data = ., aes(x = dev_level, y = n, fill = Copy.number.results)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Copy.number.results ~ .)

### test relationship between aneuploidy and blastocyst day

blast_day_aneuploidy_tabs <- group_by(dt[is_arrested == FALSE], !(Copy.number.results == "Euploid"), Day.tubed) %>%
  summarize(n=n()) %>%
  as.data.table() %>%
  setnames(., c("is_aneuploid", "day_tubed", "n")) %>%
  pivot_wider(names_from = day_tubed, values_from = n, values_fill = 0) %>%
  .[, 1:4]

chisq.test(blast_day_aneuploidy_tabs[, 2:ncol(blast_day_aneuploidy_tabs)])

dt_blast <- dt[is_arrested == FALSE]
dt_blast[, is_aneuploid := !(Copy.number.results == "Euploid")]

m1 <- lmer(data = dt_blast,
           formula = Day.tubed ~ (1 | Patient.ID) + is_aneuploid)

m1_summary <- suppressWarnings(margins_summary(m1))
m1_summary$p

