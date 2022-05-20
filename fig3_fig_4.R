library(readxl)
library(tidyverse)
library(data.table)
library(lme4)
library(GLMMadaptive)
library(cowplot)

# load data and rename columns
dt <- read_excel(path = "~/Dropbox/papers/2022_handyside/obo_deid.xlsx", sheet = 2) %>%
  as.data.table() %>%
  setnames(., make.names(colnames(.)))

# rename aneuploidy categories and set the order of levels
dt[, "Copy number results" := factor(Aneuploidy.category)]
dt[`Copy number results` == 0, `Copy number results` := "Euploid"]
dt[`Copy number results` == 1, `Copy number results` := "Full\nonly"]
dt[`Copy number results` == 2, `Copy number results` := "Full\nplus"]
dt[`Copy number results` == 3, `Copy number results` := "Intermed.\nonly"]
dt[`Copy number results` == 4, `Copy number results` := "Intermed.\nplus"]

dt$`Copy number results` <- factor(dt$`Copy number results`,
                                   levels = c("Full\nonly", "Full\nplus", "Intermed.\nplus", "Intermed.\nonly", "Euploid"))

dev_levels_dt <- dt[!duplicated(Developmental.stage.code), c("Grade.at.biopsy", "Stage.cell.number.at.arrest", "Developmental.stage.code")] %>%
  setorder(., Developmental.stage.code)
dev_levels <- c(dev_levels_dt[1:3]$Stage.cell.number.at.arrest, dev_levels_dt[4:13]$Grade.at.biopsy)

dt[, dev_level := unname(dev_levels[dt$Developmental.stage.code])]
dt$dev_level <- factor(dt$dev_level, unname(dev_levels))
dt$dev_level <- factor(dt$dev_level, levels = c(levels(dt$dev_level), " "))
dt$dev_level <- factor(dt$dev_level, levels = c(levels(dt$dev_level)[1:3], " ", levels(dt$dev_level)[4:13]))

# Figure 4A
p1_a <- ggplot(data = dt, aes(x = dev_level, y = ..count.., fill = `Copy number results`)) +
  geom_bar(position = "fill") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Proportion of embryos") +
  xlab("\n\n") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values =c("#e78ac3", "#fc8d62", "#66c2a5", "#a6d854", "#8da0cb")) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  geom_vline(xintercept = 4)

color_func <- colorRampPalette(c("#f7fbff", "#08306b"))

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in rev(1:(length(lut)-1))) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# Figure 4B
p1_b <- ggplot(data = dt, aes(x = dev_level, y = ..count.., fill = factor(Total.no.of.aneuploidies))) +
  geom_bar(position = "fill") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ylab("Proportion of embryos") +
  xlab("\n\n") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = color_func(22), name = "Num. of \naneuploid chroms.") +
  annotate("text", x = 2, y = -0.11, label = "\n\nStage at arrest\n(for arrested embryos)") +
  annotate("text", x = 9.5, y = -0.11, label = "\n\nGrade at biopsy\n(for developing embryos)") +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  geom_vline(xintercept = 4)

color.bar(color_func(22), min = 0, max = 22, nticks = 2, title = "Number of\naneuploid\nchroms.")

plot_grid(p1_a, p1_b, nrow = 2, align = "hv", axis = "lr", labels = c("A.", "B."))

####

# instead of embryo grades, try stratifying by day at sampling blastocysts

dt[, dev_level_v2 := dev_level]
dt[dev_level_v2 %in% c("AA", "AB", "BA", "BB", "BC", "CB", "CC", "DC", "CD", "DD"), dev_level_v2 := paste0("d", Day.tubed)]

dt$dev_level_v2 <- factor(dt$dev_level_v2, levels = c("Early", "Mid", "Late", " ", "d5", "d6", "d7"))

# Figure S3

p1_a_v2 <- ggplot(data = dt[!is.na(dev_level_v2)], aes(x = dev_level_v2, y = ..count.., fill = `Copy number results`)) +
  geom_bar(position = "fill") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Proportion of embryos") +
  xlab("\n\n") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values =c("#e78ac3", "#fc8d62", "#66c2a5", "#a6d854", "#8da0cb")) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  geom_vline(xintercept = 4)

p1_b_v2 <- ggplot(data = dt[!is.na(dev_level_v2)], aes(x = dev_level_v2, y = ..count.., fill = factor(Total.no.of.aneuploidies))) +
  geom_bar(position = "fill") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ylab("Proportion of embryos") +
  xlab("\n\n") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = color_func(22), name = "Num. of \naneuploid chroms.") +
  annotate("text", x = 2, y = -0.11, label = "\n\nStage at arrest\n(for arrested embryos)") +
  annotate("text", x = 6, y = -0.11, label = "\n\nDay at biopsy\n(for developing embryos)") +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  geom_vline(xintercept = 4)

plot_grid(p1_a_v2, p1_b_v2, nrow = 2, align = "hv", axis = "lr", labels = c("A.", "B."))

####

# for embryos with different copy number profiles, 
# what is the probability that they do versus do not arrest?

dt$`Copy number results` <- relevel(dt$`Copy number results`, ref = "Euploid")
setnames(dt, "Copy number results", "Copy.number.results")
dt[, is_arrested := (Developmental.stage.code %in% 1:3)]

### SAMPLE PROPORTION OF BLASTOCYSTS EQUAL TO OBSERVED PROPORTION OF ARRESTED EMBRYOS ###
# 612 of 622 blastocysts were analyzed with PGT-A (98.4%)
# 297 of 610 arrested embryos were analyzed with PGT-A (48.7%)
# 297 / 610 = x / 622
# x = 303

set.seed(1)
dt_downsample <- rbind(dt[is_arrested == FALSE][sample(1:612, 303, replace = FALSE)],
                       dt[is_arrested == TRUE])

m1 <- mixed_model(
  is_arrested ~ Copy.number.results + Age.at.egg.collection,
  random = ~ 1 | Patient.ID,
  data = dt_downsample,
  family = binomial())

m0 <- mixed_model(
  is_arrested ~ Copy.number.results,
  random = ~ 1 | Patient.ID,
  data = dt_downsample,
  family = binomial())

anova(m1, m0)

m0 <- mixed_model(
  is_arrested ~ 0 + Copy.number.results,
  random = ~ 1 | Patient.ID,
  data = dt_downsample,
  family = binomial())

coefs <- marginal_coefs(m0, std_errors = TRUE)$coef_table %>%
  as.data.table(keep.rownames = TRUE)
coefs[, ll_95ci := Estimate + qnorm(0.025) * Std.Err]
coefs[, ul_95ci := Estimate + qnorm(0.975) * Std.Err]

coefs[, rn := c("Euploid", "Full\nonly", "Full\nplus", "Intermed.\nplus", "Intermed.\nonly")]

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

dt_downsample$Copy.number.results <- factor(dt_downsample$Copy.number.results, 
                                            levels = c("Euploid", "Full\nonly", "Full\nplus", "Intermed.\nonly", "Intermed.\nplus"))

p2_data <- ggplot(data = dt_downsample, aes(x = Copy.number.results, y = ..count.., fill = is_arrested)) +
  geom_bar(position = "fill") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#8da0cb", "#e78ac3")) +
  ylim(0, 1) +  
  xlab("Copy number result") +
  ylab("Proportion of embryos arrested")

coefs[, est_prob := logit2prob(Estimate)]
coefs[, ll_prob := logit2prob(ll_95ci)]
coefs[, ul_prob := logit2prob(ul_95ci)]

p2_model <- ggplot(data = coefs[1:5,], aes(x = rn, 
                               y = logit2prob(Estimate), 
                               ymin = logit2prob(ll_95ci), 
                               ymax = logit2prob(ul_95ci), 
                               color = rn)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_point() +
  geom_errorbar(width = 0.25) +
  ylim(0, 1) +
  scale_color_manual(name = "", values =c("#8da0cb", "#e78ac3", "#fc8d62", "#a6d854", "#66c2a5")) +
  xlab("Copy number result") +
  ylab("Probability of arrest (95% CI)")

plot_grid(p2_data, p2_model)

# counts of number of aneuploid chromosomes (of different types) 
# versus arrest / no arrest and blastocyst grades

m7 <- mixed_model(
  is_arrested ~ Total.no.of.aneuploidies,
  random = ~ 1 | Patient.ID,
  data = dt_downsample,
  family = binomial())

coefs <- marginal_coefs(m7, std_errors = TRUE)$coef_table %>%
  as.data.table(keep.rownames = TRUE)

ilink <- family(m7)$linkinv

ndata <- data.table(Total.no.of.aneuploidies = 0:22)

m7_predict <- predict(m7, ndata, se.fit = TRUE, type_pred = "link")

ndata[, response := ilink(m7_predict$pred)]
ndata[, ul :=  ilink(m7_predict$pred + (2 * m7_predict$se.fit))]
ndata[, ll :=  ilink(m7_predict$pred - (2 * m7_predict$se.fit))]

p7_data <- ggplot(data = dt_downsample, aes(x = Total.no.of.aneuploidies, y = ..count.., fill = is_arrested)) +
  geom_bar(position = "fill") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#8da0cb", "#e78ac3")) +
  xlim(-1, 23) +
  ylim(0, 1) +  
  xlab("Number of aneuploid chromosomes") +
  ylab("Proportion of embryos arrested")

p7_model <- ggplot(data = ndata, aes(x = Total.no.of.aneuploidies, 
                         y = response, 
                         ymin = ll, 
                         ymax = ul)) +
  theme_bw() +
  theme() +
  geom_point() +
  geom_errorbar(width = 0.25) +
  xlim(-1, 23) +
  ylim(0, 1) +  
  xlab("Number of aneuploid chromosomes") +
  ylab("Probability of arrest (95% CI)")

# considering all embryos (euploid and aneuploid)
dt_downsample[Cell.no.post.1st.division == "5", Cell.no.post.1st.division := ">4"]

m8 <- glmer(
  is_arrested ~ 0 + Cell.no.post.1st.division + (1 | Patient.ID),
  data = dt_downsample,
  family = binomial())

m8_coefs <- summary(m8)$coef %>%
  as.data.table(keep.rownames = TRUE)

m8_coefs[, ll_95ci := Estimate + qnorm(0.025) * `Std. Error`]
m8_coefs[, ul_95ci := Estimate + qnorm(0.975) * `Std. Error`]

m8_coefs[, rn := gsub("Cell.no.post.1st.division", "", rn)]

m8_coefs$rn <- factor(m8_coefs$rn, levels = c("1", "2", "3", "4", "5", ">4"))

m8_results <- m8_coefs[!is.na(rn)]
m8_results[, est_prob := logit2prob(Estimate)]
m8_results[, ll_prob := logit2prob(ll_95ci)]
m8_results[, ul_prob := logit2prob(ul_95ci)]
m8_results[is.nan(ul_prob), ul_prob := 1]

p8_model <- ggplot(data = m8_results, aes(x = rn, 
                                        y = est_prob, 
                                        ymin = ll_prob, 
                                        ymax = ul_prob, 
                                        color = rn)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_point() +
  geom_errorbar(width = 0.25) +
  ylim(0, 1) +
  scale_color_brewer(palette = "Set2") +
  xlab("Cell number after first division") +
  ylab("Probability of arrest (95% CI)")

dt_downsample$Cell.no.post.1st.division <- factor(dt_downsample$Cell.no.post.1st.division, levels = c("1", "2", "3", "4", ">4"))

p8_data <- ggplot(data = dt_downsample[!is.na(Cell.no.post.1st.division)], aes(x = Cell.no.post.1st.division, y = ..count.., fill = is_arrested)) +
  geom_bar(position = "fill") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#8da0cb", "#e78ac3")) +
  xlab("Cell number after first division") +
  ylab("Proportion of embryos arrested")

# Figure 3

plot_grid(p2_data, p2_model,
          p7_data, p7_model,
          p8_data, p8_model,
          nrow = 3, ncol = 2,
          labels = paste0(LETTERS[1:6], "."))


