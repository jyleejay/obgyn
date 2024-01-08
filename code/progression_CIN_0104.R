library(tidyverse)
library(dplyr)
library(tidyr)
library(survival)
library(jskm)
library(cmprsk)
library(tidycmprsk)
library(ggsurvfit)
library(gtsummary)
library(data.table)
library(labelled)
library(glue)
library(naniar)
library(flextable)
theme_gtsummary_compact(set_theme = TRUE, font_size = NULL)

# import files
CIN1 <- read.csv("C:/Users/hi/Documents/CINProgression/data/CIN1_tidy_1218.csv", stringsAsFactors = F, fileEncoding = "euc-kr") %>% 
  mutate(CIN1_gr = 1) # 482

CIN2 <- read.csv("C:/Users/hi/Documents/CINProgression/data/CIN2_tidy_1218.csv", stringsAsFactors = F, fileEncoding = "euc-kr") %>% 
  mutate(CIN2_gr = 1) # 93

# combine f/u groups
entire_CIN <- CIN1 %>% 
  bind_rows(CIN2) %>% 
  mutate(CIN1_gr = as.factor(ifelse(is.na(CIN1_gr), 0, CIN1_gr)),
         CIN2_gr = as.factor(ifelse(is.na(CIN2_gr), 0, CIN2_gr)),
         initial_st = ifelse(CIN1_gr == 1, 1,
                            ifelse(CIN2_gr == 1, 2, NA))) %>% 
  select(-c(CIN1_gr, CIN2_gr))

# baseline characteristics table
tb.1 <- entire_CIN %>%
  select(event, initial_st, time, age_ca, hpvNum, starts_with("h_"), multi_inf, vacc_c, dia_pid, dia_anemia) %>% 
  mutate(event = factor(event, levels = c(0, 1, 2), labels = c("Persistent", "Progression", "Regression")),
         initial_st = factor(initial_st, levels = c(1, 2), labels = c("CIN1", "CIN2")),
         time = time,
         age_ca = factor(age_ca, levels = c("1", "2", "3", "4"), labels = c("20s", "30s", "40s", "over 50s")),
         hpvNum = factor(hpvNum, levels = c("lrNum", "Neg", "hrNum", "Unknown"), labels = c("Low risk", "Negative", "High risk", "Unknown")),
         h_16 = factor(h_16, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_18 = factor(h_18, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_52 = factor(h_52, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_58 = factor(h_58, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_31 = factor(h_31, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_33 = factor(h_33, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_45 = factor(h_45, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_35 = factor(h_35, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_39 = factor(h_39, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_51 = factor(h_51, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_56 = factor(h_56, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_59 = factor(h_59, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_66 = factor(h_66, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         h_68 = factor(h_68, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         multi_inf = factor(multi_inf, levels = c(1, 0, "Unknown"), labels = c("Yes", "No", "Unknown")),
         vacc_c = factor(vacc_c, levels = c(1, 2, 0), labels = c("Yes", "No", "Unknown")),
         dia_pid = factor(dia_pid, levels = c(1, 0), labels = c("Yes", "No")),
         dia_anemia = factor(dia_anemia, levels = c(1, 0), labels = c("Yes", "No")))

colnames(tb.1) <- c("event", "initial_st", "time", "age_ca", "hpvNum", "h_16", "h_18", "h_52", "h_58", "h_31", "h_33", "h_45", "h_35", "h_39", "h_51", "h_56", "h_59", "h_66", "h_68", "multi_inf", "vacc_c", "dia_pid", "dia_anemia")

var_label(tb.1) <- list(event = "CIN progression",
                        initial_st = "initial CIN state",
                        time = "Duration",
                        age_ca = "Diagnosed age",
                        hpvNum = "HPV virus type",
                        h_16 = "HPV high-risk type: 16",
                        h_18 = "HPV high-risk type: 18",
                        h_52 = "HPV high-risk type: 52",
                        h_58 = "HPV high-risk type: 58",
                        h_31 = "HPV high-risk type: 31",
                        h_33 = "HPV high-risk type: 33",
                        h_45 = "HPV high-risk type: 45",
                        h_35 = "HPV high-risk type: 35",
                        h_39 = "HPV high-risk type: 39",
                        h_51 = "HPV high-risk type: 51",
                        h_56 = "HPV high-risk type: 56",
                        h_59 = "HPV high-risk type: 59",
                        h_66 = "HPV high-risk type: 66",
                        h_68 = "HPV high-risk type: 68",
                        multi_inf = "HPV high-risk multiple infection",
                        vacc_c = "Vaccination",
                        dia_pid = "Pelvic inflammatory disease",
                        dia_anemia = "Anemia/Thromb")

fu_tb <- tb.1 %>%
  tbl_summary(by = event,
              type = list(c(initial_st, age_ca, hpvNum, starts_with("h_"), multi_inf, vacc_c, dia_anemia, dia_pid) ~ "categorical"),
              missing = "ifany",
              missing_text = "(Unknown)",
              statistic = list(all_continuous() ~ "{mean} ({sd})",
                               c(time) ~ "{median} ({min}, {max})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 0) %>%
  add_n() %>% 
  bold_labels() %>% 
  add_overall() %>% 
  add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9)) %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels () %>% 
  modify_caption("**Baseline characteristics**")


fu_tb %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "main_bs.docx")


# competing risk
CINdf.1 <- entire_CIN %>%
  replace_with_na(replace = list(hpvNum = "Unknown",
                                 h_16 = "Unknown",
                                 h_18 = "Unknown",
                                 h_52 = "Unknown",
                                 h_58 = "Unknown",
                                 h_31 = "Unknown",
                                 h_33 = "Unknown",
                                 h_45 = "Unknown",
                                 h_35 = "Unknown",
                                 h_39 = "Unknown",
                                 h_51 = "Unknown",
                                 h_56 = "Unknown",
                                 h_59 = "Unknown",
                                 h_66 = "Unknown",
                                 h_68 = "Unknown",
                                 multi_inf = "Unknown",
                                 vacc_c = 0)) %>% 
  select(event, initial_st, time, age_ca, hpvNum, starts_with("h_"), multi_inf, vacc_c, dia_pid, dia_anemia) %>% 
  mutate(event = factor(event, levels = c(0, 1, 2), labels = c("Persistent", "Progression", "Regression")),
         initial_st = factor(initial_st, levels = c(1, 2), labels = c("CIN1", "CIN2")),
         time = time,
         age_ca = factor(age_ca, levels = c("1", "2", "3", "4"), labels = c("20s", "30s", "40s", "over 50s")),
         hpvNum = factor(hpvNum, levels = c("lrNum", "Neg", "hrNum"), labels = c("Low risk or Negative", "Low risk or Negative", "High risk")),
         h_16 = factor(h_16, levels = c(0, 1), labels = c("No", "Yes")),
         h_18 = factor(h_18, levels = c(0, 1), labels = c("No", "Yes")),
         h_52 = factor(h_52, levels = c(0, 1), labels = c("No", "Yes")),
         h_58 = factor(h_58, levels = c(0, 1), labels = c("No", "Yes")),
         h_31 = factor(h_31, levels = c(0, 1), labels = c("No", "Yes")),
         h_33 = factor(h_33, levels = c(0, 1), labels = c("No", "Yes")),
         h_45 = factor(h_45, levels = c(0, 1), labels = c("No", "Yes")),
         h_35 = factor(h_35, levels = c(0, 1), labels = c("No", "Yes")),
         h_39 = factor(h_39, levels = c(0, 1), labels = c("No", "Yes")),
         h_51 = factor(h_51, levels = c(0, 1), labels = c("No", "Yes")),
         h_56 = factor(h_56, levels = c(0, 1), labels = c("No", "Yes")),
         h_59 = factor(h_59, levels = c(0, 1), labels = c("No", "Yes")),
         h_66 = factor(h_66, levels = c(0, 1), labels = c("No", "Yes")),
         h_68 = factor(h_68, levels = c(0, 1), labels = c("No", "Yes")),
         multi_inf = factor(multi_inf, levels = c(0, 1), labels = c("No", "Yes")),
         vacc_c = factor(vacc_c, levels = c(1, 2), labels = c("Yes", "No")),
         dia_pid = factor(dia_pid, levels = c(0, 1), labels = c("No", "Yes")),
         dia_anemia = factor(dia_anemia, levels = c(0, 1), labels = c("No", "Yes")))


colnames(CINdf.1) <- c("event", "initial_st", "time", "age_ca", "hpvNum", "h_16", "h_18", "h_52", "h_58", "h_31", "h_33", "h_45", "h_35", "h_39", "h_51", "h_56", "h_59", "h_66", "h_68", "multi_inf", "vacc_c", "dia_pid", "dia_anemia")


var_label(CINdf.1) <- list(event = "CIN progression",
                           initial_st = "initial CIN state",
                           time = "Duration",
                           age_ca = "Diagnosed age",
                           hpvNum = "HPV virus type",
                           h_16 = "HPV high-risk type: 16",
                           h_18 = "HPV high-risk type: 18",
                           h_52 = "HPV high-risk type: 52",
                           h_58 = "HPV high-risk type: 58",
                           h_31 = "HPV high-risk type: 31",
                           h_33 = "HPV high-risk type: 33",
                           h_45 = "HPV high-risk type: 45",
                           h_35 = "HPV high-risk type: 35",
                           h_39 = "HPV high-risk type: 39",
                           h_51 = "HPV high-risk type: 51",
                           h_56 = "HPV high-risk type: 56",
                           h_59 = "HPV high-risk type: 59",
                           h_66 = "HPV high-risk type: 66",
                           h_68 = "HPV high-risk type: 68",
                           multi_inf = "HPV high-risk multiple infection",
                           vacc_c = "Vaccination",
                           dia_pid = "Pelvic inflammatory disease",
                           dia_anemia = "Anemia/Thromb")

############ main plot (figure 1)
# separated plot by group
# CIN1 dataset
cin1_dt <- CINdf.1 %>% filter(initial_st == "CIN1")

# fit the model
cin1_cuminc <- cuminc(Surv(time, event) ~ 1, data = cin1_dt)

# build a risk table
cin1_risktable <- 
  cin1_cuminc %>% 
  tidy_cuminc(times = c(0, 365, 730, 1095, 1460)) |> 
  select(outcome, time, n.risk, cum.event) %>%
  {distinct(., time, n.risk) |> 
      mutate(outcome = "Numbers at Risk") |> 
      rename(stat = n.risk) |> 
      bind_rows(select(., outcome, time, stat = cum.event))} |> 
  ggplot(aes(x = time, y = factor(outcome), label = stat)) +
  geom_text(size = 3.5) +
  labs(y = NULL, x = NULL) +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.border = element_blank()
  )

# build a cuminc plot
# cum incidence estimate for plotting
# progression
cin1_ci_est_prog <- cin1_cuminc %>% 
  tbl_cuminc(
    times = 1461,
    outcome = "Progression",
    statistic = "{estimate}%",
    label_header = "**{time/365.25}-year cuminc**")

value_cin1_prog <- cin1_ci_est_prog$table_body$stat_1

# regression
cin1_ci_est_reg <- cin1_cuminc %>% 
  tbl_cuminc(
    times = 1461,
    outcome = "Regression",
    statistic = "{estimate}%",
    label_header = "**{time/365.25}-year cuminc**")

value_cin1_reg <- cin1_ci_est_reg$table_body$stat_1

# median time to event based on data set instead of analysed model data
# CIN1_mediantime <- cin1_dt %>%
#   select(event, time) %>% 
#   group_by(event) %>%
#   summarise(Median = median(time)) %>%
#   filter(event!="Persistent") %>% 
#   mutate(prog = ifelse(event == "Progression", Median, NA),
#          reg = ifelse(event == "Regression", Median, NA)) %>% 
#   fill(c(prog, reg), .direction = "updown") %>%
#   slice(1) %>% 
#   select(prog, reg)

# 0104 updated, median time based on median of estimate at 4yrs
# event = progression
# at 1462 day = 0.1327216 / 50% = 0.0663608 -> nearest value (time 793, 0.06779643)
# median of estimate = 0.04892693, nearest time = 656
CIN1_median_prog <- as.data.frame(cin1_cuminc$tidy) %>%
  filter(outcome == "Progression") %>%
  select(outcome, estimate, time) %>% 
  mutate(half_est = 0.0663608,
         half_nearest = abs(estimate - half_est)) %>% 
  slice(which.min(half_nearest))

# event = regression
# at 1462 day = 0.3793619 / 50% = 0.1896809 -> nearest value (time 635, 0.1924023)
# median of estimate = 0.2003539, nearest time = 659
CIN1_median_reg <- as.data.frame(cin1_cuminc$tidy) %>%
  filter(outcome == "Regression") %>%
  select(outcome, estimate, time) %>%
  mutate(half_est = 0.1896809,
         half_nearest = abs(estimate - half_est)) %>% 
  slice(which.min(half_nearest))

# base plot
cin1_plot_base <- cin1_cuminc %>% 
  ggcuminc(outcome = c("Progression", "Regression"),  size = 1)


# data set for median time intersection
cin1_plot_prog <- cin1_plot_base$data %>% 
  select(time, outcome, estimate) %>% 
  filter(outcome == "Progression")

cin1_plot_reg <- cin1_plot_base$data %>% 
  select(time, outcome, estimate) %>% 
  filter(outcome == "Regression")

CIN1_intersect_prog <- with(approx(cin1_plot_prog$time, cin1_plot_prog$estimate, xout = CIN1_median_prog$time),
                            data.frame(x1 = c(-Inf, x, x), y1 = c(y, y, -Inf)))

CIN1_intersect_reg <- with(approx(cin1_plot_reg$time, cin1_plot_reg$estimate, xout = CIN1_median_reg$time),
                           data.frame(x1 = c(-Inf, x, x), y1 = c(y, y, -Inf)))


# tidy plot
cin1_plot <- cin1_plot_base +
  scale_x_continuous(breaks = seq(0, 365*4, by = 365), limits = c(0, 365*4)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  labs(
    title = "Cumulative incidence of CIN1 group",
    x = "Time (Day)") +
  theme(legend.position = c(0.1, 0.9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_line(data = CIN1_intersect_prog, aes(x = x1, y = y1), linetype = 3) +
  geom_line(data = CIN1_intersect_reg, aes(x = x1, y = y1), linetype = 3) +
  annotate("text", x = 1460, y = 0.14, label = value_cin1_prog, fontface = 2) +
  annotate("text", x = 1460, y = 0.39, label = value_cin1_reg, fontface = 2) +
  annotate("text", x = CIN1_median_prog$time, y = CIN1_intersect_prog$y1[1] + 0.03, label = glue("{CIN1_median_prog$time} days", .trim = F)) +
  annotate("text", x = CIN1_median_reg$time, y = CIN1_intersect_reg$y1[1] + 0.04, label = glue("{CIN1_median_reg$time} days", .trim = F))

# align and combine plots
cin1_combined <- ggsurvfit_align_plots(list(cin1_plot, cin1_risktable))

cin1_mainplot <- patchwork::wrap_plots(
  cin1_combined[[1]], 
  cin1_combined[[2]], 
  ncol = 1,
  heights = c(1, 0.2))

cin1_mainplot






# save as ppt for manually editing
# library(rvg);library(officer)
# editable_graph <- dml(ggobj = cin1_mainplot)
# doc <- read_pptx()
# doc <- add_slide(doc)
# doc <- ph_with(x = doc, editable_graph,
#                location = ph_location_type(type = "body") )
# print(doc, target = "plot.pptx")


# CIN2 dataset
cin2_dt <- CINdf.1 %>% filter(initial_st == "CIN2")

# fit the model
cin2_cuminc <- cuminc(Surv(time, event) ~ 1, data = cin2_dt)

# build a risk table
cin2_risktable <- 
  cin2_cuminc %>% 
  tidy_cuminc(times = c(0, 365, 730, 1095, 1460)) |> 
  select(outcome, time, n.risk, cum.event) %>%
  {distinct(., time, n.risk) |> 
      mutate(outcome = "Numbers at Risk") |> 
      rename(stat = n.risk) |> 
      bind_rows(select(., outcome, time, stat = cum.event))} |> 
  ggplot(aes(x = time, y = factor(outcome), label = stat)) +
  geom_text(size = 3.5) +
  labs(y = NULL, x = NULL) +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.border = element_blank()
  )

# build a cuminc plot
# cum incidence estimate for plotting
# progression
cin2_ci_est_prog <- cin2_cuminc %>% 
  tbl_cuminc(
    times = 1461,
    outcome = "Progression",
    statistic = "{estimate}%",
    label_header = "**{time/365.25}-year cuminc**")

value_cin2_prog <- cin2_ci_est_prog$table_body$stat_1

# regression
cin2_ci_est_reg <- cin2_cuminc %>% 
  tbl_cuminc(
    times = 1461,
    outcome = "Regression",
    statistic = "{estimate}%",
    label_header = "**{time/365.25}-year cuminc**")

value_cin2_reg <- cin2_ci_est_reg$table_body$stat_1

# median time to event based on data set instead of analysed model data
# event = progression
# at 1462 day = 0.1529244 / 50% = 0.0764622 -> nearest value (time 943, 0.08066566)
CIN2_median_prog <- as.data.frame(cin2_cuminc$tidy) %>%
  filter(outcome == "Progression") %>%
  select(outcome, estimate, time) %>% 
  mutate(half_est = 0.0764622,
         half_nearest = abs(estimate - half_est)) %>% 
  slice(which.min(half_nearest))

# event = regression
# at 1462 day = 0.5678942 / 50% = 0.2839471 -> nearest value (time 414, 0.2841699)
CIN2_median_reg <- as.data.frame(cin2_cuminc$tidy) %>%
  filter(outcome == "Regression") %>%
  select(outcome, estimate, time) %>%
  mutate(half_est = 0.2839471,
         half_nearest = abs(estimate - half_est)) %>% 
  slice(which.min(half_nearest))

# base plot
cin2_plot_base <- cin2_cuminc %>% 
  ggcuminc(outcome = c("Progression", "Regression"),  size = 1)

# data set for median time intersection
cin2_plot_prog <- cin2_plot_base$data %>% 
  select(time, outcome, estimate) %>% 
  filter(outcome == "Progression")

cin2_plot_reg <- cin2_plot_base$data %>% 
  select(time, outcome, estimate) %>% 
  filter(outcome == "Regression")

cin2_intersect_prog <- with(approx(cin2_plot_prog$time, cin2_plot_prog$estimate, xout = CIN2_median_prog$time),
                            data.frame(x1 = c(-Inf, x, x), y1 = c(y, y, -Inf)))

cin2_intersect_reg <- with(approx(cin2_plot_reg$time, cin2_plot_reg$estimate, xout = CIN2_median_reg$time),
                           data.frame(x1 = c(-Inf, x, x), y1 = c(y, y, -Inf)))

# tidy plot
cin2_plot <- cin2_plot_base +
  scale_x_continuous(breaks = seq(0, 365*4, by = 365), limits = c(0, 365*4)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  labs(
    title = "Cumulative incidence of CIN2 group",
    x = "Time (Day)") +
  theme(legend.position = c(0.1, 0.9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_line(data = cin2_intersect_prog, aes(x = x1, y = y1), linetype = 3) +
  geom_line(data = cin2_intersect_reg, aes(x = x1, y = y1), linetype = 3) +
  annotate("text", x = 1460, y = 0.15, label = value_cin2_prog, fontface = 2) +
  annotate("text", x = 1460, y = 0.57, label = value_cin2_reg, fontface = 2) +
  annotate("text", x = CIN2_median_prog$time, y = cin2_intersect_prog$y1[1] + 0.02, label = glue("{CIN2_median_prog$time} days", .trim = F)) +
  annotate("text", x = CIN2_median_reg$time, y = cin2_intersect_reg$y1[1] + 0.02, label = glue("{CIN2_median_reg$time} days", .trim = F))

# align and combine plots
cin2_combined <- ggsurvfit_align_plots(list(cin2_plot, cin2_risktable))

cin2_mainplot <- patchwork::wrap_plots(
  cin2_combined[[1]], 
  cin2_combined[[2]], 
  ncol = 1,
  heights = c(1, 0.2))

cin2_mainplot

# plots by variable
# CIN1 group
# fit the model
cin1_age <- cuminc(Surv(time, event) ~ age_ca, data = cin1_dt)

# plotting
cin1_plot_age <- cin1_age %>% 
  ggcuminc(outcome = c("Progression", "Regression"),  size = 1) +
  scale_x_continuous(breaks = seq(0, 365*4, by = 365), limits = c(0, 365*4)) +
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(
    title = "4-years Cumulative incidence by categorized age - CIN1 group",
    x = "Time (Day)") +
  theme(legend.position = c(0.12, 0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# hpv number
# eliminate NAs to unwirte NA into the legend
cin1_dt_hpv <- cin1_dt %>%
  filter(!is.na(hpvNum))

cin1_hpv <- cuminc(Surv(time, event) ~ hpvNum, data = cin1_dt_hpv)

# plotting
cin1_plot_hpv <- cin1_hpv %>% 
  ggcuminc(outcome = c("Progression", "Regression"),  size = 1) +
  scale_x_continuous(breaks = seq(0, 365*4, by = 365), limits = c(0, 365*4)) +
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(
    title = "4-years Cumulative incidence by HPV - CIN1 group",
    x = "Time (Day)") +
  theme(legend.position = c(0.15, 0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# # pid
# # fit the model
# cin1_pid <- cuminc(Surv(time, event) ~ dia_pid, data = cin1_dt)
# 
# # plotting
# cin1_plot_pid <- cin1_pid %>% 
#   ggcuminc(outcome = c("Progression", "Regression"),  size = 1) +
#   scale_x_continuous(breaks = seq(0, 365*4, by = 365), limits = c(0, 365*4)) +
#   scale_y_continuous(limits = c(0, 0.6)) +
#   labs(
#     title = "4-years Cumulative incidence by PID - CIN1 group",
#     x = "Time (Day)") +
#   theme(legend.position = c(0.1, 0.7),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(colour = "black"))
# 
# # anemia
# # fit the model
# cin1_anemia <- cuminc(Surv(time, event) ~ dia_anemia, data = cin1_dt)
# 
# # plotting
# cin1_plot_anemia <- cin1_anemia %>% 
#   ggcuminc(outcome = c("Progression", "Regression"),  size = 1) +
#   scale_x_continuous(breaks = seq(0, 365*4, by = 365), limits = c(0, 365*4)) +
#   scale_y_continuous(limits = c(0, 0.6)) +
#   labs(
#     title = "4-years Cumulative incidence by Anemia/Thrombocytosis - CIN1 group",
#     x = "Time (Day)") +
#   theme(legend.position = c(0.1, 0.7),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(colour = "black"))

# CIN2 group
# fit the model
cin2_age <- cuminc(Surv(time, event) ~ age_ca, data = cin2_dt)

# plotting
cin2_plot_age <- cin2_age %>% 
  ggcuminc(outcome = c("Progression", "Regression"),  size = 1) +
  scale_x_continuous(breaks = seq(0, 365*4, by = 365), limits = c(0, 365*4)) +
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(
    title = "4-years Cumulative incidence by categorized age - CIN2 group",
    x = "Time (Day)") +
  theme(legend.position = c(0.12, 0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# hpv number
# fit the model
# eliminate NAs to unwirte NA into the legend
cin2_dt_hpv <- cin2_dt %>%
  filter(!is.na(hpvNum))

cin2_hpv <- cuminc(Surv(time, event) ~ hpvNum, data = cin2_dt_hpv)

# plotting
cin2_plot_hpv <- cin2_hpv %>% 
  ggcuminc(outcome = c("Progression", "Regression"),  size = 1) +
  scale_x_continuous(breaks = seq(0, 365*4, by = 365), limits = c(0, 365*4)) +
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(
    title = "4-years Cumulative incidence by HPV - CIN2 group",
    x = "Time (Day)") +
  theme(legend.position = c(0.15, 0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))



# univariate regression
uni_tab1 <-
  tbl_uvregression(
    data = CINdf.1,
    method = tidycmprsk::crr,
    y = Surv(time = time, event = event),
    include = c("initial_st", "age_ca", "hpvNum", "multi_inf", "h_16", "h_18", "h_52", "h_58", "h_31", "h_33", "h_45", "h_35", "h_39", "h_51", "h_56", "h_59", "h_66", "h_68", "vacc_c", "dia_anemia", "dia_pid"),
    exponentiate = TRUE,
    add_estimate_to_reference_rows = TRUE  # this adds 1 to the coef row
  ) %>% 
  # this will update the em-dash in the CI row to Ref.
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels() %>%
  modify_caption("**Univariate competing-risks regression analysis**  \n  _(Event: Progression / Competing event: Regression)_")

uni_tab1

# multivariate regression
# main table: hpvNum
multi_tab1 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + hpvNum, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels() %>%
  modify_caption("**Hazard Ratios (and 95% Confidence Intervals) for CIN Progression with CIN Regression as a competing event**")

# multi_inf
multi_tab2 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + multi_inf, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels()

multi_tab2$table_body <- multi_tab2$table_body %>% filter(startsWith(variable, "multi_"))

# h_16
multi_tab3 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_16, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels()

multi_tab3$table_body <- multi_tab3$table_body %>% filter(startsWith(variable, "h_"))

# h_18
multi_tab4 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_18, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels()

multi_tab4$table_body <- multi_tab4$table_body %>% filter(startsWith(variable, "h_"))


# h_52
multi_tab5 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_52, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels()

multi_tab5$table_body <- multi_tab5$table_body %>% filter(startsWith(variable, "h_"))

# h_58
multi_tab6 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_58, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels()

multi_tab6$table_body <- multi_tab6$table_body %>% filter(startsWith(variable, "h_"))

# h_33
multi_tab7 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_33, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels()

multi_tab7$table_body <- multi_tab7$table_body %>% filter(startsWith(variable, "h_"))

# h_45
multi_tab8 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_45, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels()

multi_tab8$table_body <- multi_tab8$table_body %>% filter(startsWith(variable, "h_"))

# h_35
multi_tab9 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_35, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels()

multi_tab9$table_body <- multi_tab9$table_body %>% filter(startsWith(variable, "h_"))

# h_51
multi_tab10 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_51, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels()

multi_tab10$table_body <- multi_tab10$table_body %>% filter(startsWith(variable, "h_"))

# h_56
multi_tab11 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_56, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>%
  italicize_levels()

multi_tab11$table_body <- multi_tab11$table_body %>% filter(startsWith(variable, "h_"))

# h_66
multi_tab12 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_66, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>%
  italicize_levels()

multi_tab12$table_body <- multi_tab12$table_body %>% filter(startsWith(variable, "h_"))

# h_68
multi_tab13 <-
  crr(Surv(time, event) ~ initial_st + age_ca + dia_pid + dia_anemia + h_68, data = CINdf.1) %>% 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE) %>% 
  modify_table_styling(
    columns = ci,
    rows = reference_row %in% TRUE,
    missing_symbol = "Ref.",
    text_format = "italic") %>% 
  add_n(location = "label") %>% 
  bold_labels() %>% 
  bold_p(t = 0.05) %>%
  italicize_levels()

multi_tab13$table_body <- multi_tab13$table_body %>% filter(startsWith(variable, "h_"))

# `columns=` argument input. Select from ‘outcome’, ‘variable’, ‘var_type’, ‘row_type’, ‘var_label’, ‘N’, ‘n’, ‘N.event’, ‘n.event’,
# ‘strata’, ‘label’, ‘stat_1’, ‘statistic’, ‘df’, ‘p.value’
###### fine-gray test for p-value for 31, 39, 59
fg_tab1 <- 
  cuminc(Surv(time, event) ~ h_31, data = CINdf.1) %>%
  tbl_cuminc(label = "HPV high-risk type: 31") %>%
  bold_labels() %>% 
  add_n(location = "label") %>% 
  add_p() %>%
  bold_p(t = 0.05) %>% 
  italicize_levels() %>% 
  modify_column_hide(columns = "stat_1") %>% 
  remove_row_type(type = "level") %>% 
  modify_footnote(c(all_stat_cols(), p.value) ~ NA) %>%
  modify_table_body(
    ~ .x %>% 
      dplyr::mutate(p.value = NA)) %>% 
  modify_table_styling(
    columns = p.value,
    rows = variable %in% c("h_31"),
    missing_symbol = "NA",
    text_format = "italic",
    footnote = "No observation when event occurred")

fg_tab2 <- 
  cuminc(Surv(time, event) ~ h_39, data = CINdf.1) %>%
  tbl_cuminc(label = "HPV high-risk type: 39") %>%
  bold_labels() %>% 
  add_n(location = "label") %>% 
  add_p() %>%
  bold_p(t = 0.05) %>% 
  italicize_levels() %>% 
  modify_column_hide(columns = "stat_1") %>% 
  remove_row_type(type = "level") %>% 
  modify_footnote(c(all_stat_cols(), p.value) ~ NA) %>%
  modify_table_styling(
    columns = p.value,
    rows = variable %in% c("h_39"),
    footnote = "Gray's test")

fg_tab3 <- 
  cuminc(Surv(time, event) ~ h_59, data = CINdf.1) %>%
  tbl_cuminc(label = "HPV high-risk type: 59") %>%
  bold_labels() %>% 
  add_p() %>%
  bold_p(t = 0.05) %>% 
  italicize_levels() %>% 
  modify_column_hide(columns = "stat_1") %>% 
  remove_row_type(type = "level") %>% 
  modify_footnote(c(all_stat_cols(), p.value) ~ NA) %>%
  modify_table_styling(
    columns = p.value,
    rows = variable %in% c("h_59"),
    footnote = "Gray's test")


stacked_tbl <- tbl_stack(list(multi_tab1, multi_tab2, multi_tab3, multi_tab4, fg_tab1, multi_tab5, multi_tab6, fg_tab2, multi_tab7, multi_tab8, multi_tab9, multi_tab10, multi_tab11, fg_tab3, multi_tab12, multi_tab13)) %>% 
  modify_table_styling(
    columns = c(estimate, ci),
    rows = variable %in% c("h_39", "h_59"),
    missing_symbol = "-",
    text_format = "italic")

stacked_tbl %>%  as_flex_table() %>%
  add_footer_lines("Competing risks survival analysis.") %>% 
  flextable::save_as_docx(path = "1220.docx")
