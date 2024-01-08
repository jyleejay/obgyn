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


# CIN2
# import file
CIN2 <- read.csv("C:/Users/hi/Documents/CINProgression/data/CIN2_0108.csv", stringsAsFactors = F, fileEncoding = "euc-kr")

# duration from idx date
CIN2$date <- as.Date(CIN2$date)
CIN2$idx <- as.Date(CIN2$idx)
CIN2 <- CIN2 %>%
  group_by(ID) %>%
  mutate(duration = date - idx)

#### hpvNum  # idx 기준 전 후 3개월내의 기록
CIN2_hpvNum <- CIN2 %>%
  filter(!is.na(hpvNum)) %>%
  filter(hpvNum != 0) %>%
  arrange(ID, date) %>%
  group_by(ID) %>%
  select(ID, date, duration, idx, hpvNum, hr1:multi_inf) %>%   # multiple: as a HR group, multi_inf: as a individual num
  filter(duration < 91 & duration > -91) %>% 
  slice(which.min(abs(duration))) %>% 
  mutate(single = as.numeric(rowSums(across(h_16:h_68), na.rm = T) == 1),
         multi_inf = ifelse(is.na(multi_inf), 0, multi_inf)) %>% # multiple infection (selected HR numbers)
  mutate_at(vars(starts_with("h_")), as.numeric) %>% 
  mutate_at(vars(starts_with("h_")), ~replace_na(., 0)) %>%
  select(ID, hpvNum, h_16:multi_inf)

# pap results # idx 기준 전 후 3개월내의 min duration 기록
CIN2_pap <- CIN2 %>% 
  filter(!is.na(results_pap)) %>% 
  arrange(ID, date) %>%
  group_by(ID) %>%
  select(ID, date, duration, idx, results_pap) %>% 
  filter(duration < 91 & duration > -91) %>%
  slice(which.min(abs(duration))) %>% 
  select(ID, results_pap)

# Preoperative evaluation
# expected state(based on biopsy or pap results within 90days before a conducted day)
# referral state
before_biopsy <- CIN2 %>%
  mutate_at(vars(negative, hsil, nd), ~replace_na(., 0)) %>% 
  mutate(leep_duration = difftime(date, leep_date, units = "days"),
         hys_duration = difftime(date, hys_date, units = "days"),
         op_duration = ifelse((is.na(leep_duration)), hys_duration, leep_duration),
         # op_date = ifelse((is.na(leep_date)), hys_date, leep_date),
         biopsy_pre = ifelse(CIN3 == 1 , "CIN3",
                             ifelse(CIN2 == 1, "CIN2",
                                    ifelse(CIN1 == 1, "CIN1",
                                           ifelse(negative == 1, "negative",
                                                  ifelse(hsil == 1, "hsil",
                                                         ifelse(nd == 1, "not detected", NA))))))) %>%
  filter(!is.na(biopsy_pre)) %>% 
  filter(op_duration < 0) %>%
  group_by(ID) %>% 
  slice(which.min(abs(op_duration))) %>% 
  select(ID, biopsy_pre)

before_pap <- CIN2 %>% 
  filter(!(is.na(results_pap))) %>% 
  mutate(leep_duration = difftime(date, leep_date, units = "days"),
         hys_duration = difftime(date, hys_date, units = "days"),
         op_duration = ifelse((is.na(leep_duration)), hys_duration, leep_duration),
         pap_pre = results_pap) %>% 
  filter(op_duration < 0 & op_duration > -91) %>%
  group_by(ID) %>% 
  slice(which.min(abs(op_duration))) %>% 
  select(ID, pap_pre)

preop_result<- before_pap %>% 
  full_join(before_biopsy, by = "ID")

# after operation CIN state -> no need to make a new variable at this time
after_hys <- CIN2 %>%
  filter(hys_date == date) %>%  
  mutate_at(vars(negative, hsil, nd), ~replace_na(., 0)) %>% 
  mutate(biopsy_hys = ifelse(CIN3 == 1 , "CIN3",
                             ifelse(CIN2 == 1, "CIN2",
                                    ifelse(CIN1 == 1, "CIN1",
                                           ifelse(negative == 1, "negative",
                                                  ifelse(hsil == 1, "hsil",
                                                         ifelse(nd == 1, "not detected", NA))))))) %>%
  select(ID, biopsy_hys)

# after treatment pap results # added 1227
after_pap <- CIN2 %>% 
  filter(!(is.na(results_pap))) %>% 
  mutate(leep_duration = difftime(date, leep_date, units = "days"),
         hys_duration = difftime(date, hys_date, units = "days"),
         op_duration = ifelse((is.na(leep_duration)), hys_duration, leep_duration),
         pap_post = results_pap) %>% 
  filter(op_duration > 0) %>%   # & op_duration < 91 -> very few left, & op_duration < 181 as well..
  group_by(ID) %>% 
  slice(which.min(abs(op_duration))) %>% 
  select(ID, pap_post)

### postoperative
postop_result <- CIN2 %>% 
  filter(date == leep_date) %>%
  select(ID, biopsy_leep) %>% 
  bind_rows(after_hys) %>% 
  mutate(biopsy_post = ifelse(is.na(biopsy_leep), biopsy_hys, biopsy_leep)) %>% 
  select(ID, biopsy_post) %>% 
  left_join(after_pap, by = "ID")

#### age group
age_ca <- CIN2 %>%
  filter(entry == 1) %>%
  mutate(age_ca = ifelse(age < 30, "1",
                         ifelse(age > 29 & age < 40, "2",
                                ifelse(age > 39 & age < 50, "3",
                                       ifelse(age > 49, "4"))))) %>% 
  select(ID, age_ca)

#### vaccination     # yes 877, no 167
vacc <- CIN2 %>%
  mutate_at(vars(vac_gar, vac_cer, othervc, novc), as.numeric) %>%
  mutate_at(vars(vac_gar, vac_cer, othervc, novc), ~replace_na(., 0)) %>% 
  mutate_at(vars(vac_gar, vac_cer, othervc, novc), as.factor) %>% 
  mutate(vacc_c = ifelse(vac_gar == 1|vac_cer == 1|othervc == 1, 1,
                         ifelse(novc == 1, 2, 0))) %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  slice(1) %>% 
  select(ID, vacc_c)

# Immediate operation
imd <- CIN2 %>% 
  filter(imd_op == 1) %>%
  arrange(ID, date) %>%
  group_by(ID) %>%
  slice(1) %>%
  select(ID, imd_op)

#### anemia/thrombocytopenia
anemia <- CIN2 %>% 
  filter(dia_anemia == 1) %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  slice(1) %>% 
  select(ID, dia_anemia)

#### pelvic inflammatory disease
pid <- CIN2 %>% 
  filter(dia_pid == 1) %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  slice(1) %>% 
  select(ID, dia_pid)

# after entry leep
leep <- CIN2 %>%
  filter(!entry == 1) %>% 
  filter(leep == 1) %>%
  arrange(ID, date) %>%
  group_by(ID) %>%
  slice(1) %>%
  select(ID, leep)

# after entry hysterectomy
hy <- CIN2 %>% 
  filter(!entry == 1) %>% 
  filter(hys == 1) %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  slice(1) %>% 
  select(ID, hys)

# remove rows before entry & 90 days
CIN2_cleaned <- CIN2 %>%
  group_by(ID) %>%
  filter(date >= idx) %>% 
  filter(entry == 1 | (entry == 0 & duration > 90))

#################################################
# event: CIN3 or Cancer
prog_reg.1 <- CIN2_cleaned %>% 
  filter(!(entry == 1)) %>% 
  filter(CIN3 == 1 | dia_cancer == 1 | CIN1 ==1 | Neg == 1) %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  slice(1)

prog.1 <- prog_reg.1 %>% 
  filter(CIN3 == 1 | dia_cancer == 1) %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  slice(1) %>% 
  mutate(event = 1) %>% 
  select(ID, date, idx, event, dia_cancer)

prog.1_id <- prog.1$ID %>% unique()

# competing event: CIN1 or Negative or Death(none)
reg.1 <- prog_reg.1 %>%
  filter(CIN1 ==1 | Neg == 1) %>%
  arrange(ID, date) %>%
  group_by(ID) %>%
  slice(1) %>%
  mutate(event = 2) %>%
  select(ID, date, idx, event)

reg.1_id <- reg.1$ID %>% unique()

##### combine event ids
proNreg <- prog.1 %>% 
  bind_rows(reg.1)

##### cohort 별 ID
CIN2_id <- CIN2$ID %>% unique() # 340
proNreg_id <- proNreg$ID        # 178
no.1_id <- setdiff(CIN2_id, proNreg_id) # 162

#### cleaned data
CIN2_tidy.1 <- CIN2_cleaned %>%
  filter(ID %in% no.1_id) %>%
  arrange(ID, date) %>%
  group_by(ID) %>%
  slice(n()) %>%
  ungroup() %>% 
  mutate(event = 0) %>%
  select(ID, date, idx, event, dia_cancer) %>%
  bind_rows(proNreg) %>%
  mutate(time = difftime(date, idx, units = "days")) %>%
  left_join(CIN2_hpvNum, by = "ID") %>% 
  left_join(age_ca, by = "ID") %>% 
  left_join(vacc, by = "ID") %>% 
  left_join(anemia, by = "ID") %>% 
  left_join(pid, by = "ID")%>%
  left_join(imd, by = "ID") %>%
  left_join(leep, by = "ID") %>% 
  left_join(hy, by = "ID") %>%
  left_join(CIN2_pap, by = "ID") %>% 
  left_join(preop_result, by = "ID") %>% 
  left_join(postop_result, by = "ID") %>% 
  mutate_at(vars(event, dia_cancer, vacc_c, dia_anemia, dia_pid, imd_op, leep, hys), as.numeric) %>%
  mutate_at(vars(event, dia_cancer, vacc_c, dia_anemia, dia_pid, imd_op, leep, hys), ~replace_na(., 0)) %>%
  mutate_at(vars(starts_with("h_"), multi_inf), as.character) %>% 
  mutate_at(vars(hpvNum, starts_with("h_"), multi_inf), ~replace_na(., "Unknown")) %>% 
  mutate_at(vars(event, dia_cancer, vacc_c, dia_anemia, dia_pid, imd_op, leep, hys, starts_with("h_"), multi_inf), as.factor)


##### separate group by imd op
# f/u group       N = 93 (updated 1206)
CIN2_tidy_fu <- CIN2_tidy.1 %>% 
  filter(imd_op == 0)

fu_group_id <- CIN2_tidy_fu$ID

# immediate operation group  N = 247 (updated 1206)
CIN2_tidy_imop <- CIN2_tidy.1 %>% 
  filter(imd_op == 1)

imop_id_2 <- CIN2_tidy_imop$ID


# **************************** #
# find out complete regression #

#### Regression(CIN1)
reg_CIN1 <- prog_reg.1 %>%
  filter(CIN1 ==1 | Neg == 1) %>%
  arrange(ID, date) %>%
  group_by(ID) %>%
  slice(1) %>%
  filter(CIN1 == 1)

reg_CIN1_id <- reg_CIN1$ID  # 43

#### Regression(Neg)
reg_Neg <- prog_reg.1 %>%
  filter(CIN1 ==1 | Neg == 1) %>%
  arrange(ID, date) %>%
  group_by(ID) %>%
  slice(1) %>%
  filter(Neg == 1)

reg_Neg_id <- reg_Neg$ID  # 115

# among f/u patients (n = 93)
# progression   # n = 9
prog_n_id <- CIN2_tidy_fu %>% 
  filter(event == 1) %>% 
  select(ID) %>% unique()

# persistent   # n = 42
pers_n_id <- CIN2_tidy_fu %>% 
  filter(event == 0) %>% 
  select(ID) %>% unique()

# regression   # n = 42
reg_n_id <- CIN2_tidy_fu %>% 
  filter(event == 2) %>% 
  select(ID) %>% unique()

# complete regression   (19 - 1 patient Neg==1&CIN1==1 so goes in CIN1)
reg_Neg_id_fu <- CIN2_tidy_fu %>% 
  filter(ID %in% reg_Neg_id) %>% 
  select(ID) %>% unique()

# partial regression   # n = 24
reg_CIN1_id_fu <- CIN2_tidy_fu %>% 
  filter(ID %in% reg_CIN1_id) %>% 
  select(ID) %>% unique()

# goes in CIN1 id
both_reg <- intersect(reg_CIN1_id_fu, reg_Neg_id_fu)

# complete regression # n = 18
reg_Neg_id_fu <- reg_Neg_id_fu %>% 
  filter(!(ID %in% both_reg$ID)) %>% 
  mutate(com_reg = 1)

# find out never Neg (partial regression)
regreg <- CIN2_cleaned %>% 
  filter(ID %in% reg_CIN1_id_fu$ID) %>% 
  filter(Neg == 1) %>% 
  mutate(com_reg = 1) %>% 
  select(ID, com_reg) %>% 
  unique()

Any_neg_id <- regreg$ID %>% unique() # 6 -> complete regression

neverNeg <- CIN2_cleaned %>% 
  filter(ID %in% reg_CIN1_id_fu$ID) %>%
  filter(!(ID %in% regreg$ID))

neverNeg_id <- neverNeg$ID %>% unique() # 18 -> partial regression

# combine complete regression n = 24
comp_reg <- reg_Neg_id_fu %>% 
  bind_rows(regreg)


# initial event level
# 0 = persistent
# 1 = progression
# 2 = regression(partial + complete)
#  0  1  2 
# 42  9 42

# after separate regression into partial/complete
# 0 = persistent
# 1 = progression
# 2 = partial regression
# 3 = complete regression
#  0  1  2  3 
# 42  9 18 24
CIN2_tidy_fu_reg <- CIN2_tidy_fu %>% 
  left_join(comp_reg, by = "ID") %>% 
  mutate(event_new = ifelse(is.na(com_reg), event, com_reg+3),
         event_new = as.factor(event_new-1)) %>% 
  select(-event) %>% 
  mutate(event = event_new) %>% 
  select(-event_new, -com_reg)

# write.csv(CIN2_tidy_fu_reg, "C:/Users/hi/Documents/CINProgression/data/CIN2_completed_fu_0108.csv", row.names = F, fileEncoding = "euc-kr")

# among immediate treated patients (n = 247)
# progression   # n = 14
prog_n_id_imop <- CIN2_tidy_imop %>% 
  filter(event == 1) %>% 
  select(ID) %>% unique()

# persistent   # n = 120
pers_n_id_imop <- CIN2_tidy_imop %>% 
  filter(event == 0) %>% 
  select(ID) %>% unique()

# regression   # n = 113
reg_n_id_imop <- CIN2_tidy_imop %>% 
  filter(event == 2) %>% 
  select(ID) %>% unique()

# complete regression
reg_Neg_id_imop <- CIN2_tidy_imop %>%
  filter(ID %in% reg_Neg_id) %>% 
  select(ID) %>% unique()

# partial regression   # n = 19
reg_CIN1_id_imop <- CIN2_tidy_imop %>% 
  filter(ID %in% reg_CIN1_id) %>% 
  select(ID) %>% unique()

# goes in CIN1 id      # n = 2
both_reg_imop <- intersect(reg_CIN1_id_imop, reg_Neg_id_imop)

# complete regression # n = 94
reg_Neg_id_imop <- reg_Neg_id_imop %>% 
  filter(!(ID %in% both_reg_imop$ID)) %>% 
  mutate(com_reg = 1)

# find out never Neg (partial regression)
regreg_imop <- CIN2_cleaned %>% 
  filter(ID %in% reg_CIN1_id_imop$ID) %>% 
  filter(Neg == 1) %>% 
  mutate(com_reg = 1) %>% 
  select(ID, com_reg) %>% 
  unique()

Any_neg_id_imop <- regreg_imop$ID %>% unique() # 6 -> complete regression

neverNeg_imop <- CIN2_cleaned %>% 
  filter(ID %in% reg_CIN1_id_imop$ID) %>%
  filter(!(ID %in% regreg_imop$ID))

neverNeg_id_imop <- neverNeg_imop$ID %>% unique() # 13 -> partial regression

# combine complete regression n = 100
comp_reg_imop <- reg_Neg_id_imop %>% 
  bind_rows(regreg_imop)

# initial event level
# 0 = persistent
# 1 = progression
# 2 = regression(partial + complete)
#   0   1   2 
# 120  14 113

# after separate regression into partial/complete
# 0 = persistent
# 1 = progression
# 2 = partial regression
# 3 = complete regression
#   0   1   2   3 
# 120  14  13 100
CIN2_tidy_imop_reg <- CIN2_tidy_imop %>% 
  left_join(comp_reg_imop, by = "ID") %>%
  mutate(event_new = ifelse(is.na(com_reg), event, com_reg+3),
         event_new = as.factor(event_new-1)) %>% 
  select(-event) %>% 
  mutate(event = event_new) %>% 
  select(-event_new, -com_reg)

# write.csv(CIN2_tidy_imop_reg, "C:/Users/hi/Documents/CINProgression/data/CIN2_completed_imop_0108.csv", row.names = F, fileEncoding = "euc-kr")

# descriptive table
CIN2_tidy_all <- CIN2_tidy_fu_reg %>% 
  bind_rows(CIN2_tidy_imop_reg)

tb_basic <- CIN2_tidy_all %>%
  select(event, time, age_ca, hpvNum, starts_with("h_"), multi_inf, vacc_c, dia_pid, dia_anemia, imd_op, leep, hys) %>% 
  mutate(event = factor(event, levels = c(0, 1, 2, 3), labels = c("Persistent", "Progression", "Partial Regression", "Complete Regression")),
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
         dia_anemia = factor(dia_anemia, levels = c(1, 0), labels = c("Yes", "No")),
         imd_op = factor(imd_op, levels = c(1, 0), labels = c("Immediate operation", "Follow-up")),
         leep = factor(leep, levels = c(1, 0), labels = c("Yes", "No")),
         hys = factor(hys, levels = c(1, 0), labels = c("Yes", "No")))

colnames(tb_basic) <- c("event", "time", "age_ca", "hpvNum", "h_16", "h_18", "h_52", "h_58", "h_31", "h_33", "h_45", "h_35", "h_39", "h_51", "h_56", "h_59", "h_66", "h_68", "multi_inf", "vacc_c", "dia_pid", "dia_anemia", "imd_op", "leep", "hys")


var_label(tb_basic) <- list(event = "CIN progression",
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
                            dia_anemia = "Anemia/Thromb",
                            imd_op = "Immediate operation",
                            leep = "Loop Electrosurgical Excision Procedure",
                            hys = "Hysterectomy")

CIN2_tb <- tb_basic %>% 
  select(event, age_ca, time, hpvNum, starts_with("h_"), multi_inf, vacc_c, dia_anemia, dia_pid, imd_op) %>% 
  tbl_summary(by = imd_op,
              type = list(c(age_ca, hpvNum, starts_with("h_"), multi_inf, vacc_c, dia_anemia, dia_pid, event) ~ "categorical"),
              missing = "ifany",
              missing_text = "(Unknown)",
              statistic = list(all_continuous() ~ "{mean} ({sd})",
                               c(time) ~ "{median} ({min}, {max})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 0) %>%
  add_n() %>% 
  bold_labels() %>% 
  add_overall() %>% 
  add_p() %>% 
  bold_p(t = 0.05) %>% 
  italicize_levels () %>% 
  modify_caption("**CIN 2 group Baseline characteristics**  \n  _(including immediate treated patients)_")

# CIN2_tb %>%
#   as_flex_table() %>% 
#   flextable::save_as_docx(path = "cin2_table1.docx")

# ************************************************************** #
# main plot
# dataset
cin2_dt <- tb_basic %>%
  filter(imd_op == "Follow-up")

table(cin2_dt$event)

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
    outcome = "Complete Regression",
    statistic = "{estimate}%",
    label_header = "**{time/365.25}-year cuminc**")

value_cin2_reg <- cin2_ci_est_reg$table_body$stat_1

# event = progression
# at 1462 day = 0.1529244 / 50% = 0.0764622 -> nearest value (time 943, 0.08066566)
CIN2_median_prog <- as.data.frame(cin2_cuminc$tidy) %>%
  filter(outcome == "Progression") %>%
  select(outcome, estimate, time) %>% 
  mutate(half_est = 0.0764622,
         half_nearest = abs(estimate - half_est)) %>% 
  slice(which.min(half_nearest))

# event = regression
# at 1462 day = 0.3210141 / 50% = 0.1605071 -> nearest value (time 446, 0.1652683)
CIN2_median_reg <- as.data.frame(cin2_cuminc$tidy) %>%
  filter(outcome == "Complete Regression") %>%
  select(outcome, estimate, time) %>%
  mutate(half_est = 0.1605071,
         half_nearest = abs(estimate - half_est)) %>% 
  slice(which.min(half_nearest))

# base plot
cin2_plot_base <- cin2_cuminc %>%
  ggcuminc(outcome = c("Progression", "Complete Regression"),  size = 1)

# data set for median time intersection
cin2_plot_prog <- cin2_plot_base$data %>% 
  select(time, outcome, estimate) %>% 
  filter(outcome == "Progression")

cin2_plot_reg <- cin2_plot_base$data %>% 
  select(time, outcome, estimate) %>% 
  filter(outcome == "Complete Regression")

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
  theme(legend.position = c(0.2, 0.9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_line(data = cin2_intersect_prog, aes(x = x1, y = y1), linetype = 3) +
  geom_line(data = cin2_intersect_reg, aes(x = x1, y = y1), linetype = 3) +
  annotate("text", x = 1460, y = 0.15, label = value_cin2_prog, fontface = 2) +
  annotate("text", x = 1460, y = 0.32, label = value_cin2_reg, fontface = 2) +
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

#### gonna linetype change? switch each other?


word_export <- read_docx()
body_add_flextable(word_export, CIN2_tb)
body_add_par(word_export, value = "")
body_add_flextable(word_export, table2)
print(word_export, 'cin2_table1.docx')


