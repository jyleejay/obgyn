library(tidyverse)
library(survival)
library(zoo) 
library(msm)
library(mstate)
library(msm.stacked)

x <- get("cav")

`%!in%` = Negate(`%in%`)

#### Load data
# import file
initial_data <- read.csv("C:/Users/hi/Documents/CINProgression/data/cohort_imdtrt_hpv.csv", stringsAsFactors = F, fileEncoding = "euc-kr")

out_id <- ("R000012191|R000019648|R000020123|R000023511|R000024433|R000027126|R000001906|R000001979|R000002178|R000016279|R000015894|R000002350|R000022599|R000012191|R000019648|R000020123|R000023511|R000024433|R000027126|R000001906|R000001979|R000002178|R000016279|R00001589")

data <- initial_data %>% 
  ungroup() %>% 
  mutate(out = as.integer(str_detect(initial_data$ID, str_c(out_id, collapse = "|")))) %>% 
  filter(out == 0) %>% 
  select(-c(out))

init_counts <- data[!duplicated(initial_data$ID), ]  # 1283
ID_cleaned <- data[!duplicated(data$ID), ]  # 1276

# duration from idx date
data$date <- as.Date(data$date)
data$idx <- as.Date(data$idx)
data <- data %>%
  group_by(ID) %>%
  mutate(duration = date - idx)

#### hpvNum  # idx 기준 전 후 3개월내의 기록
data_hpvNum <- data %>%
  filter(!is.na(hpvNum)) %>%
  filter(hpvNum != 0) %>%
  arrange(ID, date) %>%
  group_by(ID) %>%
  select(ID, date, duration, idx, hpvNum, hr1:multi_inf) %>%   # multiple: as a HR group, multi_inf: as a individual num
  filter(duration < 91 & duration > -91) %>% 
  slice(1) %>% 
  select(ID, hpvNum, hr1:multi_inf)

hpvid <- data_hpvNum$ID %>% unique() # 733

#### age group
age_ca <- data %>%
  filter(entry == 1) %>%
  mutate(age_ca = ifelse(age < 30, "1",
                         ifelse(age > 29 & age < 40, "2",
                                ifelse(age > 39 & age < 50, "3",
                                       ifelse(age > 49, "4"))))) %>% 
  select(ID, age_ca)

# vaccination
vacc <- data %>%
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
imd <- data %>% 
  filter(imd_op == 1) %>%
  arrange(ID, date) %>%
  group_by(ID) %>%
  slice(1) %>%
  select(ID, imd_op)

#### anemia/thrombocytopenia
anemia <- data %>% 
  filter(dia_anemia == 1) %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  slice(1) %>% 
  select(ID, dia_anemia)

#### pelvic inflammatory disease
pid <- data %>% 
  filter(dia_pid == 1) %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  slice(1) %>% 
  select(ID, dia_pid)

#### cervical cancer
cancer <- data %>% 
  filter(dia_cancer == 1) %>% 
  filter(date == 최초진단일자) %>% 
  select(ID, date, dia_cancer)

# remove rows before entry
tidy_data <- data %>%
  group_by(ID) %>%
  filter(date >= idx) %>% 
  filter(entry == 1 | (entry == 0 & duration > 89)) %>% 
  select(ID, date, idx, CIN1:CIN3, negative) %>% 
  mutate(time = difftime(date, idx, units = "days")) %>%
  left_join(data_hpvNum, by = "ID") %>% 
  left_join(age_ca, by = "ID") %>% 
  left_join(vacc, by = "ID") %>% 
  left_join(anemia, by = "ID") %>% 
  left_join(pid, by = "ID")%>%
  left_join(imd, by = "ID") %>%
  left_join(cancer, by = c("ID", "date")) %>%
  mutate_at(vars(negative, dia_cancer, dia_anemia, dia_pid, imd_op, dia_anemia, dia_pid), as.numeric) %>% 
  mutate_at(vars(negative, dia_cancer, dia_anemia, dia_pid, imd_op, starts_with("h_"), multi_inf, dia_anemia, dia_pid), ~replace_na(., 0)) %>% 
  mutate_at(vars(negative, dia_cancer, dia_anemia, dia_pid, imd_op, starts_with("h_"), multi_inf, dia_anemia, dia_pid), as.factor)

# f/u group for msm analysis data set
tidy_data_fu <- tidy_data %>% 
  filter(imd_op == 0)

tidy_data_imop <- tidy_data %>% 
  filter(imd_op == 1) %>% 
  mutate(G1 = case_when(CIN1 == 1 & CIN2 == 0 & CIN3 == 0 ~ 1, TRUE ~ 0),
       G2 = case_when((CIN1 == 0 & CIN2 == 1 & CIN3 == 0)|(CIN1 == 1 & CIN2 == 1 & CIN3 == 0) ~ 1, TRUE ~ 0),
       G3 = case_when((CIN1 == 0 & CIN2 == 0 & CIN3 == 1)|(CIN1 == 0 & CIN2 == 1 & CIN3 == 1)|(CIN1 == 1 & CIN2 == 0 & CIN3 == 1)|(CIN1 == 1 & CIN2 == 1 & CIN3 == 1) ~ 1, TRUE ~ 0))
       
# event data set based on biopsy only
tidy_data_event <- tidy_data_fu %>%
  filter(CIN1 == 1|CIN2 == 1|CIN3 == 1|negative == 1|dia_cancer == 1) %>% 
  select(ID, date, idx, time, CIN1:CIN3, negative, dia_cancer, hpvNum, age_ca, vacc_c, dia_anemia, dia_pid)

# remove singletons      n = 406
singleton <- tidy_data %>%
  group_by(ID) %>%
  mutate(numofvs = row_number()) %>%
  filter(which.max(numofvs) == 1)

# finalized data        n = 411
# statemax is a variable that maximum observed state so far for each patient -> needed?
model_data <- tidy_data_event %>% 
  filter(ID %!in% singleton$ID) %>% 
  mutate(state = case_when(negative == 1 ~ 1,
                           CIN1 == 1 ~ 2,
                           CIN2 == 1 ~ 3,
                           CIN3 == 1 ~ 4,
                           dia_cancer == 1 ~ 5, TRUE ~ NA))


# Cervical cancer as a final state: absorbing event  # 409
final_data <- model_data %>% 
  group_by(ID) %>% 
  mutate(numofvs = row_number()) %>% 
  filter(state == 5) %>% 
  filter(numofvs == min(numofvs)) %>% 
  mutate(endtime = time) %>% 
  select(ID, endtime) %>% 
  right_join(model_data, by = "ID") %>% 
  mutate(endtime = ifelse(is.na(endtime), 9999, endtime)) %>% 
  filter(endtime >= time) %>% 
  select(-endtime)


# CIN1 as a first observed state in our study, not negative state
a <- final_data %>% 
  group_by(ID) %>% 
  filter(state == 1) %>% 
  filter(time == 0) %>% 
  right_join(final_data, by = "ID") %>% 
  mutate(negtime = ifelse(is.na(negtime), 9999, negtime)) %>% 
  filter(negtime >= time)

# < State observation >
# 1: negative
# 2: CIN1
# 3: CIN2
# 4: CIN3
# 5: cervical cancer

statetable.msm(state, ID, final_data)
#       to
# from   1   2   3   4   5
#    1   2   3   1   1   0
#    2   8 107  25   4   1
#    3   1  26  47  14   2
#    4   0   0   7  10   1

Q <- rbind(c(0, 0.25, 0, 0, 0), 
           c(0.25, 0, 0.25, 0, 0), 
           c(0, 0.25, 0, 0.25, 0), 
           c(0, 0, 0.25, 0, 0.25), 
           c(0, 0, 0, 0, 0))

Q.crude <- crudeinits.msm(state ~ time, ID, data = final_data, qmatrix = Q)

# Initial optimization if necessary
Q.crude[2, 1] = 0.002
# Q.crude[4, 3] = 0.002
Q.crude[4, 5] = 0.0015

# Initial modeling
init_msm <- msm(state ~ time, subject = ID, data = final_data, qmatrix = Q.crude, deathexact = 5)

# Maximum likelihood estimates
# 
# Transition intensities
# Baseline                          
# State 1 - State 1 -0.0085474 (-8.953e-02,-0.0008160)
# State 1 - State 2  0.0085474 ( 8.160e-04, 0.0895317)
# State 2 - State 1  0.0007314 ( 8.343e-05, 0.0064111)
# State 2 - State 2 -0.0019643 (-4.613e-03,-0.0008363)
# State 2 - State 3  0.0012329 ( 7.829e-04, 0.0019417)
# State 3 - State 2  0.0025448 ( 1.592e-03, 0.0040667)
# State 3 - State 3 -0.0053932 (-9.526e-03,-0.0030535)
# State 3 - State 4  0.0028485 ( 1.062e-03, 0.0076365)
# State 4 - State 3  0.0069290 ( 1.916e-03, 0.0250564)
# State 4 - State 4 -0.0074554 (-2.480e-02,-0.0022414)
# State 4 - State 5  0.0005264 ( 1.875e-04, 0.0014781)
# 
# -2 * log-likelihood:  563.9205
plot(init_msm, legend.pos=c(8, 1))

plot.prevalence.msm(init_msm, mintime=0, maxtime=1000)



help(optim)

# by covariates
anemia_msm <- msm(state ~ time, subject = ID, data = final_data, qmatrix = Q.crude, deathexact = 5, covariates = ~ dia_anemia)
# Transition intensities with hazard ratios for each covariate
# Baseline                           dia_anemia1                      
# State 1 - State 1 -0.0102703 (-2.966e+16,-3.556e-21)                                  
# State 1 - State 2  0.0102703 ( 3.556e-21, 2.966e+16) 8.614e+05 (2.348e-261,3.160e+272)
# State 2 - State 1  0.0010146 ( 3.525e-22, 2.921e+15) 1.108e+06 (3.014e-261,4.076e+272)
# State 2 - State 2 -0.0023441 (-2.293e+05,-2.396e-11)                                  
# State 2 - State 3  0.0013295 ( 8.191e-04, 2.158e-03) 1.532e+01 ( 6.682e-01, 3.512e+02)
# State 3 - State 2  0.0028118 ( 1.741e-03, 4.541e-03) 1.377e+01 ( 8.412e-01, 2.253e+02)
# State 3 - State 3 -0.0053771 (-8.706e-03,-3.321e-03)                                  
# State 3 - State 4  0.0025653 ( 1.082e-03, 6.080e-03) 6.134e+00 ( 3.922e-01, 9.594e+01)
# State 4 - State 3  0.0060857 ( 1.887e-03, 1.962e-02) 3.232e+00 ( 1.723e-01, 6.061e+01)
# State 4 - State 4 -0.0064053 (-2.579e-02,-1.591e-03)                                  
# State 4 - State 5  0.0003197 ( 1.794e-11, 5.697e+03) 9.625e-05 (2.639e-109,3.510e+100)
# 
# -2 * log-likelihood:  546.7068 

pid_msm <- msm(state ~ time, subject = ID, data = final_data, qmatrix = Q.crude, deathexact = 5, covariates = ~ dia_pid)
# Warning message:
#   In msm(state ~ time, subject = ID, data = final_data, qmatrix = Q.crude,  :
#            Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite.

vacc_msm <- msm(state ~ time, subject = ID, data = final_data, qmatrix = Q.crude, deathexact = 5, covariates = ~ vacc_c)
# Transition intensities with hazard ratios for each covariate
# Baseline                           vacc_c                        
# State 1 - State 1 -1.605e+01 (-3.725e+85,-6.918e-84)                               
# State 1 - State 2  1.605e+01 ( 6.918e-84, 3.725e+85) 0.013950 (1.787e-51,1.089e+47)
# State 2 - State 1  1.134e+00 ( 5.000e-85, 2.572e+84) 0.025447 (3.162e-51,2.048e+47)
# State 2 - State 2 -1.135e+00 (-2.098e+84,-6.144e-85)                               
# State 2 - State 3  1.199e-03 ( 7.501e-04, 1.915e-03) 0.651407 (2.332e-01,1.819e+00)
# State 3 - State 2  2.627e-03 ( 1.643e-03, 4.199e-03) 1.253981 (4.850e-01,3.242e+00)
# State 3 - State 3 -5.275e-03 (-8.805e-03,-3.160e-03)                               
# State 3 - State 4  2.648e-03 ( 1.067e-03, 6.569e-03) 1.722292 (4.636e-01,6.399e+00)
# State 4 - State 3  5.983e-03 ( 1.765e-03, 2.028e-02) 0.927769 (1.877e-01,4.585e+00)
# State 4 - State 4 -6.111e-03 (-2.105e-02,-1.774e-03)                               
# State 4 - State 5  1.278e-04 ( 2.778e-11, 5.883e+02) 0.002617 (2.411e-27,2.841e+21)
# 
# -2 * log-likelihood:  554.2041 

qmatrix.msm(vacc_msm, covariates=list(vacc_c=1))
# State 1                            State 2                            State 3                            State 4                           
# State 1 -7.310e-01 (-7.593e+48,-7.038e-50)  7.310e-01 ( 7.038e-50, 7.593e+48) 0                                  0                                 
# State 2  7.976e-02 ( 8.000e-51, 7.953e+47) -8.064e-02 (-2.351e+47,-2.766e-50)  8.792e-04 ( 3.422e-04, 2.259e-03) 0                                 
# State 3 0                                   3.094e-03 ( 1.331e-03, 7.192e-03) -7.017e-03 (-1.535e-02,-3.208e-03)  3.923e-03 ( 1.086e-03, 1.417e-02)
# State 4 0                                  0                                   5.668e-03 ( 1.046e-03, 3.070e-02) -5.669e-03 (-3.069e-02,-1.047e-03)
# State 5 0                                  0                                  0                                  0                                 
# State 5                           
# State 1 0                                 
# State 2 0                                 
# State 3 0                                 
# State 4  1.736e-06 ( 1.615e-30, 1.866e+18)
# State 5 0                                 
table(final_data$vacc_c)

## Optimization (Parameter fixed)
opt_msm <- msm(state ~ time, subject = ID, data = final_data, qmatrix = Q.crude, deathexact = 5, fixedpars = c(2,7))



