library(tidyverse)
library(psych)
library(dplyr)
library(tidyr)
library(data.table)
library(stringi)
library(stringr)

######################################################################################
##### create four groups based on biopsy: 1. CIN1 | 2. CIN2 | 3. CIN3 | 4. Neg  ######
##### absence of biopsy: 1. LSIL(CIN1) | 2. HSIL(over CIN2)                     ######
######################################################################################

# import files
finalCohort <- read.csv("C:/Users/hi/Documents/CINProgression/data/cohort_temp0108.csv", stringsAsFactors = F, fileEncoding = "euc-kr") # 5124

finalCohort$수술일자 <- as.Date(as.character(finalCohort$수술일자), format = "%Y%m%d")
finalCohort$최초진단일자 <- as.Date(as.character(finalCohort$최초진단일자), format = "%Y%m%d")
finalCohort$수술처방일자 <- as.Date(as.character(finalCohort$수술처방일자), format = "%Y%m%d")


##### Grouping #####
# variables based on biopsy result: CIN1, CIN2, CIN3, cancer, negative, hsil, nd(not detected)
# pap result variable
finalCohort$CIN1_p <- dplyr::case_when(finalCohort$results_pap %in% c("CIN1", "lsil", "ascus") ~ 1, TRUE ~ 0)

finalCohort$CIN3_p <- dplyr::case_when(finalCohort$results_pap %in% c("CIN3") ~ 1, TRUE ~ 0)

finalCohort$HSIL_p <- dplyr::case_when(finalCohort$results_pap %in% c("hsil") ~ 1, TRUE ~ 0)

###### pap negative
pap_results <- finalCohort[c("ID", "date", "results_pap", "results_hpv", "numofvs", "count")] %>% 
  arrange(ID, date) %>%
  group_by(ID) %>% 
  filter(!(is.na(results_pap))) %>% 
  mutate(consecutive_neg = lag(results_pap == "neg") & results_pap == "neg",
         Neg_pap = if_else(consecutive_neg == TRUE, 1, 0)) %>% 
  select(c(ID, date, Neg_pap))

finalCohort <- merge(finalCohort, pap_results, by = c("ID", "date"), all = T)

# negative: biopsy and consecutive neg of pap, except hpv result
finalCohort <- finalCohort %>%
  mutate(Neg = case_when(negative == 1 ~ 1,
                         is.na(negative) & Neg_pap == 1 ~ 1, TRUE ~ 0))

#### add biopsy data with chart review 1206 updated
library(openxlsx)
add <- read.xlsx("C:/Users/hi/Documents/CINProgression/data/preop_forR.xlsx", sheet = 1) %>% 
  mutate(date = as.character(convertToDate(date)))

finalCohort <- full_join(add, finalCohort, by = c("ID", "date")) %>%
  arrange(ID, date) %>% 
  mutate(CIN1 = ifelse(biopsy_add == "cin1"|biopsy_add == "lsil"|CIN1 == 1, 1, 0),
         CIN2 = ifelse(biopsy_add == "cin2"|CIN2 == 1, 1, 0),
         CIN3 = ifelse(biopsy_add == "cin3"|CIN3 == 1, 1, 0),
         nd = ifelse(biopsy_add == "chronic cervicitis"|nd == 1, 1, 0),
         results_pap = ifelse(is.na(pap_add), results_pap, pap_add),
         results_hpv = ifelse(is.na(hpv_add), results_hpv, "Pos"),
         hpvNum = ifelse(is.na(hpv_add), hpvNum, "hrNum"),
         h_18 = ifelse(is.na(h_18_add), h_18, h_18_add),
         h_58 = ifelse(is.na(h_58_add), h_58, h_58_add),
         h_66 = ifelse(is.na(h_66_add), h_66, h_66_add),
         multi_inf = ifelse(is.na(multi_inf_add), multi_inf, multi_inf_add)) %>% 
  select(-c(ends_with("_add"))) %>% 
  group_by(ID) %>% 
  fill(birth, .direction = "updown")

#### add hpv data with chart review 1218 updated
add_hpv <- read.xlsx("C:/Users/hi/Documents/CINProgression/data/no_hpv_1218.xlsx")

pttn_h <- ("16|18|26|31|33|35|39|45|51|52|53|56|58|59|66|68|69|70|73|82")
pttn_l <- ("6|11|30|32|34|40|42|43|44|54|55|61|62|67|72|74|81|83|84|87|90") # 72 added for this data

add_hpv <- add_hpv %>% 
  mutate(h_16_add = as.numeric(ifelse(str_detect(hpvNum_initial, "16"), 1, 0)),
         h_18_add = as.numeric(ifelse(str_detect(hpvNum_initial, "18"), 1, 0)),
         h_52_add = as.numeric(ifelse(str_detect(hpvNum_initial, "52"), 1, 0)),
         h_58_add = as.numeric(ifelse(str_detect(hpvNum_initial, "58"), 1, 0)),
         h_31_add = as.numeric(ifelse(str_detect(hpvNum_initial, "31"), 1, 0)),
         h_33_add = as.numeric(ifelse(str_detect(hpvNum_initial, "33"), 1, 0)),
         h_45_add = as.numeric(ifelse(str_detect(hpvNum_initial, "45"), 1, 0)),
         h_35_add = as.numeric(ifelse(str_detect(hpvNum_initial, "35"), 1, 0)),
         h_39_add = as.numeric(ifelse(str_detect(hpvNum_initial, "39"), 1, 0)),
         h_51_add = as.numeric(ifelse(str_detect(hpvNum_initial, "51"), 1, 0)),
         h_56_add = as.numeric(ifelse(str_detect(hpvNum_initial, "56"), 1, 0)),
         h_59_add = as.numeric(ifelse(str_detect(hpvNum_initial, "59"), 1, 0)),
         h_66_add = as.numeric(ifelse(str_detect(hpvNum_initial, "66"), 1, 0)),
         h_68_add = as.numeric(ifelse(str_detect(hpvNum_initial, "68"), 1, 0)),
         hpvNum_add = case_when(str_detect(hpvNum_initial, pttn_h) ~ "hrNum",
                            str_detect(hpvNum_initial, pttn_l) ~ "lrNum",
                            str_detect(hpvNum_initial, "low risk") ~ "lrNum")) %>% 
  mutate(multi_inf_add = as.numeric(rowSums(across(h_16_add:h_68_add), na.rm = T) >= 2))


finalCohort <- full_join(add_hpv, finalCohort, by = c("ID"), relationship = "many-to-many") %>%
  arrange(ID, date) %>% 
  mutate(hpvNum = ifelse(is.na(hpvNum_add), hpvNum, hpvNum_add),
         multi_inf = ifelse(is.na(multi_inf_add), multi_inf, multi_inf_add),
         h_16 = ifelse(is.na(h_16_add), h_16, h_16_add),
         h_18 = ifelse(is.na(h_18_add), h_18, h_18_add),
         h_52 = ifelse(is.na(h_52_add), h_52, h_52_add),
         h_58 = ifelse(is.na(h_58_add), h_58, h_58_add),
         h_31 = ifelse(is.na(h_31_add), h_31, h_31_add),
         h_33 = ifelse(is.na(h_33_add), h_33, h_33_add),
         h_45 = ifelse(is.na(h_45_add), h_45, h_45_add),
         h_35 = ifelse(is.na(h_35_add), h_35, h_35_add),
         h_39 = ifelse(is.na(h_39_add), h_39, h_39_add),
         h_51 = ifelse(is.na(h_51_add), h_51, h_51_add),
         h_56 = ifelse(is.na(h_56_add), h_56, h_56_add),
         h_59 = ifelse(is.na(h_59_add), h_59, h_59_add),
         h_66 = ifelse(is.na(h_66_add), h_66, h_66_add),
         h_68 = ifelse(is.na(h_68_add), h_68, h_68_add)) %>% 
  select(-ends_with("_add"), -hpvNum_initial)

###### covariates
# fill disease value since the diagnosed date
## aml  # 27
amlin <- finalCohort %>%
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(aml == 1) %>% 
  dplyr::mutate(earliestdd = min(최초진단일자)) %>% 
  arrange(earliestdd) %>% 
  slice(1)

amlin <- amlin$ID %>% unique()

aml <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  filter(ID %in% amlin) %>%  
  fill(최초진단일자) %>% 
  mutate(dia_aml = if_else(date >= 최초진단일자, 1, 0)) %>% 
  select(c(ID, date, dia_aml))

## anemia  # 142
anemiain <- finalCohort %>%
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(anemiaThrom == 1) %>% 
  dplyr::mutate(earliestdd = min(최초진단일자)) %>% 
  arrange(earliestdd) %>% 
  slice(1)

anemiain <- anemiain$ID

anemia <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  filter(ID %in% anemiain) %>%  
  fill(최초진단일자) %>% 
  mutate(dia_anemia = if_else(date >= 최초진단일자, 1, 0)) %>% 
  select(c(ID, date, dia_anemia))

## PID   # 2119
pidin <- finalCohort %>%
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(pid == 1) %>% 
  dplyr::mutate(earliestdd = min(최초진단일자)) %>% 
  arrange(earliestdd) %>% 
  slice(1)

pidin <- pidin$ID %>% unique()

pid <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  filter(ID %in% pidin) %>%  
  fill(최초진단일자) %>% 
  mutate(dia_pid = if_else(date >= 최초진단일자, 1, 0)) %>% 
  select(c(ID, date, dia_pid))

##### vaccination ID
# vaccination records from EMR charts updated AUG 11
gar_id <- ("R000000239|R000002229|R000003617|R000003834|R000004961|R000004964|R000006482|R000008597|R000010600|R000011826|R000020193|R000020538|R000022206|R000023310|R000027094|R000027771|R000000356|R000000525|R000000683|R000001667|R000001785|R000001817|R000002268|R000002425|R000002454|R000002646|R000002750|R000002938|R000003009|R000024143|R000000294|R000000708|R000001393|R000001680|R000001801|R000002404|R000002644|R000002894|R000011034|R000028950")

cer_id <- ("R000003834|R000008750|R000008980|R000016494|R000023807")

other_id <-("R000001098|R000003614|R000003634|R000007247|R000008759|R000014748|R000016465|R000021758|R000022159|R000024276|R000027888|R000000081|R000000581|R000001215|R000001595|R000002079|R000002100|R000002175|R000002639|R000003200|R000001533|R000002653|R000006281|R000013865|R000015597|R000024860")

no_id <- ("R000000198|R000003520|R000003555|R000003558|R000003839|R000003870|R000003902|R000003967|R000004128|R000004157|R000004489|R000004491|R000005817|R000008012|R000008877|R00001039|R000012489|R000012635|R000012734|R000013879|R000014876|R000015819|R000016648|R000016755|R000019146|R000019773|R000020564|R000021078|R000024042|R000025451|R000026919|R000027423|R00027432|R000028187|R000000608|R000000655|R000001196|R000001947|R000002181|R000002327|R000002446|R000002682|R000003338|R000003379|R000003383|R000003753|R000016166|R000027510|R000028418|R000000091|R000000297|R000000432|R000000464|R000000737|R000000983|R000001194|R000001213|R000001218|R000001276|R000001451|R000001542|R000001600|R000001628|R000001656|R000001839|R000001895|R000001926|R000002396|R000002413|R000002435|R000002604|R000002715|R000002801|R000003138|R000003290|R000003376|R000003539|R000003817|R000004560|R000007294|R000012092|R000018193|R000018767|R000020884|R000021645|R000021954|R000023099|R000023658|R000023826|R000030157")

addVc <- finalCohort %>% 
  ungroup() %>% 
  mutate(add_g = as.integer(str_detect(finalCohort$ID, str_c(gar_id, collapse = "|"))),
         add_o = as.integer(str_detect(finalCohort$ID, str_c(other_id, collapse = "|"))),
         add_c = as.integer(str_detect(finalCohort$ID, str_c(cer_id, collapse = "|"))),
         novc = as.integer((str_detect(finalCohort$ID, str_c(no_id, collapse = "|"))))) %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  slice(1) %>% 
  select(ID, add_g, add_o, add_g, add_c, novc)

add_g <- addVc %>% 
  filter(add_g == 1)

add_c <- addVc %>% 
  filter(add_c == 1)

add_o <- addVc %>% 
  filter(add_o == 1)

novc <- addVc %>% 
  filter(novc == 1)

## gardasil
garid <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(gardasil == 1) %>% 
  slice(1) %>% 
  select(ID, gardasil)

gard_id <- c(garid$ID, add_g$ID)

gar <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  filter(ID %in% gard_id) %>%
  mutate(vac_gar = 1) %>% 
  select(c(ID, date, vac_gar))

## cervarix
cervarid <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(cervarix == 1) %>% 
  slice(1) %>% 
  select(ID, cervarix)

cerv_id <- c(cervarid$ID, add_c$ID)

cer <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  filter(ID %in% cerv_id) %>%
  mutate(vac_cer = 1) %>% 
  select(c(ID, date, vac_cer))

### other vaccine
othervc <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(ID %in% add_o$ID) %>% 
  mutate(othervc = 1) %>% 
  select(ID, date, othervc)

### no vaccination
novc <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(ID %in% novc$ID) %>% 
  mutate(novc = 1) %>% 
  select(ID, date, novc)

### cancer  # 291
cancer <- finalCohort %>%
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(cancer_d == 1) %>% 
  dplyr::mutate(earliestdd = min(최초진단일자)) %>% 
  arrange(earliestdd) %>% 
  slice(1)

cancerin <- cancer$ID

cancer <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  filter(ID %in% cancerin) %>%
  fill(최초진단일자) %>% 
  mutate(dia_cancer = if_else(date >= 최초진단일자, 1, 0)) %>% 
  select(c(ID, date, dia_cancer))

### combine cohort data with covariates
finalCohort <- finalCohort %>%
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  left_join(aml, by = c("ID", "date")) %>% 
  left_join(anemia, by = c("ID", "date")) %>% 
  left_join(pid, by = c("ID", "date")) %>% 
  left_join(gar, by = c("ID", "date")) %>% 
  left_join(cer, by = c("ID", "date")) %>% 
  left_join(othervc, by = c("ID", "date")) %>% 
  left_join(novc, by = c("ID", "date")) %>% 
  left_join(cancer, by = c("ID", "date"))

### separate the cohort into biopsy cohort and pap cohort
# get rid of unnecessary rows (only left biopsy, hpv, pap result row)
finalCohort <- finalCohort %>% 
  ungroup() %>% 
  filter(!(is.na(CIN1) & is.na(CIN2) & is.na(CIN3) & is.na(cancer) & is.na(negative) & is.na(hsil) & is.na(nd) & is.na(results_pap) & is.na(results_hpv))) %>% 
  select(-c("주소", "latestD...earliestD", "hbv":"renalf", "gardasil", "cervarix", "results_hpv")) %>% 
  mutate_at(vars(CIN1, CIN2, CIN3, hsil, CIN1_p, CIN3_p, HSIL_p), as.numeric) %>% 
  mutate_at(vars(CIN1, CIN2, CIN3, hsil, CIN1_p, CIN3_p, HSIL_p), ~replace_na(., 0)) %>% 
  mutate_at(vars(CIN1, CIN2, CIN3, hsil, CIN1_p, CIN3_p, HSIL_p), as.factor)

### grouping
finalCohort <- finalCohort %>% 
  mutate(G1 = case_when(CIN1 == 1 & CIN2 == 0 & CIN3 == 0 ~ 1, TRUE ~ 0),
         G2 = case_when((CIN1 == 0 & CIN2 == 1 & CIN3 == 0)|(CIN1 == 1 & CIN2 == 1 & CIN3 == 0) ~ 1, TRUE ~ 0),
         G3 = case_when((CIN1 == 0 & CIN2 == 0 & CIN3 == 1)|(CIN1 == 0 & CIN2 == 1 & CIN3 == 1)|(CIN1 == 1 & CIN2 == 0 & CIN3 == 1)|(CIN1 == 1 & CIN2 == 1 & CIN3 == 1) ~ 1, TRUE ~ 0),
         G_hsil = case_when(G1 == 0 & G2 == 0 & G3 == 0 & hsil == 1 ~ 1, TRUE ~ 0),
         G1_p = case_when(G1 == 0 & G2 == 0 & G3 == 0 & hsil == 0 & CIN1_p == 1 & HSIL_p == 0 ~ 1, TRUE ~ 0),
         G_p_hsil = case_when(G1 == 0 & G2 == 0 & G3 == 0 & hsil == 0 & G1_p == 0 & HSIL_p == 1 ~ 1, TRUE ~ 0))

### find index date
# drop the existing count variables
finalCohort <- finalCohort %>% 
  select(-c("numofvs", "count")) %>% 
  distinct(ID, date, .keep_all = T)

# add up by ID
finalCohort$numofvs <- ave(integer(nrow(finalCohort)), finalCohort$ID, FUN = seq_along)
finalCohort <- finalCohort %>%
  group_by(ID) %>% 
  mutate(count = n_distinct(date))

# find out index date using first positive(based on biopsy)  #pap only는 뒤에 다시
finalCohort <- finalCohort %>%
  group_by(ID) %>% 
  mutate(first_occ = row_number() == which(G1 == 1 | G2 == 1 | G3 == 1)[1])   # G_hsil == 1 일단 제외...17 ids

idx_biopsy <- finalCohort %>% 
  group_by(ID) %>%
  filter(first_occ == TRUE) %>% 
  arrange(ID, date) %>% 
  slice(1) %>% 
  mutate(idx = date) %>% 
  select(ID, idx)

groupedID <- idx_biopsy$ID %>% unique() # 2349

####### exclusion: no information for grouping 2775
noinfo <- finalCohort %>% 
  filter(!(ID %in% groupedID))

noinfoID <- noinfo$ID %>% unique()

# join the index data into finalcohort
df_idx <- finalCohort %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(ID %in% groupedID) %>% 
  left_join(idx_biopsy, by = "ID") %>% 
  fill(idx) %>% 
  mutate(entry = ifelse(idx == date, 1, 0))

# remove rows before entry  ---> filter는 각 그룹 분석할 때..hpv 결과 가져와야하기 때문에  %>% filter(date >= idx)
cohort_full <- df_idx %>%
  group_by(ID)

# imd operation variable (treatment in one year)
# hysterectomy
hysin <- cohort_full %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(hys == 1) %>%
  rename(hys_date = 수술일자) %>% 
  mutate(hys_df = difftime(hys_date, idx, units = "days")) %>% 
  select(c(ID, date, entry, idx, hys_date, hys, hys_df))

hys_1yr <- hysin %>%
  filter(hys_df < 366)

# True leep
true_leep <- cohort_full %>%
  filter(leep == 1) %>% 
  rename(leep_date = 수술처방일자) %>% 
  mutate_at(vars(CIN1, CIN2, CIN3, cancer, negative, hsil, nd), ~replace_na(., 0)) %>% 
  mutate(leep_df = difftime(leep_date, idx, units = "days"),
         leep_duration = difftime(date, leep_date, units = "days"),
         biopsy_leep = ifelse(CIN3 == 1 , "CIN3",
                              ifelse(CIN2 == 1, "CIN2",
                                     ifelse(CIN1 == 1, "CIN1",
                                            ifelse(cancer == 1, "cancer",
                                                   ifelse(negative == 1, "negative",
                                                          ifelse(hsil == 1, "hsil",
                                                                 ifelse(nd == 1, "not detected", NA)))))))) %>%
  filter(leep_df < 366) %>% # 1yr from idx date
  slice(which.min(abs(leep_duration))) %>% 
  select(ID, date, CIN1, CIN2, CIN3, cancer, negative, hsil, nd, biopsy_leep, leep, leep_date) %>% 
  filter(!is.na(biopsy_leep)) %>%   # 실제 수행되지 않은 군 제외 
  select(ID, biopsy_leep, leep, leep_date)

imd_operation <- true_leep %>% 
  bind_rows(hys_1yr) %>% 
  mutate(imd_op = 1) %>% 
  select(ID, imd_op, hys_date, leep_date, hys, leep, biopsy_leep)

# merge the dataset
cohort_full <- cohort_full %>%
  select(-c(hys, leep)) %>% 
  left_join(imd_operation, by = c("ID"))

###### exclusion #######
# 1. cancer
#### cancer diagnosed in 90 days from the index date     *** 기존 180일에서 90일로 기준 수정 18 Oct
cancer_d <- cohort_full %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  filter(cancer_d == 1) %>% 
  mutate(dia_dft = difftime(최초진단일자, idx, units = "days")) %>% 
  select(c(ID, date, idx, 최초진단일자, cancer_d, dia_dft))

cancer_ex <- cancer_d %>% 
  filter(dia_dft < 91)

# 제외할 id   N = 88
cancer_ex_id <- cancer_ex$ID %>% unique()   

cohort_full <- cohort_full %>% 
  filter(!(ID %in% cancer_ex_id))

# 2. visiting duration (since the idx date) is less than 365days  
less1yr <- cohort_full %>%
  arrange(ID, date) %>% 
  group_by(ID) %>% 
  mutate(latestD = max(date),
         duration = difftime(latestD, idx, units = "days")) %>% 
  select(ID, date, latestD, idx, duration) %>% 
  filter(duration < 366)

less1yr_id <- less1yr$ID %>% unique() # N = 1021

cohort_full <- cohort_full %>%    
  filter(!(ID %in% less1yr_id))

cohort_id <- cohort_full$ID %>% unique() # 1240

# after cleaning up, make a 1st visit & last visit variable again
cohort_full <- cohort_full %>% 
  mutate(date = as.Date(date)) %>% 
  select(-c("numofvs", "count")) %>% 
  distinct(ID, date, .keep_all = T)

# add up by ID
cohort_full$numofvs <- ave(integer(nrow(cohort_full)), cohort_full$ID, FUN = seq_along)

cohort_full <- cohort_full %>% 
  group_by(ID) %>% 
  mutate(count = n_distinct(date))

# time difference btw visits
cohort_full <- cohort_full %>%
  arrange(ID, date) %>% 
  group_by(ID) %>%
  mutate(dftime = date - lag(date))

# time differ btw entry date and last visit
cohort_full <- cohort_full %>%
  arrange(ID, date) %>% 
  group_by(ID) %>%
  mutate(latestD = max(date),
         duration = difftime(latestD, idx, units = "days"))

# data has been cleaned, create a age using dob and date of visit
cohort_full$age <- trunc((cohort_full$birth %--% cohort_full$date) / years(1))
cohort_full$age <- as.numeric(cohort_full$age)

### finalize the cohort data
# make a csv file (full cohort)
write.csv(cohort_full, "C:/Users/hi/Documents/CINProgression/data/cohort_0108.csv", row.names = F, fileEncoding = "euc-kr")

# g1
gr1 <- cohort_full %>% 
  arrange(ID, date) %>% 
  filter(entry == 1) %>% 
  filter(G1 == 1) %>% 
  select(ID, date, idx)

group1_id <- gr1$ID %>% unique() # 628

group1 <- cohort_full %>% 
  filter(ID %in% gr1$ID)

# g2
gr2 <- cohort_full %>% 
  arrange(ID, date) %>% 
  filter(entry == 1) %>% 
  filter(G2 == 1) %>% 
  select(ID, date, idx)

group2_id <- gr2$ID %>% unique() # 340

group2 <- finalCohort %>% 
  filter(ID %in% gr2$ID)

# g3
gr3 <- cohort_full %>% 
  arrange(ID, date) %>% 
  filter(entry == 1) %>% 
  filter(G3 == 1) %>% 
  select(ID, date, idx)

group3_id <- gr3$ID %>% unique() # 272

group3 <- finalCohort %>% 
  filter(ID %in% gr3$ID)

#### save the csv file by group
CIN1 <- cohort_full %>%
  arrange(ID, date) %>% 
  filter(ID %in% group1_id)

write.csv(CIN1, "C:/Users/hi/Documents/CINProgression/data/CIN1_0108.csv", row.names = F, fileEncoding = "euc-kr")


CIN2 <- cohort_full %>%
  arrange(ID, date) %>% 
  filter(ID %in% group2_id)

write.csv(CIN2, "C:/Users/hi/Documents/CINProgression/data/CIN2_0108.csv", row.names = F, fileEncoding = "euc-kr")


CIN3 <- cohort_full %>%
  arrange(ID, date) %>% 
  filter(ID %in% group3_id)

write.csv(CIN3, "C:/Users/hi/Documents/CINProgression/data/CIN3_0108.csv", row.names = F, fileEncoding = "euc-kr")
