library(tidyverse)
library(psych)
library(dplyr)
library(tidyr)
library(data.table)
library(stringi)
library(stringr)
library(textclean)

# import files
cohort_c <- read.csv("C:/Users/hi/Documents/CINProgression/data/코호트.csv", stringsAsFactors = F, fileEncoding = "euc-kr")
cohort_t <- read.csv("C:/Users/hi/Documents/CINProgression/data/병리검사.csv", stringsAsFactors = F, fileEncoding = "euc-kr")
cohort_r <- read.csv("C:/Users/hi/Documents/CINProgression/data/수술기록.csv", stringsAsFactors = F, fileEncoding = "euc-kr")
cohort_s <- read.csv("C:/Users/hi/Documents/CINProgression/data/수술처방.csv", stringsAsFactors = F, fileEncoding = "euc-kr")
cohort_d <- read.csv("C:/Users/hi/Documents/CINProgression/data/진단정보.csv", stringsAsFactors = F, fileEncoding = "euc-kr")
cohort_p <- read.csv("C:/Users/hi/Documents/CINProgression/data/처방정보.csv", stringsAsFactors = F, fileEncoding = "euc-kr")

# 일자 -> as.Date
cohort_c$birth <- as.Date(as.character(cohort_c$생년월), format = "%Y%m%d")
cohort_c <- cohort_c[-c(2, 3, 6)]
cohort_c <- rename(cohort_c, c("ID" = "연구번호"))  # 5124

cohort_t$date <- as.Date(as.character(cohort_t$처방일자), format = "%Y%m%d")
cohort_t <- cohort_t[-c(5, 7, 9)]
cohort_t <- rename(cohort_t, c("ID" = "연구번호"))

cohort_r$date <- as.Date(as.character(cohort_r$수술일자), format = "%Y%m%d")
cohort_s$date <- as.Date(as.character(cohort_s$수술처방일자), format = "%Y%m%d")
cohort_d$date <- as.Date(as.character(cohort_d$진단일자), format = "%Y%m%d")
cohort_p$date <- as.Date(as.character(cohort_p$처방일자), format = "%Y%m%d")

cohort_r <- rename(cohort_r, c("ID" = "연구번호"))
cohort_s <- rename(cohort_s, c("ID" = "연구번호"))
cohort_d <- rename(cohort_d, c("ID" = "연구번호"))
cohort_p <- rename(cohort_p, c("ID" = "연구번호"))

#####  병리검사 #####
# *******PGC001 병리검사구분 == "분자병리" is HPV, "세포병리" is pap
# hpv
# HPV, hybrid capture PMO07001, PMO07002 제외
hpv_code <- paste("PGC001|PMO03001|PMO03002|PMO03003|PMO03004|PMO03005|PMO08006|PMO08007|PMO12039|PMO12040|PMO12040A|PMO12049|PMO12050|PMO12075|PMO12075A|PMO12078|PMO12078A", collapse = "|")
c_hpv <- cohort_t[grep(paste(hpv_code), cohort_t$처방코드),] %>% 
  filter(병리검사구분 != "세포병리") # 처방코드 PGC001 중 pap 제외

# pap
pap_code <- paste("PGC001|PCY10001|PCY10003|PCY10003A|PCY10004|PCY10004A|PCY10005|PCY10005A|PCY10006A|PCY10008|PCY10010|POU20001|POU21001|POU21002", collapse = "|")
c_pap <- cohort_t[grep(paste(pap_code), cohort_t$처방코드),] %>% 
  filter(병리검사구분 != "분자병리") # 처방코드 PGC001 중 hpv 제외

# biopsy
cerv_code <- cohort_t[!((cohort_t$처방코드 %in% c_pap$처방코드)|(cohort_t$처방코드 %in% c_hpv$처방코드)|(cohort_t$처방코드 %in% "PMO07001")|(cohort_t$처방코드 %in% "PMO07002")), ]
cerv_code <- cerv_code$처방코드 %>% unique()

# ---------------------
# 1. pap result
# trimming and cleaning
clean_result <- c_pap$판독결과 %>% 
  str_squish() %>% # remove heading & tailing white space
  str_to_lower()

c_pap$result = sub("comments.*", "", clean_result) %>% 
  str_replace_all("[[:punct:]]+", "") %>% 
  str_squish()

# outcomes: ascus, lsil, hsil, asc-h
# Atypical squamous cells cannot exclude HSIL (ASC-H)
# Severity: Reactive change, ASCUS, LSIL, ASC-H, Atypical glandular change, Adenocarcinoma in citu, Squamous carcinoma insitu, Cervical cancer
# hsil > lsil > atypical squamous cells > RCC=negative
# 1117 updated
# extract results
tidy_pap <- c_pap %>% 
  mutate(results_pap = case_when(str_detect(result, "highgrade squamous intraepithelial lesion") ~ "hsil",
                                 str_detect(result, "lowgrade squamous intraepithelial lesion") ~ "lsil",
                                 str_detect(result, "atypical squamous cells cannot exclude hsil") ~ "asc-h",
                                 str_detect(result, "asch") ~ "asc-h",
                                 str_detect(result, "atypical squamous cells of undetermined significance") ~ "ascus",
                                 str_detect(result, "ascus") ~ "ascus",
                                 str_detect(result, "endocervical adenocarcinoma in situ") ~ "endocervical adenocarcinoma in situ",
                                 str_detect(result, "atypical squamous cells") ~ "atypical squamous cells",
                                 str_detect(result, "atypical endocervical cells") ~ "atypical endocervical cells",
                                 str_detect(result, "favor neoplastic") ~ "favor neoplastic",
                                 str_detect(result, "atypical glandular cells") ~ "atypical glandular cells",
                                 str_detect(result, "reactive changes|reactive cellular changes") ~ "reactive",
                                 str_detect(result, "squamous cell carcinoma") ~ "squamous cell carcinoma",
                                 str_detect(result, "unsatisfactory") ~ "unsatisfactory",
                                 str_detect(result, "negative") ~ "neg",
                                 str_detect(result, "absence") ~ "absence")) %>% 
  select(ID, date, results_pap)

# ---------------------
# 2. HPV
# trimming and cleaning
c_hpv <- c_hpv %>%
  group_by(ID, date) %>% arrange(date) %>% filter(row_number() == 1)  # 중복데이터 제거(업로드 상의 문제로 중복된 데이터 생성된 것으로 보임)

hpv_clean_result <- c_hpv$판독결과 %>% 
  str_squish() %>% # remove heading & tailing white space
  str_to_lower()

c_hpv$result = sub("comments.+$|comment].+$", "", hpv_clean_result) # remove unnecessary text

# 오타 수정
c_hpv$result <- gsub("negaitve", "negative", c_hpv$result)
c_hpv$result <- gsub("postive", "positive", c_hpv$result)


hpv_full <- c_hpv %>%
  dplyr::mutate(results_hpv = case_when(stri_detect_fixed(result, "(-)") ~ "Neg",
                                    stri_detect_fixed(result, "(+)") ~ "Pos",
                                    stri_detect_fixed(result, "(++)") ~ "Pos",
                                    stri_detect_fixed(result, "(+++)") ~ "Pos",
                                    stri_detect_fixed(result, "not detected") ~ "Neg",
                                    stri_detect_fixed(result, "positive") ~ "Pos",
                                    stri_detect_fixed(result, "negative") ~ "Neg",
                                    stri_detect_fixed(result, "감염되었") ~ "Pos",
                                    stri_detect_fixed(result, "음성") ~ "Neg"))

hpv_full["results_hpv"][is.na(hpv_full["results_hpv"])] <- "invalid"  ### [결과: 이미지참조]로 알 수 없음 or invalid

# trimming
hpv_full$result <- hpv_full$result %>% 
  str_squish() %>% 
  str_replace_all("\\(\\+\\+\\+\\)|\\+\\+\\+|\\(\\+\\+\\)|\\+\\+|\\(\\+\\)|\\+", "POS") %>% 
  str_replace_all("[[:punct:]]+", "") %>% 
  str_squish() %>% 
  str_replace(" ", "")

hpv_full$result <- gsub(" ", "", hpv_full$result) ### remove all white space
hpv_full$result <- gsub("임상적의의.+$|14종의.+$|.*\\진단실", "", hpv_full$result) # remove unnecessary text

# hpv result patterns
pttn_p <- "positive[0-9]+|highriskgroup[0-9]+|[0-9]+번|subtype[0-9]+group[0-9]+|pcr[0-9]+|result[0-9]+"

pttn_o <- ("positiveother|othertype|positive6other|분류되지않은")

pttn_n <- ("notdetected|negative|음성")

pttn_pp <- ("[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS$|[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS$|[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS$|[0-9]{1,2}?POS[0-9]{1,2}?POS[0-9]{1,2}?POS$|[0-9]{1,2}?POS[0-9]{1,2}?POS$|[0-9]+POS$|pcr[0-9]+POS|result[0-9]+POSadditional") ######## 101.5 등 필요없는 결과는 추출한 후 지우기

pttn_hrp <- "positivep[0-9]+p[0-9]+p[0-9]|positivep[0-9]+p[0-9]+|positivep[0-9]+"


# ----------------------------
# get specific result by patterns *** separate dataset (num_p, num_o) / unnest (pttn_pp) / unnest (pttn_hrp)
# 1) num
hpv <- hpv_full[-c(2:5)] %>% 
  mutate(result = gsub("hpvlowrisktypepositive", "", result),
         num_p = str_extract(result, str_c(pttn_p, collapse = "|")),
         num_p = gsub("p13358|p2565966|p3353968", "", num_p), # ignor <- ("p13358|p2565966|p3353968")  # remove these patterns from num_e
         num_o = str_extract(result, str_c(pttn_o, collapse = "|")),
         num_n = str_extract(result, str_c(pttn_n, collapse = "|")),
         pos = str_extract_all(result, str_c(pttn_pp, collapse = "|")),
         pos = gsub("\\d{3,}POS", NA, pos),  # get rid of unnecessary results: formatted like 142.32 and text
         pos = gsub("pcr|additional|result", "", pos),
         hrp = str_extract_all(result, str_c(pttn_hrp, collapse = "|")))

# 2) pp
hpv_pp <- hpv %>% unnest(pos)

# 3) hrp
hpv_hrp <- hpv %>% unnest(hrp)

# count hpv num 16|18|52|58|31|33|45|35|39|51|56|59|66|68
count <- hpv %>% 
  filter(results_hpv == "Pos") %>% 
  mutate(h_16 = as.numeric(ifelse(str_detect(num_p, "16") | str_detect(num_o, "16"), 1, 0)),
         h_18 = as.numeric(ifelse(str_detect(num_p, "18") | str_detect(num_o, "18"), 1, 0)),
         h_52 = as.numeric(ifelse(str_detect(num_p, "52") | str_detect(num_o, "52"), 1, 0)),
         h_58 = as.numeric(ifelse(str_detect(num_p, "58") | str_detect(num_o, "58"), 1, 0)),
         h_31 = as.numeric(ifelse(str_detect(num_p, "31") | str_detect(num_o, "31"), 1, 0)),
         h_33 = as.numeric(ifelse(str_detect(num_p, "33") | str_detect(num_o, "33"), 1, 0)),
         h_45 = as.numeric(ifelse(str_detect(num_p, "45") | str_detect(num_o, "45"), 1, 0)),
         h_35 = as.numeric(ifelse(str_detect(num_p, "35") | str_detect(num_o, "35"), 1, 0)),
         h_39 = as.numeric(ifelse(str_detect(num_p, "39") | str_detect(num_o, "39"), 1, 0)),
         h_51 = as.numeric(ifelse(str_detect(num_p, "51") | str_detect(num_o, "51"), 1, 0)),
         h_56 = as.numeric(ifelse(str_detect(num_p, "56") | str_detect(num_o, "56"), 1, 0)),
         h_59 = as.numeric(ifelse(str_detect(num_p, "59") | str_detect(num_o, "59"), 1, 0)),
         h_66 = as.numeric(ifelse(str_detect(num_p, "66") | str_detect(num_o, "66"), 1, 0)),
         h_68 = as.numeric(ifelse(str_detect(num_p, "68") | str_detect(num_o, "68"), 1, 0))) %>% 
  select(ID, date, starts_with("h_")) %>%
  mutate(multi_inf = as.numeric(rowSums(across(h_16:h_68), na.rm = T) >= 2))

hpvsingle_inf <- count %>% 
  filter(multi_inf == 0) %>% 
  select(-multi_inf)

hpvmulti_inf <- count %>% 
  filter(multi_inf == 1) %>% 
  select(ID, date, multi_inf)

multi_16 <- count %>%
  filter(multi_inf == 1) %>% 
  filter(h_16 == 1) %>%
  mutate(multi_16 = 1) %>% 
  select(ID, date, multi_16)

multi_18 <- count %>%
  filter(multi_inf == 1) %>% 
  filter(h_18 == 1) %>%
  mutate(multi_18 = 1) %>% 
  select(ID, date, multi_18)

multi_1618 <- count %>%      # 16&18
  filter(multi_inf == 1) %>% 
  filter((h_16 == 1)&(h_18 == 1)) %>%
  mutate(multi_1618 = 1) %>% 
  select(ID, date, multi_1618)

# hpv positive number grouping: high risk / low risk
pttn_h <- ("16|18|26|31|33|35|39|45|51|52|53|56|58|59|66|68|69|70|73|82|p1|p2|p3")
pttn_l <- ("6|11|30|32|34|40|42|43|44|54|55|61|62|67|72|74|81|83|84|87|90|other|positiveother|분류")

# extract hpv_num
hpv_num <- hpv %>% 
  mutate(hr_p = as.factor(ifelse(results_hpv == "Pos" & str_detect(num_p, pttn_h), "1", 0)),
         lr_o = as.factor(ifelse(results_hpv == "Pos" & str_detect(num_o, pttn_l), "1", 0)),
         lr_p = as.factor(ifelse(results_hpv == "Pos" & str_detect(num_p, pttn_l), "1", 0)),
         hr_1_n = as.integer(results_hpv == "Pos" & (str_detect(num_p, str_c("16|18", collapse = "|")))),
         hr_2_n = as.integer(results_hpv == "Pos" & (str_detect(num_p, str_c("52|58", collapse = "|")))),
         hr_3_n = as.integer(results_hpv == "Pos" & (str_detect(num_p, str_c("31|33|45", collapse = "|")))),
         hr_4_n = as.integer(results_hpv == "Pos" & (str_detect(num_p, str_c("35|39|51|56|59|66|68", collapse = "|"))))) %>% 
  select(c(ID, date, hr_p:hr_4_n))

# hpv_pp
hpv_pp <- hpv_pp %>% 
  mutate(hr_pp = as.factor(ifelse(results_hpv == "Pos" & str_detect(pos, pttn_h), "1", 0)),
         lr_pp = as.factor(ifelse(results_hpv == "Pos" & str_detect(pos, pttn_l), "1", 0)),
         hr_1_p = as.integer(results_hpv == "Pos" & (str_detect(pos, str_c("16|18", collapse = "|")))),
         hr_2_p = as.integer(results_hpv == "Pos" & (str_detect(pos, str_c("52|58", collapse = "|")))),
         hr_3_p = as.integer(results_hpv == "Pos" & (str_detect(pos, str_c("31|33|45", collapse = "|")))),
         hr_4_p = as.integer(results_hpv == "Pos" & (str_detect(pos, str_c("35|39|51|56|59|66|68", collapse = "|"))))) %>% 
  select(c(ID, date, hr_pp:hr_4_p))

# hpv_hrp
hpv_hrp <- hpv_hrp %>% 
  mutate(hr_hrp = as.factor(ifelse(results_hpv =="Pos" & str_detect(hrp, pttn_h), "1", 0)),
         lr_hrp = as.factor(ifelse(results_hpv =="Pos" & str_detect(hrp, pttn_l), "1", 0)),
         hr_1_h = as.integer(results_hpv == "Pos" & (str_detect(hrp, str_c("16|18", collapse = "|")))),
         hr_2_h = as.integer(results_hpv == "Pos" & (str_detect(hrp, "52"))),
         hr_3_h = as.integer(results_hpv == "Pos" & (str_detect(hrp, str_c("31|45", collapse = "|")))),
         hr_4_h = as.integer(results_hpv == "Pos" & (str_detect(hrp, str_c("51|p3", collapse = "|"))))) %>% 
  select(c(ID, date, hr_hrp:hr_4_h))

# combine'em
list_df <- list(hpv_num, hpv_pp, hpv_hrp)
combine_hpvn <- list_df %>% reduce(full_join, by = c("ID", "date"))

# make variables: HRgroup, LRgroup, negative, invalid, hpvNum
tidy_hpv <- hpv %>%
  select(ID, date, results_hpv) %>%
  left_join(combine_hpvn, by = c("ID", "date")) %>% 
  mutate(hrN = as.factor(ifelse((hr_p == 1)|(hr_pp == 1)|(hr_hrp == 1), "1", 0)),
         lrN = as.factor(ifelse((lr_o == 1)|(lr_p == 1)|(lr_pp == 1)|(lr_hrp == 1), "1", 0)),
         Neg = as.factor(ifelse(results_hpv == "Neg", "1", 0)),
         inval = as.factor(ifelse((results_hpv == "invalid"), "1", 0)),
         hpvNum = case_when(lrN == 1 & is.na(hrN) ~ "lrNum",
                            hrN == 1 ~ "hrNum",
                            Neg == 1 ~ "Neg",
                            inval == 1 ~ "invalid"))

#### HR positive number grouping - last step
hgr <- tidy_hpv %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  filter(results_hpv == "Pos") %>% 
  select(ID, date, hr_1_n, hr_2_n, hr_3_n, hr_4_n, hr_1_p, hr_2_p, hr_3_p, hr_4_p, hr_1_h, hr_2_h, hr_3_h, hr_4_h) %>% 
  filter(!((is.na(hr_1_n) & is.na(hr_2_n) & is.na(hr_3_n) & is.na(hr_4_n) & (is.na(hr_1_p) & is.na(hr_2_p) & is.na(hr_3_p) & is.na(hr_4_p) & (is.na(hr_1_h) & is.na(hr_2_h) & is.na(hr_3_h) & is.na(hr_4_h)))))) %>%
  mutate(multi = as.numeric(rowSums(across(hr_1_n:hr_4_h), na.rm = T) >= 2))

HRgr <- hgr %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  mutate(hr_1 = as.factor(ifelse((hr_1_n == 1|hr_1_h == 1|hr_1_p == 1), "1", 0)),
         hr_2 = as.factor(ifelse((hr_2_n == 1|hr_2_h == 1|hr_2_p == 1), "1", 0)),
         hr_3 = as.factor(ifelse((hr_3_n == 1|hr_3_h == 1|hr_3_p == 1), "1", 0)),
         hr_4 = as.factor(ifelse((hr_4_n == 1|hr_4_h == 1|hr_4_p == 1), "1", 0))) %>%
  mutate_at(vars(hr_1:hr_4), as.numeric) %>% 
  mutate_at(vars(hr_1:hr_4), ~replace_na(., 0)) %>% 
  mutate(multiple = as.numeric(rowSums(across(hr_1:hr_4), na.rm = T) >= 2)) %>% 
  select(c(ID, date, hr_1:multiple)) %>% 
  filter(!(hr_1 == 0 & hr_2 == 0 & hr_3 == 0 & hr_4 == 0))

# 우선순위 적용
HRgr_gr <- HRgr %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  filter(multiple == 1) %>% 
  mutate(hr_1_gr = as.factor(ifelse(hr_1 == 1, "1", 0)),
         hr_2_gr = as.factor(ifelse((hr_1 == 0) & (hr_2 == 1), "1", 0)),
         hr_3_gr = as.factor(ifelse((hr_1 == 0) & (hr_2 == 0) & (hr_3 == 1), "1", 0)),
         hr_4_gr = as.factor(ifelse((hr_1 == 0) & (hr_2 == 0) & (hr_3 == 0) & (hr_4 == 1), "1", 0))) %>% 
  select(c(ID, date, hr_1_gr:hr_4_gr))

HRgroup <- dplyr::left_join(HRgr, HRgr_gr, by = c("ID", "date"))

# multiple인 것도 합치기
HRgroup <- HRgroup %>% 
  arrange(ID, date) %>% 
  group_by(ID) %>%
  mutate(hr1 = as.factor(ifelse((hr_1 == 1)|(hr_1_gr == 1), "1", 0)),
         hr2 = as.factor(ifelse(((hr_1 == 0)|(hr_1_gr == 0)) & ((hr_2 == 1)|(hr_2_gr == 1)), "1", 0)),
         hr3 = as.factor(ifelse(((hr_1 == 0)|(hr_1_gr == 0)) & ((hr_2 == 0)|(hr_2_gr == 0)) & ((hr_3 == 1)|(hr_3_gr == 1)), "1", 0)),
         hr4 = as.factor(ifelse(((hr_1 == 0)|(hr_1_gr == 0)) & ((hr_2 == 0)|(hr_2_gr == 0)) & ((hr_3 == 0)|(hr_3_gr == 0)) & ((hr_4 == 1)|(hr_4_gr == 1)), "1", 0))) %>% 
  select(c(ID, date, hr1:hr4, multiple))


# finalized hpv dataset
tidy_hpv <- tidy_hpv %>% 
  left_join(HRgroup, by = c("ID", "date")) %>% 
  select(ID, date, results_hpv, hpvNum, hr1:hr4, multiple) %>%  # hr1~4 multiple
  left_join(hpvsingle_inf, by = c("ID", "date")) %>% 
  left_join(hpvmulti_inf, by = c("ID", "date")) %>%   # hr number multiple
  left_join(multi_16, by = c("ID", "date")) %>%
  left_join(multi_18, by = c("ID", "date")) %>%
  left_join(multi_1618, by = c("ID", "date"))

# -------------
# 3. biopsy
# 대표검체 종류
#### 제외확실: breast / skin / colon / stomach / lymph node / fetus / kidney / lung / peritoneum / urinary / fetus / anus
#### 제외: Placenta 태반 / vagina 질 / ovary 난소 / fallopian tube 자궁관, 난관 / vulva

c_bio <- cohort_t %>% 
  filter(처방코드 %in% cerv_code)

# replace roman numerals into arabic numbers
c_bio$cleanN <- gsub("\\bI\\b|\\bⅠ\\b", "1", c_bio$판독결과)
c_bio$cleanN <- gsub("\\bII\\b|\\bⅡ\\b", "2", c_bio$cleanN)
c_bio$cleanN <- gsub("\\bIII\\b|\\bⅢ\\b", "3", c_bio$cleanN)
c_bio$cleanN <- gsub("cervoccal", "cervical", c_bio$cleanN)

# trimming and cleaning
c_bio$clean_result <- c_bio$cleanN %>% 
  str_squish() %>% # remove heading & tailing white space
  str_to_lower() %>% 
  str_squish

c_bio$result = sub(".*diagnosis]|.*diagnosis>>|.*diagnosis >>", "diagnosis", c_bio$clean_result) %>% 
  str_replace_all("[[:punct:]]+", "") %>% 
  str_squish()

### re-trimming: delete unnecessary result
c_bio$result <- sub(" cf.*", "", c_bio$result)

### re-trimming clean_result for filtering
c_bio$clean_result <- c_bio$clean_result %>%
  str_replace_all("[[:punct:]]+", "") %>%
  str_squish() # remove white space again


############## 연구에 필요한 결과
# Uterine cervix 자궁경부 / uterus 자궁   *Endometrium은 자궁경부검사와 함께 진행되어 경부관련 결과는 같은 날짜에 수행된 자궁경부검사에서 확인할 수 있으므로 제외
c_bio <- c_bio[-c(2:5)] %>% 
  mutate(
    uterinec = case_when(
      str_detect(c_bio$clean_result, pattern = "대표검체 uterine cervix") ~ 1),
    uterus = case_when(
      str_detect(c_bio$clean_result, pattern = "대표검체 uterus") ~ 1)) %>% 
  filter(uterinec == 1 | uterus == 1) %>%  # cerv related results
  select(-c(uterinec, uterus))

# ---------------------------------------------------------------------------
#### possible outcome
# 1. CIN1 or 2 or 3
# 2. LSIL or HSIL
# 3. CIN3: Adenocarcinoma in situ(AIS) 
# 3. Cancer: invasive adenocarcinoma or Adenocarcinoma or Squamous cell carcinoma(SCC) or adenosquamous carcinoma
# 4. Negative or Free from tumor or No tumor present or No epithelial lesion
# 5. not detected

# CIN formatted
# Define a function to extract CIN values
extract_cin <- function(text) {
  # Define a regular expression pattern to match CIN-like patterns
  pattern <- "cin\\d|cin\\s\\d|cin\\d\\d|cin\\s\\d\\d|cervical\\sintraepithelial\\sneoplasia\\s\\d|cervical\\sintraepithelial\\sneoplasia\\d|cin\\sgrade\\s\\d"
  
  # Use the regmatches function to extract matching words
  matches <- regmatches(text, gregexpr(pattern, text))
  
  # Flatten the list of matches into a vector
  cin_words <- unlist(matches)
  
  # If there are matches, return the first match; otherwise, return NA
  if (length(cin_words) > 0) {
    return(cin_words[1])
  } else {
    return(NA)
    
  }
}

c_bio <- c_bio %>% 
  mutate(CIN = sapply(c_bio$result, extract_cin))

c_bio$CIN <- c_bio$CIN %>%
  str_replace_all(" ", "") %>%
  str_squish()

# table(c_bio$CIN)

c_bio$CIN <- gsub("cervicalintraepithelialneoplasia", "cin", c_bio$CIN)
c_bio$CIN <- gsub("grade", "", c_bio$CIN)
c_bio$CIN <- gsub("22", "2", c_bio$CIN)
c_bio$CIN <- gsub("23|32", "3", c_bio$CIN)

# cin1 cin2 cin3 
# 496  994 1038 

CIN <- c_bio %>% 
  filter(!(is.na(CIN)))


# grade-formatted
# Define a function to extract grade values
extract_grade <- function(text2) {
  # Define a regular expression pattern to match CIN-like patterns
  pattern2 <- "grade\\d|grade\\s\\d|hsil|high grade squamous|lsil|low grade squamous|lowgrade squamous|highgrade squamous|squamous dysplasia mild|condyloma|flat condylomaab|flat condyloma|condylomab|mild dysplasia focal|squamous dysplasia moderate|no remaining high"
  
  # Use the regmatches function to extract matching words
  matches2 <- regmatches(text2, gregexpr(pattern2, text2))
  
  # Flatten the list of matches into a vector
  grade_words <- unlist(matches2)
  
  # If there are matches, return the first match; otherwise, return NA
  if (length(grade_words) > 0) {
    return(grade_words[1])
  } else {
    return(NA)
  }
}


c_bio <- c_bio %>% 
  mutate(GRD = sapply(c_bio$result, extract_grade))

grade <- c_bio %>%
  filter(!(is.na(GRD))) %>% 
  mutate(CIN1_g = ifelse(str_detect(GRD, "grade1|grade 1|neoplasia 1|neoplasia1|lsil|low grade squamous|lowgrade squamous|squamous dysplasia mild|condyloma|flat condylomaab|flat condyloma|condylomab|mild dysplasia focal"), 1, 0),
         HSIL_g = ifelse(str_detect(GRD, "hsil|high grade squamous|highgrade squamous|squamous dysplasia moderate"), 1, 0),
         CIN2_g = ifelse(str_detect(GRD, "grade2|grade 2|neoplasia 2|neoplasia2"), 1, 0),
         CIN3_g = ifelse(str_detect(GRD, "grade3|grade 3|neoplasia 3|neoplasia3"), 1, 0),
         neg_g = ifelse(str_detect(GRD, "no remaining high"), 1, 0)) %>% 
  select(ID, date, CIN1_g:neg_g)

c_bio <- c_bio %>% 
  left_join(grade, by = c("ID", "date"))

# why <- c_bio %>%
#   filter(is.na(GRD)) %>% 
#   filter((CIN1_g==1)|(CIN2_g==1)|(CIN3_g==1)|(HSIL_g==1)|(neg_g==1))

# cancer
# 추후 adenocarcinoma in situ는 CIN3로, no tumor due to adenocarcinoma~~ 같은 형식도 있으므로 또 다시 정제작업 필요
extract_cancer <- function(text3) {
  # Define a regular expression pattern to match CIN-like patterns
  pattern3 <- "squamous cell carcinoma in situ|squamous carcinoma in situ|adenocarcinoma in situ|adenocarcinoma|invasive adenocarcinoma|squamous cell carcinoma|adenosquamous carcinoma|carcinoma in situ|adenosquamous carcioma|no tumor due to adenocarcinoma|no remaining adenocarcinoma in situ"
  
  # Use the regmatches function to extract matching words
  matches3 <- regmatches(text3, gregexpr(pattern3, text3))
  
  # Flatten the list of matches into a vector
  cancer_words <- unlist(matches3)
  
  # If there are matches, return the first match; otherwise, return NA
  if (length(cancer_words) > 0) {
    return(cancer_words[1])
  } else {
    return(NA)
  }
}

c_bio <- c_bio %>% 
  mutate(cancer = sapply(c_bio$result, extract_cancer))

cancer <- c_bio %>% 
  filter(!(is.na(cancer)))

# negative or free tumor
extract_neg <- function(text4) {
  # Define a regular expression pattern to match CIN-like patterns
  pattern4 <- "free from tumor|no epitherlial lesion|no epithelial|no tumor|negative|no tumor due to adenocarcinoma|no evidence of malignancy|no remaining adenocarcinoma in situ"
  
  # Use the regmatches function to extract matching words
  matches4 <- regmatches(text4, gregexpr(pattern4, text4))
  
  # Flatten the list of matches into a vector
  neg_words <- unlist(matches4)
  
  # If there are matches, return the first match; otherwise, return NA
  if (length(neg_words) > 0) {
    return(neg_words[1])
  } else {
    return(NA)
  }
}

c_bio <- c_bio %>% 
  mutate(neg = sapply(c_bio$result, extract_neg))

neg <- c_bio %>% 
  filter(!(is.na(neg)))

##########
# NA format(not detected)
extract_nd <- function(text5) {
  pattern5 <- "chronic cervicitis|chronic cervictis|chronic cervicits|endocervical polyp|chronic endocervicitis|mucoid materials insufficient|chronic inflammation|chronic inflammationa|uterine cervix polypectomy|favor reactive change|suspicious atrophy|favor atrophic change|fragments endocervical cells|inflammatory exudate|acute endocervicitis|atrophic cervicitis|cervical polyp|dexocervical resected margin involvement absent|dysplastic cells|cronic cervicitis|nonneoplastic endocervical tissue|fragments of atrophic"
  
  # Use the regmatches function to extract matching words
  matches5 <- regmatches(text5, gregexpr(pattern5, text5))
  
  # Flatten the list of matches into a vector
  nd_words <- unlist(matches5)
  
  # If there are matches, return the first match; otherwise, return NA
  if (length(nd_words) > 0) {
    return(nd_words[1])
  } else {
    return(NA)
  }
}


c_bio <- c_bio %>% 
  mutate(nd = sapply(c_bio$result, extract_nd))

nd <- c_bio %>% 
  filter(!(is.na(nd)))

# 제외
extract_out <- function(text6) {
  pattern6 <- "endometrial polyp|endometrium polypectomy|uterus myomectom|leiomyoma|paratubal cyst|proliferative endometrium|secretory endometrium|submucosal leiomyoma|endometrial hyperplasia focal atypiab|fallopian tube left cystectomy|fallopian tubes bilateral unremarkable|findings suggestive leiomyoma submucosal|foreign body reactio|lymph node dissection|ovaries bilateral unremarkabl|suggestive endometrial polyp|consistent leiomyom|endometrial hyperplasia without atypi|endometriosis consistent|fibrosis stroma|leiomyoma degeneration|mixed phase pattern composed secretory proliferative glands|secretory endometrium glandular stromal breakdow|unremarkable endometriumc|uterine myometrium myomectomy leiomyoma|uterus labeled myoma myomectomy"
  
  # Use the regmatches function to extract matching words
  matches6 <- regmatches(text6, gregexpr(pattern6, text6))
  
  # Flatten the list of matches into a vector
  out_words <- unlist(matches6)
  
  # If there are matches, return the first match; otherwise, return NA
  if (length(out_words) > 0) {
    return(out_words[1])
  } else {
    return(NA)
  }
}

c_bio <- c_bio %>% 
  mutate(out = sapply(c_bio$result, extract_out))

out <- c_bio %>% 
  filter(!(is.na(out))) ####### 확인 후 out할 자료만 뺄 것...다른 형식으로 진단명 나온 사람들 다수 포함되어 있음!!!

# finalized the results
# 1. ...in situ는 CIN3이기에 cancer column에서 제외
# 2. no tumor due to adenocarcinoma~~는 negative이므로...
# 3. see comment file 결과 반영
# 4. 정제한 결과 확인
cin3_pttn <- c("no remaining adenocarcinoma in situ", "no remaining adenocarcinoma in situ", "adenocarcinoma in situ", "squamous carcinoma in situ", "squamous cell carcinoma in situ", "carcinoma in situ")

cin3 <- c_bio[-c(4, 5)] %>%
  filter(cancer %in% cin3_pttn) %>% 
  mutate(cin3_p = 1) %>% 
  select(ID, date, cin3_p)

c_bio <- c_bio %>% 
  mutate(cancer = gsub("adenocarcinoma in situ|squamous carcinoma in situ|squamous cell carcinoma in situ|carcinoma in situ|no tumor due to adenocarcinoma|no remaining adenocarcinoma in situ", NA, cancer))

biopsy <- left_join(c_bio, cin3, by = c("ID", "date"))

##### out인 것인지, NA인지 혹은 다른 결과 추출할 수 있는 것인지 확인하기
# stillis.na <- biopsy[-c(4, 5)] %>%
#   filter(is.na(CIN)&is.na(CIN3)&is.na(GRD)&is.na(cancer)&is.na(neg)&is.na(nd)&is.na(out))

#### combine see comment file results added na0929 data as well
library(xlsx)
comments <- read.xlsx("C:/Users/hi/Documents/CINProgression/material/see_comments_0917.xlsx", sheetName = "forR")
combine <- left_join(biopsy, comments, by = c("ID", "date"))

tidy_biopsy <- combine %>% 
  mutate(CIN1 = (ifelse((CIN == "cin1")|(CIN1_g == 1), 1, 0)),
         CIN2 = (ifelse((CIN == "cin2")|(CIN2_g == 1), 1, 0)),
         CIN3 = (ifelse((CIN == "cin3")|(CIN3_g == 1)|(cin3_p == 1), 1, 0)),
         cancer = (ifelse((!is.na(cancer)), 1, 0)),
         negative = (ifelse((!is.na(neg))|(neg_g==1), 1, 0)),
         nd = (ifelse((!(is.na(nd))), 1, 0)),
         hsil = (ifelse(HSIL_g == 1, 1, 0))) %>% 
  mutate_at(vars(CIN1, CIN2, CIN3, negative, hsil, nd), as.numeric) %>% 
  mutate_at(vars(CIN1, CIN2, CIN3, negative, hsil, nd), ~replace_na(., 0)) %>% 
  mutate_at(vars(CIN1, CIN2, CIN3, negative, hsil, nd), as.factor) %>% 
  select(ID, date, CIN1, CIN2, CIN3, cancer, negative, hsil, nd, result_c)


#########################
biopsy_result_c <- tidy_biopsy %>% 
  filter(CIN1 == 0) %>% 
  filter(CIN2 == 0) %>% 
  filter(CIN3 == 0) %>% 
  filter(cancer == 0) %>% 
  filter(negative == 0) %>% 
  filter(nd == 0) %>% 
  mutate(CIN1_c = (ifelse((result_c == "CIN1")|(result_c == "LSIL"), 1, 0)),
         CIN2_c = (ifelse(result_c == "CIN2", 1, 0)),
         CIN3_c = (ifelse(result_c == "CIN3", 1, 0)),
         cancer_c = (ifelse(result_c == "cervical cancer", 1, 0)),
         negative_c = (ifelse(result_c == "negative", 1, 0))) %>% 
  select(ID, date, CIN1_c:negative_c)


# biopsy_result_c_no_info <- biopsy_result_c %>% filter(hsil == 0) %>% unique() # 146
# biopsy_result_c_hsil <- biopsy_result_c %>% filter(hsil == 1) %>% unique() # 58


tidy_biopsy <- left_join(tidy_biopsy, biopsy_result_c, by = c("ID", "date")) %>% 
  mutate(CIN1 = ifelse((CIN1 == 1)|(CIN1_c == 1), 1, 0),
         CIN2 = ifelse((CIN2 == 1)|(CIN2_c == 1), 1, 0),
         CIN3 = ifelse((CIN3 == 1)|(CIN3_c == 1), 1, 0),
         cancer = ifelse((cancer == 1)|(cancer_c == 1), 1, 0),
         negative = ifelse((negative == 1)|(negative_c == 1), 1, 0)) %>% 
  select(-c(result_c:negative_c))


# combine results by tests
df_list <- list(tidy_biopsy, tidy_pap, tidy_hpv)
comb_dt <- df_list %>%  reduce(full_join, by = c("date", "ID")) %>% 
  arrange(ID, date)

# Num of visits
dtCounts <- comb_dt %>% distinct(ID, date, .keep_all = T)

# add up by ID
dtCounts$numofvs <- ave(integer(nrow(dtCounts)), dtCounts$ID, FUN = seq_along)
dtCounts <- dtCounts %>%
  group_by(ID) %>% 
  mutate(count = n_distinct(date))

# 1st visit & last visit
dt_count <- dtCounts %>%
  group_by(ID) %>%
  dplyr::mutate(earliestD = min(date), latestD = max(date)) %>% 
  arrange(ID, date)

# time difference btw visits
dt_count <- dt_count %>%
  group_by(ID) %>%
  mutate(dftime = date - lag(date))

# time differ btw 1st and last visit
dt_count <- dt_count %>%
  group_by(ID) %>%
  dplyr::mutate(latestD - earliestD) %>% 
  arrange(ID, date)

# 2. 수술 기록 , 수술 처방
# 수술기록(모든 수술코드 hysterectomy)
hyst <- cohort_r %>%
  filter(!(is.na(수술코드))) %>% 
  mutate(hys = 1) %>% 
  select(ID, date, 수술일자, hys)

# 수술처방(leep)
leep <- cohort_s %>% 
  filter(수술처방코드 == "JR4262") %>% 
  mutate(leep = 1) %>% 
  select(ID, date, 수술처방일자, leep)


# 3. diagnosis data
d_cleaned <- cohort_d %>%
  filter(진단상태.확정.Rule.Out. == "C") %>% 
  mutate(cancer_d = ifelse(str_detect(진단코드, "C53"), 1, 0),
         hbv = ifelse(str_detect(진단명, "HBV"), 1, 0),
         hepatitis = ifelse(str_detect(진단명, "hepatitis|Hepatitis"), 1, 0),
         aml = ifelse(str_detect(진단코드, "C8|C9"), 1, 0), # 혈액암
         anemiaThrom = ifelse(str_detect(진단명, "anemia|hemoglobinuria|ITP|thromb|Thromb"), 1, 0),
         pid = ifelse(str_detect(진단코드, "N701.000|N709|N72|N730|N731|N738|N739|N76|N77"), 1, 0),  # pelvic inflammatory disease
         esrd = ifelse(str_detect(진단코드, "N185"), 1, 0),
         renalf = ifelse(str_detect(진단명, "Renal failure"), 1, 0)) %>% 
  select(-c(2, 3, 5:11)) %>% 
  filter(cancer_d == 1|hbv == 1|hepatitis == 1|aml == 1|anemiaThrom == 1|pid == 1|esrd == 1|renalf == 1)


# 4. prescription data
p_cleaned <- cohort_p %>%
  group_by(ID, date) %>%
  mutate(cervpunch = ifelse(str_detect(처방코드, "GC85"), 1, 0),
         gardasil = ifelse(str_detect(처방코드, "DV-HP|DV-9|5FA|VHDV|XDV"), 1, 0),
         cervarix = ifelse(str_detect(처방코드, "DV-JH|DV-JS|XD-HPVR"), 1, 0),
         std = ifelse(str_detect(처방코드, "LPD|PMO|VHLPD|VPH"), 1, 0),
         uro = ifelse(str_detect(처방코드, "LMF|LMR"), 1, 0)) %>% 
  select(-c(2:5))



###### make a final data set #######
dflist <- list(hyst, leep, d_cleaned, p_cleaned)

combine <- dflist %>%
  reduce(full_join, by = c("date", "ID")) %>% 
  arrange(ID, date)

tidydata <- cohort_c %>% 
  left_join(dt_count, by = "ID") ################ merge or join combine df

tidydata <- merge(tidydata, combine, by = c("ID", "date"), all = T)

write.csv(tidydata, "C:/Users/hi/Documents/CINProgression/data/cohort_temp0108.csv", row.names = F, fileEncoding = "euc-kr")
