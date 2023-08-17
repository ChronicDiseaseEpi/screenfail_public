## count participants
library(tidyverse)

## function ----
## Can calculate combined standard deviation of groups as per
## https://handbook-5-1.cochrane.org/chapter_7/table_7_7_a_formula_for_combining_groups.htm
CombMean <- function(n1, m1, n2, m2){
  top <- n1*m1 + n2*m2
  top/(n1 + n2)
}

CombSd <- function(n1, m1, n2, m2, s1, s2){
  top = (n1 - 1)*s1^2 + (n2 - 1)*s2^2 + (n1*n2/(n1 + n2)) * (m1^2 + m2^2 - 2*m1*m2)
  bottom = n1 + n2 - 1
  res = (top/bottom)^0.5
  res
}

CombSdVectorised <- function(n, m, s){
  ## get rows
  myrows <- length(n)
  myrows_now <- myrows
  while (myrows_now >= 2) {
    ## select first two values
    n1 <- n[1]
    n2 <- n[2]
    s1 <- s[1]
    s2 <- s[2]
    m1 <- m[1]
    m2 <- m[2]
    ## replace first value with combination of first two values and drop second value
    new_s1 <- CombSd(n1, m1, n2, m2, s1, s2)
    new_m1 <- CombMean(n1, m1, n2, m2)
    new_n1 <- n1 + n2
    s[1] <- new_s1
    m[1] <- new_m1
    n[1] <- new_n1
    s <- s[-2]
    m <- m[-2]
    n <- n[-2]
    # print(s)
    # print(m)
    ## recalculate the length
    myrows_now <- length(s)
  }
  # print(sum(n))
  s
}

## Read data ----
gsk <- read_csv("Data/full_baseline_table.csv")
viv <- read_csv("Data/baseline_table_trial_status.csv")
mdls <- readRDS("Scratch_data/model_with_data.Rds")

## process data ----
## number of trials each model
mdls$n_trials <- map_int(mdls$cfs, ~ sum(!duplicated(.x$nct_id)))
tot <- bind_rows(gsk = gsk, viv = viv, .id = "repo") %>% 
  mutate(status = if_else(status == "fail", "sf", status))

## drop como_count data where is incorrect
crct <- mdls$cfs[mdls$model == "cc"][[1]] %>% 
  distinct(nct_id) %>% 
  pull(nct_id)

## drop ethnicity data where is incorrect
crct_e <- mdls$cfs[mdls$model == "eth"][[1]] %>% 
  distinct(nct_id) %>% 
  pull(nct_id)
tot <- tot %>% 
  mutate(across(contains("como"), ~ if_else(nct_id %in% crct, .x, NA_real_)))
tot <- tot %>% 
  mutate(across(c(white, asian, black, ethnic_na, other, native ), ~ if_else(nct_id %in% crct_e, .x, NA_real_)),
         native = if_else(!is.na(white) & is.na(native), 0, native ))

## not sure whey these are not in it
tot <- tot %>% 
  mutate(across(c(age_mean, age_sd), ~ round(.x, 1))) %>% 
  select(-ethnic_na)

## add index condition
indx <- mdls %>% 
  select(cfs) %>% 
  unnest(cfs) %>% 
  distinct(cond, nct_id)
tot <- tot %>% 
  inner_join(indx)
## order is important here summarise SD before mean as it depends on means
## also como_ct_mean has to come later
## for processing for table impute age mean and sd for pulmonary hypertension
## note remove this later
tot <- tot %>% 
  group_by(nct_id) %>% 
  mutate(across(c(age_mean, age_sd), ~ if_else(nct_id == "NCT00125918", mean(.x, na.rm = TRUE), .x))) %>% 
  ungroup()
tot_smry <- tot %>% 
  group_by(cond, status) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            trials_como = sum(!duplicated(nct_id[!is.na(como_ct_mean)])),
            trials_eth = sum(!duplicated(nct_id[!is.na(white)])),
            age_sd = CombSdVectorised(n, age_mean, age_sd),
            age_mean = weighted.mean(age_mean, n, na.rm = TRUE),
            n_como = sum(n * !is.na(como_ct_mean)),
            n_eth = sum(n * !is.na(white)),
            across(c(n, male, como_ct_0:native), sum, na.rm = TRUE)) %>% 
  ungroup()
tot_smry2 <- tot %>% 
  mutate(cond = "All") %>% 
  group_by(cond, status) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            trials_como = sum(!duplicated(nct_id[!is.na(como_ct_mean)])),
            trials_eth = sum(!duplicated(nct_id[!is.na(white)])),
            age_sd = CombSdVectorised(n, age_mean, age_sd),
            age_mean = weighted.mean(age_mean, n, na.rm = TRUE),
            n_como = sum(n * !is.na(como_ct_mean)),
            n_eth = sum(n * !is.na(white)),
            across(c(n, male, como_ct_0:native), sum, na.rm = TRUE)) %>% 
  ungroup()
tot_smry <- bind_rows(tot_smry,
                      tot_smry2)
rm(tot_smry2)

MakeMean <- function(x, n){
  a <- round(100*x/n, 1)
  if_else(n == 0, "", paste0(x, " (", a, "%)"))
}
MakeBrac <- function(x, y) {
  if_else(x == 0, "0",  paste0(x, " (", y, ")"))
}
tot_smry <- tot_smry %>% 
  mutate(across(c(como_ct_0:como_ct_2plus), ~ MakeMean(.x, n_como)),
         across(c(white:native), ~ MakeMean(.x, n_eth)),
         age = paste0(round(age_mean), " (", round(age_sd), ")"),
         male = MakeMean(male, n),
         participants = n,
         trials_como = MakeBrac(trials_como, n_como),
         trials_eth = MakeBrac(trials_eth, n_eth)) %>% 
  select(cond, status, trials, participants, male, age, trials_como, como_ct_0:como_ct_2plus,trials_eth, white:native) %>% 
  mutate(age = if_else(cond == "Hypertension, Pulmonary" & status == "sf", "", age))
write_csv(tot_smry, "Outputs/Table1_new.csv")

tot <- tot %>% 
  mutate(age_mean = paste0(age_mean, " (", age_sd, ")"),
         como_ct_mean = paste0(como_ct_mean, " (", como_ct_sd, ")")) %>% 
  rename(age = age_mean,
         como_ct = como_ct_mean) %>% 
  select(-age_sd, -como_ct_sd) %>% 
  mutate(across(c(male, asian, black, white, other, native), ~ paste0(.x,
                      " (",
                      round(100*.x/n, 1),
         "%)")))
tot <- tot %>% 
  arrange(nct_id, status)

tot_out <- tot %>% 
  mutate_all(as.character()) %>% 
  mutate_all(function(x) if_else(x == "NA (NA%)", "", as.character(x)))
tot_out <- tot_out %>% 
  select(-trial, -repo)
tot_out <- tot_out %>% 
  mutate(status = if_else(status == "sf", "No", "Yes")) %>% 
  rename(Randomised = status)
names(tot_out)[-1] <- str_to_sentence(names(tot_out)[-1])
names(tot_out)[1] <- "NCT ID"
names(tot_out) <- str_replace(names(tot_out), "Como_ct", "Comorbidity count")
names(tot_out) <- str_replace(names(tot_out), "_", " ")
names(tot_out) <- str_replace(names(tot_out), "2plus", ">= 2")
tot_out <- tot_out[c("NCT ID", "Randomised", "N", "Age", "Male", "Comorbidity count", 
                     "Comorbidity count total", "Comorbidity count 0", "Comorbidity count 1", 
                     "Comorbidity count >= 2", "White", "Asian", "Black",
                     "Native",  "Other")]


write_csv(tot_out, "Outputs/table1_long_version_for_sa.csv", na = "")


