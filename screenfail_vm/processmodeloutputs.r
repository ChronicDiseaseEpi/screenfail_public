# processmodeloutputs.R

## read in multivariate models
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
# library(brms)
library(tidybayes)
library(broom)
library(cmdstanr)

## set working directory if on linux
if(R.version$os == "linux-gnu") setwd("~/Documents/screenfail")
set_cmdstan_path("/opt/cmdstanr/cmdstan-2.29.2/")

## read in models ----
modname <- readRDS("model_with_data.Rds")
cmpr_cond <- map(paste0("mod", 1:9), ~ readRDS(paste0(.x, ".Rds")))
names(cmpr_cond) <- modname$model
cond <- map(paste0("mod_cond", 1:9), ~ readRDS(paste0(.x, ".Rds")))
names(cond) <- modname$model
smpl <- map(paste0("mod_smpl", 5), ~ readRDS(paste0(.x, ".Rds")))
names(smpl) <- modname$model[[5]]

## summarise posteriors and count number of trials per model ----
SummaryFx <- function(x) {
  a <- summarise_draws(x, ~quantile(.x, probs = c(0.025, 0.975)), .cores = 16)
  b <- summarise_draws(x, .cores = 16)  
  bind_cols(b, a[,-1])
}
cmpr_cond_smry <- map(cmpr_cond, SummaryFx) %>% 
  bind_rows(.id = "model")
cond_smry <- map(cond, SummaryFx) %>% 
  bind_rows(.id = "model")
smpl_smry <- map(smpl, SummaryFx) %>% 
  bind_rows(.id = "model")

## Limit the "smpl" and "cmpr_cond" to b and sd
cmpr_cond_smry <- cmpr_cond_smry %>% 
  filter(str_detect(variable, "^(b|sd)"))
smpl_smry <- smpl_smry  %>% 
  filter(str_detect(variable, "^(b|sd)"))
cond_smry <- cond_smry %>% 
  filter(str_detect(variable, "^(b|sd|r)"), 
         ! (str_count(variable, "_") == 1 & str_detect(variable, "^r")))

## for the b; the terms in the square brackets are the coefficients (eg intercept, age, sex, etc); can get from asce model and others
modname$terms <- map(modname$mods_cmpr, ~ .x$data$X %>% colnames())
TermsLevelsExtract <- function(x) {
  x %>% 
    mutate(term = case_when(
      str_detect(variable, "^(b|sd)") ~ str_extract(variable, "\\[[0-9]{1,2}\\]"),
      str_detect(variable, "^r") ~ str_extract(variable, "[0-9]{1,2}\\["),
      TRUE ~ ""
    ),
    term = str_remove_all(term, "\\[|\\]") %>% as.integer(),
    lvl = case_when(
      str_detect(variable, "^b") ~ "0",
      str_detect(variable, "^(sd|r)") ~ str_extract(variable, "_[0-9]{1,2}"),
      TRUE ~ ""
    ),
    lvl = str_remove_all(lvl, "_") %>% as.integer(),
    variable_orig = variable,
    variable = str_extract(variable, "[a-z]{1,2}"))
}
cond_smry <- TermsLevelsExtract(cond_smry)
smpl_smry <- TermsLevelsExtract(smpl_smry)
cmpr_cond_smry <- TermsLevelsExtract(cmpr_cond_smry)

## label levels
smpl_smry <- smpl_smry %>% 
  mutate(lvl_c = c("pop", "nct_id")[lvl+1])
cond_smry <- cond_smry %>% 
  mutate(lvl_c = c("pop", "cond", "nct_id")[lvl+1])
cmpr_cond_smry <- cmpr_cond_smry %>% 
  mutate(lvl_c = c("pop", "cmpr", "cond", "nct_id")[lvl+1])

## label terms
smry <- bind_rows(smpl = smpl_smry,
                  cond = cond_smry,
                  cmpr_cond = cmpr_cond_smry,
                  .id = "hier") %>% 
  group_by(model) %>% 
  nest() %>% 
  ungroup()
smry <- smry %>% 
  inner_join(modname %>% select(model, terms))
smry$data <- map2(smry$data, smry$terms, function(mydf, trm){
  mydf %>% 
    mutate(term_c = trm[term] %>% str_sub(5))
})
smry <- smry %>% 
  select(-terms) %>% 
  unnest(data)
trls <- modname %>% 
  select(model, cfs) %>% 
  unnest(cfs) %>% 
  distinct(model, nct_id) %>% 
  count(model) %>% 
  rename(trials = n)
smry <- smry %>% 
  left_join(trls)

## now need to pull the Rs as samples
asce <- cond$asce$draws(format = "df")
asce <- asce[, smry %>% 
               filter(model == "asce", hier == "cond") %>% 
               pull(variable_orig)] %>% 
  as_tibble() %>% 
  gather(key = "variable_orig", "value")
asce <- asce %>% 
  inner_join(smry %>% 
               filter(model == "asce", hier == "cond", lvl_c %in% c("pop", "cond"), variable %in% c("b", "r")) %>% 
               select(variable_orig, variable, lvl, lvl_c, term, term_c))
asce <- asce %>% 
  mutate(iter = rep(1:4000, nrow(asce)/4000))
asce <- asce %>% 
  mutate(cond_n = case_when(
    variable == "b" ~ 0L,
    variable == "r" ~ str_extract(variable_orig, "\\[[0-9]{1,2}") %>% str_remove("\\[") %>% as.integer(),
    TRUE ~ NA_integer_)
  )
asce <- asce %>% 
  # filter(iter == 1, term_c == "Intercept")  %>% 
  group_by(iter, term_c) %>% 
  mutate(b = value[1]) %>% 
  filter(!variable == "b") %>% 
  ungroup() %>% 
  mutate(value  = b+value) %>% 
  select(iter, term_c, cond_n, value)
asce <- asce %>% 
  group_by(term_c, cond_n) %>% 
  summarise(est = mean(value),
            q2_5 = quantile(value, probs = c(0.025)),
            q97_5 = quantile(value, probs = c(0.975))) %>% 
  ungroup()
conds <- modname$cfs[modname$model == "asce"][[1]] %>% 
  pull(cond) %>% 
  unique() %>% 
  sort()
asce <- asce %>% 
  mutate(cond_c = conds[cond_n])
write_csv(asce, "Outputs/asce_by_cond.csv")

## now pulled R can drop them from smry
smry <- smry %>% 
  filter(!variable == "r")
write_csv(smry, "Outputs/overall_outputs.csv")
## diagnostics; no divergent transitions for main models
MakeDiag <- function(models){
  my_diag <- map(models, ~ .x$diagnostic_summary())
  my_diag <- map(my_diag, ~ map(.x, function(y) {
    names(y) <- paste0("chain", 1:4)
    y
  }))
  my_diag <- map(my_diag, as_tibble)
  bind_rows(my_diag, .id = "model")
}
cond_diag <- MakeDiag(cond)
cmpr_cond_diag <- MakeDiag(cmpr_cond)
smpl_diag <- MakeDiag(smpl)

write_csv(cond_diag, "Outputs/cond_diag.csv")
write_csv(cmpr_cond_diag, "Outputs/cmpr_cond_diag.csv")
write_csv(smpl_diag, "Outputs/smpl_diag.csv")

## Pull interaction estimate for men and women
asce_int <- readRDS("mod_cond8.Rds")
asce_int <- asce_int$draws(format = "df")
asce_int <- asce_int[, smry %>% 
               filter(model == "asce_int", hier == "cond", term_c %in% c("black", "male:black")) %>% 
               pull(variable_orig)] %>% 
  as_tibble() %>% 
  gather(key = "variable_orig", "value")
asce_int <- asce_int %>% 
  inner_join(smry %>% 
               filter(model == "asce_int", hier == "cond", lvl_c %in% c("pop"), variable %in% c("b"), term_c %in% c("black", "male:black")) %>% 
               select(variable_orig, variable, lvl, lvl_c, term, term_c))
asce_int <- asce_int %>% 
  mutate(iter = rep(1:4000, nrow(asce_int)/4000))
asce_int <- asce_int %>% 
  mutate(term_c = case_when(
    term_c == "male:black" ~ "nter",
    term_c == "black" ~ "black_women",
    TRUE ~ NA_character_))
asce_int <- asce_int %>% 
  select(term_c, iter, value) %>% 
  spread(term_c, value)
asce_int <- asce_int %>%
  select(-iter) %>% 
  mutate(black_men = black_women + nter)

asce_int_smry <- SummaryFx(asce_int)
asce_int_smry %>% 
  transmute(across(c(mean, median, `2.5%`, `97.5%`), exp))
