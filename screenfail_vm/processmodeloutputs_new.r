# processmodeloutputs.R

## read in multivariate models
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(brms)
library(tidybayes)
library(broom)

## read in models ----
model <- readRDS("FromVM/model_object.Rds")[[1]]

## summarise posteriors and count number of trials per model ----
cmpr_cond_smry <- summary(model)
cmpr_cond_smry <- cmpr_cond_smry$fixed
cmpr_cond_smry <- as_tibble(cmpr_cond_smry, rownames = "term", .name_repair = "universal")
names(cmpr_cond_smry) <- names(cmpr_cond_smry) %>% 
  str_to_lower() %>% 
  str_replace_all("\\.+", "_")
cmpr_cond_smry <- cmpr_cond_smry %>% 
  mutate(across(c(estimate, l_95_ci, u_95_ci), ~ round(exp(.x), 2)))
cmpr_cond_smry <- cmpr_cond_smry %>% 
  filter(!term == "termIntercept") %>% 
  mutate(term = str_remove(term, "^term"),
         res = paste0(estimate, " (", l_95_ci, " to ", u_95_ci, ")"))
cmpr_cond_smry <- cmpr_cond_smry %>% 
  mutate(term_f = factor(term,
                         levels = c("age10",
                                    "male",
                                    "como_ct",
                                    "asian",
                                    "black",
                                    "nativ",
                                    "other"),
                         labels = c("Age (decades)",
                                    "Male",
                                    "Comorbidity count",
                                    "Asian",
                                    "Black",
                                    "Indigenous",
                                    "Other"))) %>% 
  arrange(term_f) %>% 
  select(term_f, res)
write_csv(cmpr_cond_smry, "Outputs/wider_priors_sa.csv")
