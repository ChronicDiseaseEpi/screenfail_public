# 08_new_plots_tables
library(tidyverse)
mydf <- read_csv("screenfail_vm/Outputs/overall_outputs.csv")

## Produce model summary for ASCE model for different levels fo nesting ----
asce <- mydf %>% 
  filter(model == "asce") %>% 
  rename(est = mean, lci = `2.5%`, uci = `97.5%`) %>% 
  mutate(across(c(est, lci, uci), function(x) x %>% 
                  round(2) %>% 
                  formatC(format = "f", digits = 2) ),
         res = paste0(est, " (", lci, " to ", uci, ")") )
tbl_mod_hier <- asce %>% 
  mutate(hier = case_when(
    hier == "smpl" ~ "Trial",
    hier == "cond" ~ "Trial and condition",
    hier == "cmpr_cond" ~ "Trial, condition and treatment"
  ),
  lvl_c = case_when(
    lvl_c == "nct_id" ~ "Between   Trial SD",
    lvl_c == "cond" ~ "Between  condition SD",
    lvl_c == "cmpr" ~ "Between comparison SD"
  )) 

lkp <- c("Intercept",
         "Age (decades)",
         "Male", 
         "Comorbidity count", 
         "Asian", 
         "Black or African American", 
         "American Indian or Alaska Native or Native Hawaiian or Other Pacific Islander", 
         "Other")
names(lkp) <- c("Intercept",
                "age10", 
                  "male", 
                  "como_ct", 
                  "asian", 
                  "black", 
                  "nativ", 
                  "other")
sds <- tbl_mod_hier %>% 
  filter(variable == "sd") %>% 
  select(hier, term_c, res, lvl_c) %>% 
  mutate(term_c = factor(term_c, levels = names(lkp), labels = (lkp))) %>% 
  spread(lvl_c, res) %>% 
  arrange(term_c) %>% 
  select(term_c, everything())
sds
write_csv(sds, "Outputs/Table3.csv")
