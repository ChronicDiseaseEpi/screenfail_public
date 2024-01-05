library(tidyverse)
library(brms)
library(tidybayes)
library(broom)

## Creates Stan scripts and data structures for running in virtual machine

## Functions ----
Brm <- function(...){
  ## Instead of running Brm model, create code and data for running in VM
  a <- brms::make_stancode(...)
  b <- brms::make_standata(...)
  list(code = a, data = b)
}

## set working directory if on linux
if (R.version$os == "linux-gnu") setwd("~/Documents/screenfail")

## if on linux set moment_match to FALSE
my_moment_match <- !R.version$os == "linux-gnu"

## Read in data ----

## Read in trial metadata
metadata1 <- read_csv("Data/hte_metadata.csv")
metadata2 <- read_csv("Data/new_metadata.csv") %>% 
  select(nct_id, cmpr, cond = condition) %>% 
  mutate(cond = case_when(
    cond %in% c("Diabetes Mellitus, Type 2, nephropathy",
                "Diabetic Nephropathies",
                "Hypertension; Diabetic Nephropathies") ~ "Diabetes",
    cond == "Heart failure, Congestive and Microalbuminuria" ~ "Heart failure",
    cond %in% c("Rhinitis, Allergic, Perennial and Seasonal", "Rhinitis, Allergic, Seasonal", 
                "Rhinitis, Vasomotor", "Rhinitis, Allergic, Perennial") ~ "Rhinitis",
    TRUE ~ cond
  ))
metadata <- bind_rows(metadata1, metadata2)
rm(metadata1, metadata2)

## identify which trials comorbidity count is unreliable
gsk <- read_csv("Data/full_baseline_table.csv")
viv <- read_csv("Data/baseline_table_trial_status.csv")
tot <- bind_rows(gsk = gsk, viv = viv, .id = "repo") %>% 
  mutate(status = if_else(status == "fail", "sf", status))
unreliable <- tot %>% 
  filter(como_ct_total == 0) %>% 
  distinct(nct_id) %>% 
  pull(nct_id)
## identify further unreliable ones with low comorbidity count (very low, suggests not recorded in screen failure)
unreliable2 <- read_csv("nct_id
NCT01159912
NCT00118703
NCT00051558
NCT00117325
NCT01181895
NCT00968708
NCT00827242")
unreliable <- union(unreliable, unreliable2$nct_id)
reliable <- setdiff(tot$nct_id, unreliable)

# pull in results from vivli trials (log reg models)
res_viv <- read_csv("Data/modelcofs.csv") # %>% filter(!is.na(p.value))
res_sa_viv <- read_csv("Data/modelcofs_sa.csv") # %>% filter(!is.na(p.value))
res_sa2_viv <- read_csv("Data/modelcofs_sa2.csv") # %>% filter(!is.na(p.value))
nct_id_viv <- res_viv$nct_id %>% unique() %>% sort()


# pull in results from gsk trials (log reg models)
res_gsk <- read_csv("Data/modelcofs_gsk.csv") # %>% filter(!is.na(p.value))
res_sa_gsk <- read_csv("Data/modelcofs_sa_gsk.csv") # %>% filter(!is.na(p.value))
res_sa2_gsk <- read_csv("Data/modelcofs_sa2_gsk.csv") # %>% filter(!is.na(p.value))
nct_id_gsk <- res_gsk$nct_id %>% unique() %>% sort()

# pull in results from variance-covariance matrices for vivli and gsk
res_vcov_viv <- read_csv("Data/modelvarcovs.csv")
res_vcov_gsk <- read_csv("Data/modelvarcovs_gsk.csv")

## Transform and harmonise data from different sources ----
# bind results together
res <- bind_rows(viv = res_viv, gsk = res_gsk, .id = "repo")
res_sa <- bind_rows(res_sa_viv, res_sa_gsk) # sensitivity analysis: only those with >= 1 como 
res_sa2 <- bind_rows(res_sa2_viv, res_sa2_gsk) # sensitivity analysis: only those with >= 2 como 
res_vcov <- bind_rows(res_vcov_viv, res_vcov_gsk)
rm(res_viv, res_gsk, res_sa_viv, res_sa2_gsk, res_vcov_viv, res_vcov_gsk)

## examine completeness; add this as a variables
cmplt <- res %>% 
  filter(!is.na(estimate)) %>%
  mutate(term = if_else(term %in% c("asian", "black", "other", "nativ"), "eth", term)) %>% 
  distinct(nct_id, term) %>%
  count(nct_id, term) %>% 
  filter(term %in% c("age10", "male", "como_ct", "eth")) %>% 
  spread(term, n, fill = 0L) %>% 
  select(nct_id, age10, male, everything())
## set single sex trials to "male ==1" as male is NOT missing for these - note two vivli trials
## still both missing age:- NCT02183064 and NCT01335464
cmplt <- cmplt %>% 
  mutate(male = if_else(nct_id %in% c("NCT02183064", "NCT00670501", "NCT00861757", "NCT00384930", 
                       "NCT00827242", "NCT00855582", "NCT00848081", "NCT00970632", "NCT01335464", 
                       "NCT00670319"), 1L, male))


map(cmplt %>% select(-nct_id), table)
## drop any missing models
res <- res %>% 
  inner_join(cmplt) %>% 
  filter(
    !(model %in% c("as", "as_int")     & (age10 == 0 | male == 0)),
    !(model %in% c("asc", "asc_int")   & (age10 == 0 | male == 0 | como_ct == 0)),
    !(model %in% c("asce", "asce_int") & (age10 == 0 | male == 0 | como_ct == 0 | eth == 0)),
    !(model %in% c("cc") & (como_ct == 0)),
    !(model %in% c("eth") & (eth == 0)))
cmplt <- cmplt %>% 
  mutate(tot = paste0(age10, male, como_ct, eth),
         tot = case_when(
           tot == "1111" ~ "sexcometh",
           tot == "1001" ~ "eth",
           tot == "1010" ~ "com",
           tot == "1011" ~ "cometh",
           tot == "1101" ~ "sexeth",
           tot == "1110" ~ "sexcom"
         ))
permdl <- res %>% 
  filter(!is.na(estimate), !term == "(Intercept)") %>% 
  group_by(model) %>% 
  summarise(terms = paste(term %>% unique() %>% sort(), collapse = ", "),
            n = sum(!duplicated(nct_id))) %>% 
  ungroup()

## drop terms where have particular ones not featuring in trial for race/ethnicity
res <- res %>% 
  filter(!is.na(estimate)) 

## drop trials where comorbidity is unreliable for the comorbidity analyses
res <- res %>% 
  filter(model %in% c("as", "as_int", "eth") |
         nct_id %in% reliable)

## limit it to trials in tot, 5 dropped; two with missing age sex data,
## three others not in tot "NCT00384930" "NCT01335477" "NCT01769378". 1 BI and two Eli Lilly
## leaves 52 trials
setdiff(res$nct_id, tot$nct_id)
res <- res %>% 
  semi_join(tot %>% select(nct_id))

## Create data in block diagonal structure to allow fitting of MVN models ----
## limit to trials in res
res_vcov %>% 
  anti_join(res %>% select(model, nct_id)) %>%
  distinct(model, nct_id) %>% 
  count(model)
res_vcov <- res_vcov %>% 
  semi_join(res %>% select(model, nct_id, rows = term)) %>% 
  semi_join(res %>% select(model, nct_id, cols = term))

## Expand data so is symmetrical
res_vcov <- bind_rows(res_vcov,
                       res_vcov %>% rename(cols = rows, rows = cols)) 
res_vcov <- bind_rows(res_vcov,
                      res_vcov %>% 
                         distinct(model, nct_id, rows) %>% 
                         mutate(cols = rows,
                                values = 1)) %>% 
  arrange(nct_id, model)
## add in standard errors
res_vcov <- res_vcov %>% 
  inner_join(res %>% 
               select(model, nct_id, rows = term, rows_se = std.error) ) %>% 
  inner_join(res %>% 
               select(model, nct_id, cols = term, cols_se = std.error) )
## calculate variance-covariance from correlation and standard errors
res_vcov <- res_vcov %>% 
  mutate(cvrnc = values*rows_se*cols_se) %>% 
  select(-rows_se, -cols_se, -values)
## identify models where parameters not estimated;

## join data for fitting block diagonal models
res_nst <- res %>% 
  select(model, nct_id, term, estimate, std.error) %>% 
  nest(est = c(term, estimate, std.error))
vcv_nst <- res_vcov %>% 
  select(model, nct_id, rows, cols, cvrnc ) %>% 
  nest(vcv = c(rows, cols, cvrnc))
rm(res_vcov)
## drop ones without variance covariance; 
res_nocov <- res_nst %>% 
  anti_join(vcv_nst)
res_nst <- res_nst %>% 
  inner_join(vcv_nst)
rm(vcv_nst)

## add in comparison and condition
res_nst <- res_nst %>% 
  inner_join(metadata)
res <- res %>% 
  inner_join(metadata)
res_sa <- res_sa %>% 
  inner_join(metadata)
res_sa2 <- res_sa2 %>% 
  inner_join(metadata)

## simplify condition and treatment comparison names
res <- res %>% 
  mutate(cond = if_else(cond == "Rhinitis, allergic", "Rhinitis", cond))
res_sa <- res_sa %>% 
  mutate(cond = if_else(cond == "Rhinitis, allergic", "Rhinitis", cond))
res_sa2 <- res_sa2 %>% 
  mutate(cond = if_else(cond == "Rhinitis, allergic", "Rhinitis", cond))
res_nst <- res_nst %>% 
  mutate(cond = if_else(cond == "Rhinitis, allergic", "Rhinitis", cond))

## Select which multivariate models to fit ----
## all should be TRUE as tests of model matching etc
res_nst$test <- map2_lgl(res_nst$est, res_nst$vcv, ~ all(.x$term %in% .y$rows) & all(.y$rows %in% .x$term))
if (!all(res_nst$test)) warning("Mismatch")
res_nst$vcv <- map(res_nst$vcv, ~ {
  new <-  .x %>% spread(cols, cvrnc) %>% 
    select(-rows)
  new <- as.matrix(new)
  rownames(new) <- colnames(new)
  new
})

## Rounding errors have caused asymmetry so recreate as symmetrical matrix
## errors here some are not symmetrical review.
res_nst$symm <- map_lgl(res_nst$vcv, matrixcalc::is.symmetric.matrix)
res_nst %>% 
  count(symm)
res_nst$vcv <- map(res_nst$vcv, function(x) {
  dnames <- rownames(x)
  upper <- x[upper.tri(x)]
  dg <- diag(x)
  y <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  y[upper.tri(y)] <- upper
  y <- t(y)
  y[upper.tri(y)] <- upper
  diag(y) <- dg
  if (any(abs(x - y) > 1E-10)) warning("larger difference than expected")
  rownames(y) <- dnames
  colnames(y) <- dnames
  y
})
res_nst$symm <- map_lgl(res_nst$vcv, matrixcalc::is.symmetric.matrix)
res_nst %>% 
  count(symm)
res_nst$pd <- map_lgl(res_nst$vcv, matrixcalc::is.positive.definite)
res_nst %>% 
  count(pd)
res_nst %>% 
  filter(!pd)

## there are two trials with missing models for age and sex but 
## have como count and ethnicity "NCT02183064" and"NCT01335464"
## should be dropped from all as very small number of screen failures
## already NOT in tot
res_nst <- res_nst %>% 
  filter(!nct_id %in% c("NCT02183064", "NCT01335464"))
res <- res %>% 
  filter(!nct_id %in% c("NCT02183064", "NCT01335464"))
res_nst %>% 
  count(pd)
saveRDS(res, "Scratch_data/res.Rds")

## rearrange matrix so matches terms
res_nst$vcv <- map2(res_nst$est, res_nst$vcv, ~ .y[.x$term, .x$term])
res_nst$test2 <- map2_lgl(res_nst$est, res_nst$vcv, ~ all(.x$term == rownames(.y))) 
if (!all(res_nst$test2)) warning("Mismatch")

## print random sample to check
PrintRand <- function(){
  print(list(i <- sample(seq_along(res_nst$model), 1),
             res_nst$nct_id[i],
             res_nst$model[i],
  res_nst$est[[i]],
  res_nst$vcv[[i]]) )
}
PrintRand()

## create vector of y's and block diagonal
res_nst <- res_nst %>% 
  select(model, nct_id, cfs = est, cvs = vcv, cond, cmpr) %>%
  nest(data = c(nct_id, cond, cmpr, cfs, cvs))
res_nst2 <- res_nst %>% 
  filter(model == "eth") %>% 
  mutate(model = "eth_33")
res_nst2$data[[1]] <- res_nst2$data[[1]] %>% 
  filter(nct_id %in% res_nst$data[res_nst$model == "asce"][[1]]$nct_id)
res_nst <- bind_rows(res_nst,
                     res_nst2)
# arrange into block diagonal matrices
CreateDataLong <- function(a) {
  a <- a %>% 
    arrange(nct_id)
  cvs <- a %>% 
    select(nct_id, cvs)
  bd <- Matrix::bdiag(cvs$cvs)
  b <- a %>% 
    select(-cvs) 
  list(cfs = b %>% unnest(cfs),
       bd = bd)
}
res_nst$data <- map(res_nst$data, CreateDataLong)
res_nst$cfs <- map(res_nst$data, ~ .x[[1]])
res_nst$bd <- map(res_nst$data, ~ .x[[2]])
res_nst$data <- NULL

## Fit MVN models, time consuming, run on VM ----
## 
myform <- estimate ~ 0 + term  + (0 + term | nct_id) + (0 + term | cmpr) + (0 + term | cond) + fcor(covs) 
myform2 <- estimate ~ 0 + term  + (0 + term | nct_id) +                    (0 + term | cond) + fcor(covs) 
myform3 <- estimate ~ 0 + term  + (0 + term | nct_id)                                     + fcor(covs) 

## default prior for intercept on logistic regression model in rstanarm is 0, 1.4
## change to 0 and 2 as this has better coverage of low values and high values which are plausible 
## for screen failure
plot(density(plogis(rnorm(10000, 0, 2))), col = "red", xlim = c(0, 1))
lines(density(plogis(rnorm(10000, 0, 1.4))), col = "blue", xlim = c(0, 1))

mypriors <- c(
  prior(constant(1), class = "sigma"),
  set_prior("normal(0, 2)", class = "b", coef = "termIntercept"),
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(0, 1)", class = "sd"))
mypriors2 <- c(
  prior(constant(1), class = "sigma"),
  set_prior("normal(0, 4)", class = "b", coef = "termIntercept"),
  set_prior("normal(0, 2)", class = "b"),
  set_prior("normal(0, 2)", class = "sd"))

res_nst$mods_cmpr_cond <- pmap(list(res_nst$cfs, res_nst$bd, res_nst$model), function(x, y, z) {
  print(z)
  mod_res <- Brm(
    formula = myform,
    data = x, 
    data2 = list(covs = y),
    prior = mypriors,
    control = list(adapt_delta = 0.99),
    cores = 4)
})

## Save to run on new VM where can use brms directly ----
saveRDS(list(res_nst = res_nst %>% filter(model == "asce"), 
             myform = myform, 
             myform2 = myform2, 
             myform3 = myform3,
             mypriors = mypriors,
             mypriors2 = mypriors2), "Scratch_data/newVM.Rds")

## run a single model to obtain model objects
# onemod <- brm(myform, data = res_nst$cfs[[1]], data2 = list(covs = res_nst$bd[[1]]), prior = mypriors, chains = 0)
# saveRDS(onemod, "Scratch_data/one_model_compile.Rds")
onemod <- readRDS("Scratch_data/one_model_compile.Rds")
prior_summary(onemod)

## save as models for 
res_nst$mods_cond <- pmap(list(res_nst$cfs, res_nst$bd, res_nst$model), function(x, y, z) {
  print(z)
  mod_res <- Brm(
    formula = myform2,
    data = x, 
    data2 = list(covs = y),
    prior = mypriors,
    control = list(adapt_delta = 0.99),
    cores = 4)
})

res_nst$mods_smpl <- pmap(list(res_nst$cfs, res_nst$bd, res_nst$model), function(x, y, z) {
  print(z)
  mod_res <- Brm(
    formula = myform3,
    data = x, 
    data2 = list(covs = y),
    prior = mypriors,
    control = list(adapt_delta = 0.99),
    cores = 4)
})

## univariate normal likelihood models - create filtered datasets which can be used for Brm ----
res_sa <- res_sa %>% 
  filter(!is.na(p.value))
res_sa2 <- res_sa2 %>% 
  filter(!is.na(p.value))

dt_lst <- list(
  age_agesex = res %>% filter(term == "age10" & model == "as"), # age and sex models
  sex_agesex = res %>% filter(term == "male" & model == "as"), # age and sex models
  como_agesexcomo = res %>% filter(term == "como_ct" & model == "asc"), # age sex como models
  como01_agesexcomo = res_sa %>% filter(term == "como_ct" & model == "asc"),# age sex como models
  como2p_agesexcomo = res_sa2 %>% filter(term == "como_ct" & model == "asc"), # age sex como models
  age_agesexinter = res %>% filter(term == "age10" & model == "as_int"), # age sex interaction models
  sex_agesexinter = res %>% filter(term == "male" & model == "as_int"), # age sex interaction models
  black_agesexcomoeth = res %>% filter(term == "black" & model == "asce"), # age sex como ethnicity models
  asian_agesexcomoeth = res %>% filter(term == "asian" & model == "asce"), # age sex como ethnicity models
  nativ_agesexcomoeth = res %>% filter(term == "nativ" & model == "asce"), # age sex como ethnicity models
  other_agesexcomoeth = res %>% filter(term == "other" & model == "asce") # age sex como ethnicity models
)

## save multivariate for running on the VM
saveRDS(res_nst, "Scratch_data/model_with_data.Rds")
