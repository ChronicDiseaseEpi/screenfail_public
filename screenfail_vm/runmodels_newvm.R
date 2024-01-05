## code for running models on new VM where can run using brms directly and do not need to create stancode and standata locally
library(dplyr)
library(brms)
library(purrr)

res <- readRDS("newVM.Rds")
list2env(res, envir = .GlobalEnv)
rm(res)
## save as models for 
res_nst$mods_cond <- pmap(list(res_nst$cfs, res_nst$bd, res_nst$model), function(x, y, z) {
  print(z)
  mod_res <- brm(
    formula = myform2,
    data = x, 
    data2 = list(covs = y),
    prior = mypriors2,
    control = list(adapt_delta = 0.99),
    cores = 4)
})

x <- res_nst$mods_cond
SummaryFx <- function(x) {
  a <- summarise_draws(x, ~quantile(.x, probs = c(0.025, 0.975)), .cores = 16)
  b <- summarise_draws(x, .cores = 16)  
  bind_cols(b, a[,-1])
}
y <- SummaryFx(x)
saveRDS(y, "model_summaries.Rds")
saveRDS(y, "model_object.Rds")