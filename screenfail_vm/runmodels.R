#runmodels
library(cmdstanr)
library(dplyr)

set_cmdstan_path("/opt/cmdstanr/cmdstan-2.29.2/")
setwd("Documents/screenfail/")
mydf <- readRDS("model_with_data.Rds")

# count <- 0
# lapply(mydf$mods, function(x) {
#     print(count)
#     count <<- count + 1
#     mod <- cmdstan_model(write_stan_file(x$code))
#     fit <- mod$sample(data = x$data, num_cores = 4)
#     fit$save_object(file = paste0("mod", count, ".Rds"))
#   })

## full model
# 3, 4,5, 6, 7, 9
# for(i in seq_along(mydf$mods)) {
 # for(i in c(5)) {
 #  print(i)
 #  mod <- cmdstan_model(write_stan_file(mydf$mods_cmpr_cond[[i]]$code))
 #  fit <- mod$sample(data = mydf$mods_cmpr_cond[[i]]$data, num_cores = 4, adapt_delta = 0.99, max_treedepth = 15)
 #  fit$save_object(file = paste0("mod", i, ".Rds"))
 # }

for(i in c(6:9)) {
  print(i)
  mod <- cmdstan_model(write_stan_file(mydf$mods_cond[[i]]$code))
  fit <- mod$sample(data = mydf$mods_cond[[i]]$data, num_cores = 4, adapt_delta = 0.99, max_treedepth = 15)
  fit$save_object(file = paste0("mod_cond", i, ".Rds"))
}

# for(i in 5) {
#   print(i)
#   mod <- cmdstan_model(write_stan_file(mydf$mods_smpl[[i]]$code))
#   fit <- mod$sample(data = mydf$mods_smpl[[i]]$data, num_cores = 4, adapt_delta = 0.99, max_treedepth = 15)
#   fit$save_object(file = paste0("mod_smpl", i, ".Rds"))
# }