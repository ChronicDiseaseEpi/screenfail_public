## 01_trial_level_data.R

library(tidyverse)

aact <- readRDS("Data/aact_extract.Rds")
gsk <- read_csv("Data/gsk_sf_null.csv")
vivli <- read_csv("Data/vivli_sf_null.csv")

## first get list of included trials from baseline tables
gsk <- gsk %>%
  select(nct_id, n_rand, n_fail) %>%
  rename(n_sf = n_fail) %>%
  mutate(source = "gsk")
vivli <- vivli %>%
  select(nct_id, n_rand, n_sf) %>%
  mutate(source = "vivli")
trials <- bind_rows(gsk, vivli)

# create tables of relevant variables
fac <- trials %>% 
  inner_join(aact$facilities) %>%
  distinct(city, country, nct_id) %>%
  count(nct_id) %>% 
  rename(facilities = n)

cntr <- trials %>% 
  inner_join(aact$facilities) %>%
  distinct(country, nct_id) %>%
  count(nct_id) %>% 
  rename(countries = n)

mask <- trials %>% 
  inner_join(aact$designs) %>%
  select(nct_id, masking)

studies <- trials %>%
  select(nct_id) %>%
  inner_join(aact$studies %>% 
               select(nct_id, start_date, completion_date, phase, enrollment, number_of_arms))

# combine interesting tables
table <- trials %>%
  left_join(studies) %>%
  left_join(fac) %>%
  left_join(cntr) %>%
  left_join(mask)
  
table1 <- table %>%
  mutate(duration_yrs = round(as.numeric(difftime(completion_date, start_date, units = "days")/365.25), digits = 1)) %>%
  mutate(n_screen = n_rand + n_sf) %>%
  mutate(pc_sf = round(n_sf/n_screen*100, digits = 1)) %>%
  select(-c(n_rand)) %>% 
  arrange(pc_sf)

write_csv(table1, "Outputs/trial_details.csv")


mydf <- read_csv("Outputs/trial_details.csv")


plot1 <- ggplot(mydf, aes(x = start_date, pc_sf, size = n_screen, colour = phase)) +
  geom_point() +
  scale_x_date("Year trial started") +
  scale_y_continuous("Percentage screen failure") +
  scale_colour_manual("Phase", values = c("purple", "red", "blue")) + 
  scale_size("Individuals screened") 
plot1
pdf("Outputs/Figure1.pdf")
plot1
dev.off()
