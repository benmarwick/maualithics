library(tidyverse)

# import the data
ma_flakes <- readxl::read_excel(here::here("data/UW Fieldschool Mau A 2015 post_excavation_season_April_2016.xlsx"),
                                col_names = FALSE)

# transpose so that columns are variables, rows are individual artefacts
# and set col data types so we can compute on them
ma_flakes_tr <-
ma_flakes %>%
  pivot_longer(-`...1`)  %>%
  pivot_wider(names_from = `...1`,
              values_from = value) %>%
  mutate_all(parse_guess) %>%
  # fix a few bad guesses
  mutate(across(c(mass,
                thickness_at_0.25_width,
                thickness_at_0.75_width,
                dorsal_cortex_perc,
                dorsal_scar_initiation_count),
         parse_number))

# export this data frame to use with our analysis
write_csv(ma_flakes_tr,
          here::here("data/ma_flakes_tr.csv"))

# explore dorsal cortex distributions, as one possible varible
# to divide the assemblage into groups, we need to explore a few
# Weight could be another?
ggplot(ma_flakes_tr) +
  aes(dorsal_cortex_perc) +
  geom_histogram()




