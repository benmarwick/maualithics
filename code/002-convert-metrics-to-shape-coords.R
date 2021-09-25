# here we import the CSV file of caliper measurements, then compute on the
# measurements to generate coordinates for a shape outline

library(tidyverse)

# import
ma_flakes_tr <- read_csv(here::here("data/ma_flakes_tr.csv"))

# custom function is from here, we have a local version in 999-functions.R
# so we can work offline if we need to
# https://gist.github.com/benmarwick/fac5754131c1fcc14d32b1658215f0e3

source(here::here("code/999-functions.R"))

polys <- polygon_shape_by_variables(data = ma_flakes_tr,       # our data frame
                                    x = 'width_at_0.50_length',                # col for coord for x-axis of poly centers
                                    y = 'max_dimension',       # col for coord for y-axis of poly centers

                                    sf = 10,    # scale factor
                                    top_horizontal_measurement =    'platform_width', # cols to use to compute poly verts
                                    main_vertical_measurement =     'percussion_length',
                                    upper_horizontal_measurement =  'width_at_0.25_length',
                                    middle_horizontal_measurement = 'width_at_0.50_length',
                                    bottom_horizontal_measurement = 'width_at_0.75_length')

# take a look to ensure we are on the right track
p <- ggplot() +
  geom_polygon(data=polys,
               aes(x=x,
                   y=y,
                   group = fill),
               colour="black",
               alpha = 0.5) +
  coord_equal() +
  theme_bw()
p

# make ready to convert to TPS format (thin plate spline)
polys_lmk <-
polys %>%
  group_by(fill) %>%
  mutate(landmark.number = row_number(fill),
         specimen.ids  = fill) %>%
  ungroup() %>%
  select(-fill)

# export as TPS file
create_tps(polys_lmk,
           here::here("data/landmarks.tps"))

