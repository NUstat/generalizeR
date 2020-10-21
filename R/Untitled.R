cc
cc$level

county_data <- read.csv("/Users/katiecoburn/Dropbox/Generalizer Project/data/grf19_lea_county.csv") %>%
  tibble() %>%
  select(LEAID, NAME_COUNTY19) %>%
  clean_names()

test <- full_join(county_data, cc, by = "leaid") %>%
  mutate(ncessch = factor(ncessch))

new_data <- distinct(test, ncessch, .keep_all = TRUE)

