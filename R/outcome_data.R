set.seed(3592)

stratum_sample_sizes = c(20, 4, 5, 11)
total_sample_size = stratum_sample_sizes %>% sum()

trial_data = output$recruit_data %>%
  group_by(Stratum) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = stratum_sample_sizes) %>%
  mutate(samp = purrr::map2(.x = data, .y = n,
                            .f = function(.x, .y) slice_sample(.data = .x, n = .y))) %>%
  select(-data, -n) %>%
  unnest(samp) %>%
  mutate(
    trial = 1,
    treatment = rbinom(total_sample_size, 1, 0.8),
    outcome = ifelse(treatment == 1,
                     1*total + 2*pct_black_or_african_american + 3*pct_female + 4*pct_free_and_reduced_lunch,
                     4*total + 3*pct_black_or_african_american + 2*pct_female + 1*pct_free_and_reduced_lunch)
  ) %>%
  select(-Stratum, -rank, -ncessch)

selection_covariates = c("total", "pct_black_or_african_american", "pct_female", "pct_free_and_reduced_lunch")
