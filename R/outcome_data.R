inference_pop = cc %>%
  filter(st == "TX") %>%
  filter(charter == 1) %>%
  filter(g_10_offered == "Yes")

selection_covariates = c("total", "pct_black_or_african_american", "pct_female", "pct_free_and_reduced_lunch")

output = stratify(inference_pop, guided = FALSE, n_strata = 4, variables = selection_covariates, idnum = "ncessch")

#set.seed(3592)

stratum_sample_sizes = c(20, 4, 5, 11)
total_sample_size = stratum_sample_sizes %>% sum()

beta = c(0, 4, 8, 12, 16)*4
delta = c(0.5, 1, 2, 3, 4)*4

target_data = output$recruit_data %>%
  mutate(
    total = scale(total),
    pct_female = scale(pct_female),
    pct_black_or_african_american = scale(pct_black_or_african_american),
    pct_free_and_reduced_lunch = scale(pct_free_and_reduced_lunch)
  )

trial_data = target_data %>%
  group_by(Stratum) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = stratum_sample_sizes) %>%
  mutate(samp = purrr::map2(.x = data, .y = n,
                            .f = function(.x, .y) slice_sample(.data = .x, n = .y))) %>%
  dplyr::select(-data, -n) %>%
  unnest(samp)

treatment_data = trial_data %>%
  group_by(Stratum) %>%
  nest() %>%
  ungroup() %>%
  mutate(samp = purrr::map(.x = data,
                            .f = function(.x) slice_sample(.data = .x, prop = 0.5))) %>%
  dplyr::select(-data) %>%
  unnest(samp)

target_data = target_data %>%
  mutate(trial = ifelse(ncessch %in% trial_data$ncessch, 1, 0),
         treatment = case_when(
           ncessch %in% treatment_data$ncessch ~ 1,
           ncessch %in% trial_data$ncessch ~ 0
           ),
         outcome = beta[1] + delta[1]*treatment +
                   (beta[2] + delta[2]*treatment)*total +
                   (beta[3] + delta[3]*treatment)*pct_female +
                   (beta[4] + delta[4]*treatment)*pct_black_or_african_american +
                   (beta[5] + delta[5]*treatment)*pct_free_and_reduced_lunch +
                   rnorm(1, 0, 10)*(1+treatment)
         ) %>%
  dplyr::select(-Stratum, -rank, -ncessch)

weighting("outcome", "treatment", "trial", selection_covariates, target_data)
