inference_pop = cc %>%
  filter(st == "TX") %>%
  filter(charter == 1) %>%
  filter(g_10_offered == "Yes")

selection_covariates = c("total", "pct_black_or_african_american", "pct_female")

output = stratify(inference_pop, guided = FALSE, n_strata = 4, variables = selection_covariates, idnum = "ncessch")

set.seed(3592)

stratum_sample_sizes = c(1, 37, 1, 1) #c(20, 4, 5, 11)
total_sample_size = stratum_sample_sizes %>% sum()

beta = 0
delta = c(0.5, 1, 1, 1, 0) * 90295422

target_data = output$recruit_data

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
           ncessch %in% trial_data$ncessch ~ 0)
         ) %>%
  dplyr::select(-Stratum, -rank, -ncessch)

total_mean = target_data %>%
  filter(trial == 1) %>%
  pull(total) %>%
  mean()

total_sd = target_data %>%
  filter(trial == 1) %>%
  pull(total) %>%
  sd()

epsilon = rnorm(length(target_data$total), 0, 10^4)

target_data = target_data %>%
  mutate(total = (total-total_mean)/total_sd,
         outcome = delta[1]*treatment + (beta + delta[2]*treatment)*total +
                                        (beta + delta[3]*treatment)*pct_black_or_african_american +
                   epsilon
         )

control = target_data %>%
  filter(treatment == 0)

target_data %>%
  filter(treatment == 0) %>%
  pull(outcome) %>%
  var()

control_var = control %>%
  pull(outcome) %>%
  var()

target_data = target_data %>%
  mutate(outcome = outcome / control_var)

# mean in control group should be close to 0
# mean in treatment group should be close to 0.5 (change delta[1])

target_data %>%
  filter(treatment == 1) %>%
  pull(outcome) %>%
  mean()

test = weighting("outcome", "treatment", "trial", selection_covariates, target_data, is_data_disjoint = FALSE)
hist(test$weights[test$weights != 0])

assess("trial", selection_covariates, target_data, is_data_disjoint = FALSE)

treatment_mean = target_data %>%
  filter(treatment == 1) %>%
  select(outcome) %>%
  pull(outcome) %>%
  mean()

treatment_var = target_data %>%
  filter(treatment == 1) %>%
  select(outcome) %>%
  pull(outcome) %>%
  var()

control_mean = target_data %>%
  filter(treatment == 0) %>%
  select(outcome) %>%
  pull(outcome) %>%
  mean()

control_var = target_data %>%
  filter(treatment == 0) %>%
  select(outcome) %>%
  pull(outcome) %>%
  var()



fit = lm(outcome ~ total + pct_black_or_african_american + pct_female + pct_free_and_reduced_lunch, data = control)
fit %>% summary()


