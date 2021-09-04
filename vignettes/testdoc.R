library(devtools)
library(tidyverse)
library(generalizeRdata)
library(generalizeR)
library(ggthemes)
library(gridExtra)
library(usmap)
library(janitor)
library(grid)
library(cluster)
library(ClusterR)
library(crayon)
library(kableExtra)
library(crayon)
library(tidyverse)
library(viridisLite)
library(viridis)
library(ggrepel)
library(ggnewscale)

# 1. Get your inference population (technically step 2)
selection_vars <- c("pct_female", "pct_black_or_african_american", "pct_free_and_reduced_lunch")
id_vars <- "ncessch"

inference_pop <- cc %>%
  filter(charter == 1) %>%
  filter(g_10_offered == "Yes") %>%
  select(ncessch, all_of(selection_vars), all_of(id_vars))

inference_pop <- inference_pop %>%
  mutate(factor = as.factor(sample(1:5, nrow(inference_pop), replace = T)))

#2. let's generate a simple sample from this

test <- stratify(inference_pop, guided = T)
stratify_results <- stratify(inference_pop, guided = F, n_strata = 4, variables = selection_vars, idnum = id_vars,
                             verbose = F)








# assess test -------------------------------------------------------------



sample <- tibble(ncessch = c(stratify_results$recruitment_lists[[1]]$ncessch[1:20],
                             stratify_results$recruitment_lists[[2]]$ncessch[1:4],
                             stratify_results$recruitment_lists[[3]]$ncessch[1:6],
                             stratify_results$recruitment_lists[[4]]$ncessch[1:11]))

#3 let's try assess_wrap

gen_results <- assess_wrap(sample,inference_pop,join_var = "ncessch")

# sample_trial <- sample %>% mutate(trial = 1) %>%
#   right_join(inference_pop) %>%
#   mutate(trial = ifelse(is.na(trial),0,1))
#
# assess_object <- assess(trial = "trial", selection_covariates = selection_vars,
#                         data = sample_trial, is_data_disjoint = TRUE)




