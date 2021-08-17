# Left truncation simulations
# Setting 1: binary treatment with both arms subject to LT

library(tidyverse)
library(survival)

library(doParallel)
registerDoParallel(cores = detectCores())
print(detectCores())

CENSORING_RATE <- 0.5
BETA <- log(0.8)

GenerateFullDataset <- function(n,
                                p_trt,
                                entry_rate,
                                p_preindex,
                                dependence,
                                baseline_hazard,
                                beta = BETA,
                                censoring_rate = CENSORING_RATE) {
  # Args:
  #   n: sample size
  #   p_trt: P(trt = 1)
  #   entry_rate: rate of entry event, calculated using CalculateEntryRate in order to satisfy given p_lt
  #   p_preindex: probability of entry before index date (i.e. entry = 0)
  #   dependence: parameter controlling dependence of survival time on entry time, on the log hazard ratio scale 
  #   baseline_hazard: baseline hazard of death, constant
  # Output: simulated dataset
  
  data <- data.frame(entry = rexp(n, rate = entry_rate)) %>%
    mutate(preindex_entry = rbinom(n, 1, p_preindex),
           entry = entry * (1 - preindex_entry),
           trt = rbinom(n, 1, p_trt),
           lambda = baseline_hazard * exp(beta * trt + dependence * entry),
           survival_time = rexp(n, lambda),
           censoring_time = rexp(n, lambda * censoring_rate / (1 - censoring_rate)),
           time = pmin(survival_time, censoring_time),
           status = as.numeric(censoring_time > survival_time),
           is_lt = (entry > time))
  
  return(data)
}


CalculateEntryRate <- function(n,
                               p_trt,
                               p_lt,
                               p_preindex,
                               dependence,
                               baseline_hazard,
                               nreps = 1000) {
  # Given a desired LT proportion with other simulation parameters, computes the entry_rate parameter
  # that will yield this.
  # Args:
  #   n: sample size
  #   p_trt: P(trt = 1)
  #   p_lt: desired LT proportion, can range from (0, 1 - p_preindex)
  #   p_preindex: probability of entry before index date (i.e. entry = 0)
  #   dependence: parameter controlling dependence of survival time on entry time, on the log hazard ratio scale 
  #   baseline_hazard: baseline hazard of death, constant
  #   nreps: number of simulated datasets to compute empirical LT
  # Output: rate of entry event that will result in p_lt
  
  lt_fn <- function(er) {
    
    lt_obs <- replicate(nreps, 
                        mean(GenerateFullDataset(n, p_trt, er, p_preindex, dependence, baseline_hazard)$is_lt))
    
    return(mean(lt_obs) - p_lt)
  }

  er_root <- tryCatch(uniroot(lt_fn, interval = c(1e-4, 10))$root, error = function(e) e)
  
  if (inherits(er_root, "error")) {
    er_root <- uniroot(lt_fn, interval = c(1e-8, 1e8))$root
  }
  
  return(er_root)
}

FitModels <- function(data) {
  # Given simulated dataset, fits model to un-truncated data to determine true parameter, and estimates based on 
  # truncated data.
  # Args:
  #   data: data.frame created by GenerateFullDataset
  # Output: data.frame of model outputs
  
  truth <- data %>% 
    mutate(type = "truth")
  
  naive <- data %>% 
    filter(!is_lt) %>% 
    mutate(type = "naive")
  
  risk_set_adjustment <- naive %>% 
    mutate(type = "risk_set_adjustment")
  
  baseline_dfs <- map_dfr(c(0.25, .5, .75),
                          function(q) {
                            naive %>% 
                              filter(entry <= quantile(entry, q)) %>% 
                              mutate(type = paste0("baseline", as.character(q * 100)))
                          })
  
  all_data <- list(truth, 
                   naive,
                   risk_set_adjustment,
                   baseline_dfs) %>% 
    bind_rows()
  
  models <- map_dfr(unique(all_data$type),
                    function(type_) {
                      data <- filter(all_data, type == type_)
                      
                      if (type_ %in% c("truth", "naive")) {
                        cox_formula <- as.formula(Surv(time, status) ~ trt)
                      } else {
                        cox_formula <- as.formula(Surv(entry, time, status) ~ trt)
                      }
                      
                      model <- coxph(cox_formula, 
                                     control = coxph.control(timefix = FALSE, iter.max = 1000), 
                                     data = data) %>% 
                        broom::tidy(conf.int = TRUE) %>% 
                        mutate(type = type_,
                               n_used = nrow(data))
                      
                      return(model)
                    })
  
  return(models)
}


ComputeErrorMetrics <- function(sim_models) {
  # Computes error metrics (percent bias and coverage) for estimators, compared to the truth based on the full
  # non-truncated dataset. 
  # Args:
  #   sim_models: data.frame containing regression model outputs for each model fit across simulation iterations
  # Output: single row data.frame containing percent bias and coverage for models fit to LT data 
  
  truth <- sim_models %>%
    filter(type == "truth") %>% 
    transmute(truth = estimate, 
              iter)
  
  metrics <- sim_models %>%
    filter(type != "truth") %>% 
    select(type, estimate, conf.low, conf.high, iter, n_used) %>%
    inner_join(truth, by = "iter") %>% 
    mutate(rel_bias_i = exp(estimate) / exp(truth),
           abs_bias_i = abs(exp(estimate) - exp(truth)),
           cover_i = 1*(truth >= conf.low & truth <= conf.high)) %>% 
    group_by(type) %>%
    summarise(rel_bias = mean(rel_bias_i, na.rm = TRUE),
              se_rel_bias = sd(rel_bias_i, na.rm = TRUE) / sqrt(sum(!is.na(rel_bias_i))),
              abs_bias = mean(abs_bias_i, na.rm = TRUE),
              se_abs_bias = sd(abs_bias_i, na.rm = TRUE) / sqrt(sum(!is.na(abs_bias_i))),
              coverage = mean(cover_i, na.rm = TRUE),
              se_coverage = sd(cover_i, na.rm = TRUE) / sqrt(sum(!is.na(cover_i))),
              n_used = mean(n_used, na.rm = TRUE)) %>% 
    pivot_wider(id_cols = type,
                names_from = type,
                values_from = c(rel_bias, se_rel_bias, abs_bias, se_abs_bias, coverage, se_coverage, n_used))
  
  return(metrics)  
}

  
RunSimulationSetting <- function(n,
                                 p_trt,
                                 p_lt,
                                 p_preindex,
                                 dependence,
                                 baseline_hazard,
                                 niter = 500) {
  # Given simulation config parameters, runs niter simulation iterations as follows:
  # 1) Determine required entry_rate parameter to produce p_lt using CalculateEntryRate()
  # 2) Repeat niter times:
  #      * Generate data with GenerateFullDataset()
  #      * Compute true parameter and fit estimators with FitModels()
  # 3) Combine into single row data.frame with ComputeErrorMetrics()
  # Args:
  #   n: expected sample size observed; transformed to sample size that would be observed without any LT
  #   p_trt: P(trt = 1)
  #   p_lt: desired LT proportion, can range from (0, 1 - p_preindex)
  #   p_preindex: probability of entry before index date (i.e. entry = 0)
  #   dependence: parameter controlling dependence of survival time on entry time, on the log hazard ratio scale 
  #   baseline_hazard: baseline hazard of death, constant
  #   niter: number of simulation iterations
  # Output: single row data.frame of sim parameters and results
  
  # Sample size that would be observed without any LT, so that expected observed sample size is the original n
  n <- n / (1 - p_lt)
  
  sim_params <- data.frame(n = n,
                           p_trt = p_trt,
                           p_lt = p_lt,
                           p_preindex = p_preindex, 
                           hr_dep = round(exp(dependence), 2), 
                           baseline_hazard = baseline_hazard)
  
  # Set rate of entry in order to satisfy given p_lt
  entry_rate <- CalculateEntryRate(n, p_trt, p_lt, p_preindex, dependence, baseline_hazard)
  
  sim_models <- NULL
  for (i in 1:niter) {
    data_i <- GenerateFullDataset(n, p_trt, entry_rate, p_preindex, dependence, baseline_hazard)
    fits_i <- tryCatch(FitModels(data_i), error = function(e) e)
    
    if (!inherits(fits_i, "error")) {
      fits_i <- fits_i %>% 
        mutate(iter = i)
      sim_models[[i]] <- fits_i
    }
  }
  sim_models <- sim_models %>%
    bind_rows()
  
  out <- bind_cols(sim_params, ComputeErrorMetrics(sim_models))
  return(out)
}


RunGridSimsParallel <- function(grid) {
  # Runs simulations and collects results over a parameter grid, with parallel loop
  # Args:
  #   grid: a data.frame with columns (n, p_trt, p_lt, p_preindex, dependence, baseline_hazard, niter) with each row 
  #         corresponding to a different configuration
  # Output: data.frame that collects simulation results as returned by RunSimulationSetting for each grid config,
  #         and then transforms to long format
  
  GridSim <- function(params) {
    RunSimulationSetting(n = params$n, 
                         p_trt = params$p_trt, 
                         p_lt = params$p_lt, 
                         p_preindex = params$p_preindex,
                         dependence = params$dependence, 
                         baseline_hazard = params$baseline_hazard,
                         niter = params$niter)
  }
  
  # Filter out impossible sim configurations
  grid <- grid %>%
    filter(p_preindex + p_lt <= 0.9)
  
  results <- foreach(ix = 1:nrow(grid)) %dopar% {
    print(grid[ix,])
    GridSim(grid[ix,])
  }
  
  results <- results %>%
    bind_rows() %>%
    pivot_longer(cols = -c(n, p_trt, p_lt, p_preindex, hr_dep, baseline_hazard),
                 names_to = c("metric", "model"),
                 names_pattern = "(rel_bias|abs_bias|coverage|n_used|se_rel_bias|se_abs_bias|se_coverage)_(.+)",
                 values_to = "value")
  
  return(results)
}


if (T) {
  grid <- expand.grid(n = c(300, 1000),
                      p_trt = 0.5,
                      p_lt = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                      p_preindex = c(0.2),
                      dependence = log(c(1, 1.01, 1.05, 1.10)),
                      baseline_hazard = c(1 / 12), # time scale is in months
                      niter = 500)
  
  results <- RunGridSimsParallel(grid)
  save(results, file = "SimResults_allrw.RData")
}

