# Left truncation simulations
# Setting 3: single real world arm (subject to LT)

library(tidyverse)
library(survival)

library(doParallel)
registerDoParallel(cores = detectCores())
print(detectCores())

source("survmean.R") # RMST function from survival package
source("mackenzie.R") # Mackenzie survival estimation functions

CENSORING_RATE <- 0.5

GenerateFullDataset <- function(n,
                                entry_rate,
                                p_preindex,
                                dependence,
                                baseline_hazard,
                                censoring_rate = CENSORING_RATE) {
  # Args:
  #   n: sample size
  #   entry_rate: rate of entry event, calculated using CalculateEntryRate in order to satisfy given p_lt
  #   p_preindex: probability of entry before index date (i.e. entry = 0)
  #   dependence: parameter controlling dependence of survival time on entry time, on the log hazard ratio scale 
  #   baseline_hazard: baseline hazard of death, constant
  # Output: simulated dataset
  
  data <- data.frame(entry = rexp(n, rate = entry_rate)) %>%
    mutate(preindex_entry = rbinom(n, 1, p_preindex),
           entry = entry * (1 - preindex_entry),
           lambda = baseline_hazard * exp(dependence * entry),
           survival_time = rexp(n, lambda),
           censoring_time = rexp(n, lambda * censoring_rate / (1 - censoring_rate)),
           time = pmin(survival_time, censoring_time),
           status = as.numeric(censoring_time > survival_time),
           is_lt = (entry > time))

  return(data)
}


CalculateEntryRate <- function(n,
                               p_lt,
                               p_preindex,
                               dependence,
                               baseline_hazard,
                               nreps = 1000) {
  # Given a desired LT proportion with other simulation parameters, computes the entry_rate parameter
  # that will yield this.
  # Args:
  #   n: sample size
  #   p_lt: desired LT proportion, can range from (0, 1 - p_preindex)
  #   p_preindex: probability of entry before index date (i.e. entry = 0)
  #   dependence: parameter controlling dependence of survival time on entry time, on the log hazard ratio scale 
  #   baseline_hazard: baseline hazard of death, constant
  #   nreps: number of simulated datasets to compute empirical LT
  # Output: rate of entry event that will result in p_lt
  
  lt_fn <- function(er) {
    
    lt_obs <- replicate(nreps, 
                        mean(GenerateFullDataset(n, er, p_preindex, dependence, baseline_hazard)$is_lt))
    
    return(mean(lt_obs) - p_lt)
  }
  
  er_root <- tryCatch(uniroot(lt_fn, interval = c(1e-4, 10))$root, error = function(e) e)
  
  if (inherits(er_root, "error")) {
    er_root <- uniroot(lt_fn, interval = c(1e-8, 1e8))$root
  }
  
  return(er_root)
}


FitModels <- function(data, baseline_hazard) {
  # Given simulated dataset, fits model to un-truncated data to determine true parameter, and estimates based on 
  # truncated data.
  # Args:
  #   data: data.frame created by GenerateFullDataset
  #   baseline_hazard: baseline hazard of death, constant - used here to determine restriction time for RMST
  # Output: data.frame of model outputs
  
  # This choice of tau covers 90% of observed times for entry = 0 subjects and >90% for entry > 0 subjects
  rmst_tau <- -log(1 - 0.9) / baseline_hazard
  
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
                        km_formula <- as.formula(Surv(time, status) ~ 1)
                      } else {
                        km_formula <- as.formula(Surv(entry, time, status) ~ 1)
                      }
                      
                      model <- survfit(km_formula, 
                                       timefix = FALSE, 
                                       data = data) %>%
                        survmean(rmean = rmst_tau) %>%
                        t() %>%
                        as.data.frame %>%
                        mutate(type = type_,
                               n_used = nrow(data))
                      
                      return(model)
                    })
  
  models <- models %>% 
    rename(conf.low = "0.95LCL",
           conf.high = "0.95UCL")
  
  return(models)
}

FitMackenzieEsts <- function(data) {
  # Fits Mackenzie estimators of median survival andprobability of LT
  # Args:
  #   data: data.frame created by GenerateFullDataset
  # Output: data.frame of true parameters and Mackenzie estimates
  
  truth <- survfit(Surv(time, status) ~ 1, 
                   timefix = FALSE, 
                   data = data) %>%
    broom::glance() %>%
    transmute(median,
              prob_lt = mean(data$is_lt),
              type = "truth")
  
  data_lt <- data %>% 
    filter(!is_lt)
  
  median_est <- EstimateMedianSurvival(data_lt)
  prob_lt_est <- EstimatePropLT(data_lt)
    
  mackenzie <- tribble(~median,    ~prob_lt,    ~type,
                       median_est, prob_lt_est, "mackenzie")
  
  output <- bind_rows(truth, mackenzie)
  return(output)
}

ComputeMackenzieErrorMetrics <- function(mk_outs) {
  
  truth <- mk_outs %>%
    filter(type == "truth") %>% 
    transmute(median_truth = median, 
              prob_lt_truth = prob_lt,
              iter)
  
  metrics <- mk_outs %>%
    filter(type != "truth") %>% 
    select(median, prob_lt, iter) %>%
    inner_join(truth, by = "iter") %>% 
    mutate(median_rel_bias_i = median / median_truth,
           median_abs_bias_i = abs(median - median_truth),
           prob_lt_abs_bias_i = abs(prob_lt - prob_lt_truth)) %>% 
    summarise(median_rel_bias_mackenzie = mean(median_rel_bias_i, na.rm = TRUE),
              se_median_rel_bias_mackenzie = sd(median_rel_bias_i, na.rm = TRUE) / sqrt(sum(!is.na(median_rel_bias_i))),
              median_abs_bias_mackenzie = mean(median_abs_bias_i, na.rm = TRUE),
              se_median_abs_bias_mackenzie = sd(median_abs_bias_i, na.rm = TRUE) / sqrt(sum(!is.na(median_abs_bias_i))),
              prob_lt_abs_bias_mackenzie = mean(prob_lt_abs_bias_i, na.rm = TRUE),
              se_prob_lt_abs_bias_mackenzie = sd(prob_lt_abs_bias_i, na.rm = TRUE) / sqrt(sum(!is.na(prob_lt_abs_bias_i)))) 
    
  return(metrics)
}

ComputeErrorMetrics <- function(sim_models) {
  # Computes error metrics (percent bias and coverage) for estimators, compared to the truth based on the full
  # non-truncated dataset. 
  # Args:
  #   sim_models: data.frame containing regression model outputs for each model fit across simulation iterations
  # Output: single row data.frame containing percent bias and coverage for models fit to LT data 
  
  truth <- sim_models %>%
    filter(type == "truth") %>% 
    transmute(median_truth = median, 
              rmst_truth = rmean,
              iter)
  
  metrics <- sim_models %>%
    filter(type != "truth") %>% 
    select(type, rmean, median, conf.low, conf.high, iter, n_used) %>%
    inner_join(truth, by = "iter") %>% 
    mutate(median_rel_bias_i = median / median_truth,
           median_abs_bias_i = abs(median - median_truth),
           median_cover_i = 1*(median_truth >= conf.low & median_truth <= conf.high),
           rmst_rel_bias_i = rmean / rmst_truth,
           rmst_abs_bias_i = abs(rmean - rmst_truth)) %>% 
    group_by(type) %>%
    summarise(median_rel_bias = mean(median_rel_bias_i, na.rm = TRUE),
              se_median_rel_bias = sd(median_rel_bias_i, na.rm = TRUE) / sqrt(sum(!is.na(median_rel_bias_i))),
              median_abs_bias = mean(median_abs_bias_i, na.rm = TRUE),
              se_median_abs_bias = sd(median_abs_bias_i, na.rm = TRUE) / sqrt(sum(!is.na(median_abs_bias_i))),
              median_coverage = mean(median_cover_i, na.rm = TRUE),
              se_median_coverage = sd(median_cover_i, na.rm = TRUE) / sqrt(sum(!is.na(median_cover_i))),
              rmst_rel_bias = mean(rmst_rel_bias_i, na.rm = TRUE),
              se_rmst_rel_bias = sd(rmst_rel_bias_i, na.rm = TRUE) / sqrt(sum(!is.na(rmst_rel_bias_i))),
              rmst_abs_bias = mean(rmst_abs_bias_i, na.rm = TRUE),
              se_rmst_abs_bias = sd(rmst_abs_bias_i, na.rm = TRUE) / sqrt(sum(!is.na(rmst_abs_bias_i))),
              n_used = mean(n_used, na.rm = TRUE)) %>% 
    pivot_wider(id_cols = type,
                names_from = type,
                values_from = c(median_rel_bias, se_median_rel_bias, median_abs_bias, se_median_abs_bias, 
                                median_coverage, se_median_coverage, n_used,
                                rmst_rel_bias, se_rmst_rel_bias, rmst_abs_bias, se_rmst_abs_bias))
  
  return(metrics)  
}
  
RunSimulationSetting <- function(n,
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
  #   p_lt: desired LT proportion, can range from (0, 1 - p_preindex)
  #   p_preindex: probability of entry before index date (i.e. entry = 0)
  #   dependence: parameter controlling dependence of survival time on entry time, on the log hazard ratio scale 
  #   baseline_hazard: baseline hazard of death, constant
  #   niter: number of simulation iterations
  # Output: single row data.frame of sim parameters and results
  
  print(c(n, p_lt, p_preindex, dependence, baseline_hazard, niter))
  
  # Sample size that would be observed without any LT, so that expected observed sample size is the original n
  n <- n / (1 - p_lt)
  
  sim_params <- data.frame(n = n,
                           p_lt = p_lt,
                           p_preindex = p_preindex, 
                           hr_dep = round(exp(dependence), 2), 
                           baseline_hazard = baseline_hazard)
  
  # Set rate of entry in order to satisfy given p_lt
  entry_rate <- CalculateEntryRate(n, p_lt, p_preindex, dependence, baseline_hazard)
  
  sim_models <- NULL
  mk_outs <- NULL
  for (i in 1:niter) {
    data_i <- GenerateFullDataset(n, entry_rate, p_preindex, dependence, baseline_hazard)
    fits_i <- tryCatch(FitModels(data_i, baseline_hazard), error = function(e) e)
    mk_i <- tryCatch(FitMackenzieEsts(data_i), error = function(e) e)
      
    if (!inherits(fits_i, "error") & !inherits(mk_i, "error")) {
      fits_i <- fits_i %>% mutate(iter = i)
      mk_i <- mk_i %>% mutate(iter = i)
      
      sim_models[[i]] <- fits_i
      mk_outs[[i]] <- mk_i
    }
  }
  sim_models <- sim_models %>% bind_rows()
  mk_outs <- mk_outs %>% bind_rows()
    
  out <- bind_cols(sim_params, 
                   ComputeErrorMetrics(sim_models), 
                   ComputeMackenzieErrorMetrics(mk_outs))
  return(out)
}
  

RunGridSimsParallel <- function(grid) {
  # Runs simulations and collects results over a parameter grid, with parallel loop
  # Args:
  #   grid: a data.frame with columns (n, p_lt, p_preindex, dependence, baseline_hazard, niter) with each row 
  #         corresponding to a different configuration
  # Output: data.frame that collects simulation results as returned by RunSimulationSetting for each grid config,
  #         and then transforms to long format
  
  GridSim <- function(params) {
    RunSimulationSetting(n = params$n, 
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
    pivot_longer(cols = -c(n, p_lt, p_preindex, hr_dep, baseline_hazard),
                 names_to = c("metric", "model"),
                 names_pattern = "(median_rel_bias|median_abs_bias|rmst_rel_bias|rmst_abs_bias|se_median_rel_bias|se_median_abs_bias|se_rmst_rel_bias|se_rmst_abs_bias|se_median_coverage|median_coverage|prob_lt_abs_bias|se_prob_lt_abs_bias|n_used)_(.+)",
                 values_to = "value")
  
  return(results)
}


if (T) {
  grid <- expand.grid(n = c(300, 1000),
                      p_lt = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                      p_preindex = c(0.2),
                      dependence = log(c(1, 1.01, 1.05, 1.10)),
                      baseline_hazard = c(1 / 12), # time scale is in months
                      niter = 500)
  
  results <- RunGridSimsParallel(grid)
  save(results, file = "SimResults_singlerw.RData")
}

