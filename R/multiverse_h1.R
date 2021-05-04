

library(tidyverse)
library(parameters)

simit <- function(targetedness_tp, start_date, toxic_type) {
  
  text_type <- "transcr"
  print(start_date)
  print(targetedness_tp)
  
  analysis_dat <- readRDS(paste0("robust/", text_type, "/", start_date, "/log_dat_", start_date, ".rds"))
  
  if(toxic_type == "continous"){
    analysis_dat$toxic <- analysis_dat$TOXICITY
  }
  
  if(targetedness_tp == "+1m"){
    
    analysis_dat <- analysis_dat %>% 
      mutate(targetedness_gen = targetedness_hl_1m)
    
  } else if (targetedness_tp == "+500k"){
    
    analysis_dat <- analysis_dat %>% 
      mutate(targetedness_gen = targetedness_hl_500k)
    
  } else if (targetedness_tp == "+100k"){
    
    analysis_dat <- analysis_dat %>% 
      mutate(targetedness_gen = targetedness_hl_100k)
    
  } else if (targetedness_tp == "+50k"){
    
    analysis_dat <- analysis_dat %>% 
      mutate(targetedness_gen = targetedness_hl_50k)
    
  } else if (targetedness_tp == "+10k"){
    
    analysis_dat <- analysis_dat %>% 
      mutate(targetedness_gen = targetedness_hl_10k)
    
  } else if (targetedness_tp == "+5k"){
    
    analysis_dat <- analysis_dat %>% 
      mutate(targetedness_gen = targetedness_hl_5k)
    
  } else if (targetedness_tp == "+1k"){
    
    analysis_dat <- analysis_dat %>% 
      mutate(targetedness_gen = targetedness_hl_1k)
    
  } 
  
  
  
  
  
  if(targetedness_tp != "numeric"){
    mod1_glm <- glm(targetedness_gen ~ toxic, data = analysis_dat, family = binomial())
    
    mod2_glm <- glm(targetedness_gen ~ toxic  + affil  + advertiser_binary + age + gender + max_reached_log + 
                      add_runs_log + num_days_log + closeness_to_election_log + spd_lower_bound_log,
                    data = analysis_dat, family = binomial())
    
    mod3_glm <- glm(targetedness_gen ~ advertiser_binary*toxic  + affil  + age + gender + max_reached_log + 
                      add_runs_log + num_days_log + closeness_to_election_log + spd_lower_bound_log, 
                    data = analysis_dat, family = binomial())
    
    mod4_glm <- glm(targetedness_gen ~ affil*toxic  + advertiser_binary + age + gender + max_reached_log + 
                      add_runs_log + num_days_log + closeness_to_election_log + spd_lower_bound_log, 
                    data = analysis_dat, family = binomial())
    
    
    
    
    glm1 <- parameters::model_parameters(mod1_glm,
                                         df_method="wald",
                                         exponentiate = T) %>%
      mutate(mod = "mod1")
    
    glm2 <- parameters::model_parameters(mod2_glm,
                                         df_method="wald",
                                         exponentiate = T) %>%
      mutate(mod = "mod2")
    
    glm3 <- parameters::model_parameters(mod3_glm,
                                         df_method="wald",
                                         exponentiate = T) %>%
      mutate(mod = "mod3")
    
    glm4 <- parameters::model_parameters(mod4_glm,
                                         df_method="wald",
                                         exponentiate = T) %>%
      mutate(mod = "mod4")
    
    
    glm_dat <- glm1 %>%
      bind_rows(glm2) %>%
      bind_rows(glm3) %>%
      bind_rows(glm4) %>%
      mutate(mod_type = "glm") %>%
      mutate(targetedness_type = targetedness_tp) %>%
      mutate(tox_type = "continous") %>%
      mutate(text_type = "transcription") %>%
      mutate(hypothesis = "h1")
    
    
    saveRDS(glm_dat, file = paste0("robust/fin/", start_date, "/glm_", targetedness_tp ,"_", toxic_type, "_transcr.rds"))
    
  } else {
    
    mod1_lm <- lm(targetedness_num ~ toxic, data = analysis_dat)
    
    mod2_lm <- lm(targetedness_num ~ toxic + affil  + advertiser_binary + age + gender + max_reached_log + 
                    add_runs_log + num_days_log + closeness_to_election_log + spd_lower_bound_log, 
                  data = analysis_dat)
    
    mod3_lm <- lm(targetedness_num ~ advertiser_binary*toxic  + affil +  
                    max_reached_log + age +
                    gender + add_runs_log  + closeness_to_election_log +
                    num_days_log + spd_lower_bound_log, 
                  data = analysis_dat)
    
    mod4_lm <- lm(targetedness_num ~ affil*toxic + 
                    max_reached_log + age + advertiser_binary+ 
                    gender + add_runs_log  + closeness_to_election_log +
                    num_days_log + spd_lower_bound_log,
                  data = analysis_dat)
    
    
    
    
    lm1 <- parameters::model_parameters(mod1_lm, standardize = "refit") %>%
      mutate(mod = "mod1")
    
    lm2 <- parameters::model_parameters(mod2_lm, standardize = "refit") %>%
      mutate(mod = "mod2")
    
    lm3 <- parameters::model_parameters(mod3_lm, standardize = "refit") %>%
      mutate(mod = "mod3")
    
    lm4 <- parameters::model_parameters(mod4_lm, standardize = "refit") %>%
      mutate(mod = "mod4")
    
    
    lm_dat <- lm1 %>%
      bind_rows(lm2) %>%
      bind_rows(lm3) %>%
      bind_rows(lm4) %>%
      mutate(mod_type = "lm") %>%
      mutate(targetedness_type = targetedness_tp) %>%
      mutate(tox_type = "continous") %>%
      mutate(text_type = "transcription") %>%
      mutate(hypothesis = "h1")
    
    saveRDS(lm_dat, file = paste0("robust/fin/", start_date, "/lm_", targetedness_tp ,"_", toxic_type, "_transcr.rds"))
    
  }
  
}



library(furrr)

plan(multiprocess, workers = 4)

targetedness_tp <- c("+1m", "+500k", "+100k", "+50k", "+10k", "+5k", "+1k", "numeric")
dates <- c("2020-04-03", "2020-05-03", "2020-06-03", "2020-07-03", "2020-08-03", "2020-09-03", "2020-10-03")
toxic_type <- c("binary", "continous")

expand_grid(targetedness_tp, dates) %>% 
  split(1:nrow(.)) %>%
  future_walk(~{simit(.x$targetedness_tp, .x$dates, .x$toxic_type)})

