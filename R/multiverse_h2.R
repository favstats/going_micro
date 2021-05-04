

library(tidyverse)
library(parameters)

simit <- function(start_date, toxic_type) {
  
  
  text_type <- "transcr"
  print(start_date)
  
  
  
  analysis_dat <- readRDS(paste0("robust/", text_type, "/", start_date, "/log_dat_", start_date, ".rds"))
  
  
  
  if(toxic_type != "continous"){
    
    mod1_glm3 <- glm(toxic ~ advertiser_binary, data = campaign_dat, family = binomial())
    
    mod2_glm3 <- glm(toxic ~ affil  + age + gender + max_reached_log +# targetedness_num +
                       add_runs_log + advertiser_binary  + num_days_log+ 
                       closeness_to_election_log +
                       spd_lower_bound_log , data = campaign_dat, family = binomial())
    
    
    
    
    
    glm1 <- parameters::model_parameters(mod1_glm3,
                                         df_method="wald",
                                         exponentiate = T) %>%
      mutate(mod = "mod1")
    
    glm2 <- parameters::model_parameters(mod2_glm3,
                                         df_method="wald",
                                         exponentiate = T) %>%
      mutate(mod = "mod2")
    
    
    
    glm_dat <- glm1 %>%
      bind_rows(glm2) %>%
      mutate(mod_type = "glm") %>%
      # mutate(targetedness_type = targetedness_tp) %>%
      mutate(tox_type = "binary") %>%
      mutate(text_type = "transcription") %>%
      mutate(hypothesis = "h2")
    
    
    saveRDS(glm_dat, file = paste0("robust/fin/", start_date, "/glm_", "h2" ,"_binary_transcr.rds"))
    
  } else {
    
    analysis_dat$toxic <- analysis_dat$TOXICITY
    
    mod1_lm3 <- lm(toxic ~ advertiser_binary, data = campaign_dat)
    
    mod2_lm3 <- lm(toxic ~ affil  + age + gender + max_reached_log +# targetedness_num +
                     add_runs_log + advertiser_binary  + num_days_log+ 
                     closeness_to_election_log +
                     spd_lower_bound_log , data = campaign_dat)
    
    
    
    
    lm1 <- parameters::model_parameters(mod1_lm3, standardize = "refit") %>%
      mutate(mod = "mod1")
    
    lm2 <- parameters::model_parameters(mod2_lm3, standardize = "refit") %>%
      mutate(mod = "mod2")
    
    
    lm_dat <- lm1 %>%
      bind_rows(lm2) %>%
      mutate(mod_type = "lm") %>%
      # mutate(targetedness_type = targetedness_tp) %>%
      mutate(tox_type = "continous") %>%
      mutate(text_type = "transcription") %>%
      mutate(hypothesis = "h2")
    
    
    saveRDS(lm_dat, file = paste0("robust/fin/", start_date, "/lm_", "h2" ,"_continous_transcr.rds"))
    
  }
  
}


library(furrr)
plan(multiprocess, workers = 4)

dates <- c("2020-04-03", "2020-05-03", "2020-06-03", "2020-07-03", "2020-08-03", "2020-09-03", "2020-10-03")
toxic_type <- c("continous", "binary")


expand_grid(toxic_type, dates) %>% 
  split(1:nrow(.)) %>% 
  future_walk(~{simit(.x$dates, .x$toxic_type)})

