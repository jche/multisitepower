
# get cluster results down into a single file

require(tidyverse)
require(glue)

dir <- "case_study"
fhead <- "res2_"

all_files <- list.files(dir)

max_runID <- 0
all_dfs <- all_files[str_starts(all_files, glue("{fhead}"))] %>%
  map(function(fname) {
    df <- read_csv(glue("{dir}/{fname}")) %>%
      mutate(runID = runID + max_runID)
    max_runID <<- max(df$runID)
    df
  })
res <- bind_rows(all_dfs)

rm(all_files)
rm(all_dfs)

if (T) {
  foo <- res %>% 
    mutate(runID = rep(1:1000, each=60))
  write_csv(foo, "case_study/case_study2_results.csv")
}

if (F) {
  # for full_study: subset into smaller studies
  #  - have already: J=25, ICC=tau=0.2, tx_sd=0.2, 0.3
  
  hits <- res %>% 
    filter(J==25, 
           ICC==0.2,
           tau==0.2)
           # tx_sd==0.2)
  
  # Get ATEs per method (long data)
  ATEhats <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("ATEhat"), names_to = "method", values_to = "ATEhat") %>%
    dplyr::select(ATE, sid, method, ATEhat) %>%
    mutate(method = str_sub(method, 8)) %>%
    ungroup()
  
  # Get SEs per method
  SEs <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("SE_"), names_to = "method", values_to = "SE") %>%
    pull(SE)
  
  # get quantile estimates per method
  q5s <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("q0.05"), names_to = "method", values_to = "q5") %>%
    pull(q5)
  q10s <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("q0.1"), names_to = "method", values_to = "q10") %>%
    pull(q10)
  q95s <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("q0.95"), names_to = "method", values_to = "q95") %>%
    pull(q95)
  
  # join data
  tidy_results <- ATEhats %>%
    mutate(SE = SEs,
           q5 = q5s,
           q10 = q10s,
           q95 = q95s)
  
  if (F) {
    temp1 <- read_csv("final_sims/example_full.csv")
    temp2 <- read_csv("final_sims/example2_full.csv")
    
    tidy_results %>% 
      bind_rows(temp1 %>% filter(runID <= 400)) %>%
      bind_rows(temp2 %>% filter(runID <= 400)) %>% 
      write_csv(glue("{dir}/{fhead}_txsd.csv"))
  }
  
  write_csv(tidy_results, glue("{dir}/{fhead}_full.csv"))
}