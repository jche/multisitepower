
# get cluster results down into a single file

require(tidyverse)
require(glue)

dir <- "results_sree"
# dir <- "results2"   # results2 is ONLY J=20 cases
fhead <- "sree_sims"

all_files <- list.files(dir)

max_runID <- 0
all_dfs <- all_files[str_detect(all_files, glue("{fhead}-"))] %>%
  map(function(fname) {
    df <- read_csv(glue("{dir}/{fname}")) %>%
      mutate(runID = runID + max_runID)
    max_runID <<- max(df$runID)
    df
  })
res <- bind_rows(all_dfs)

if (F) {
  max_runID <- 0
  all_dfs_overall <- all_files[str_detect(all_files, glue("{fhead}_overall"))] %>%
    map(function(fname) {
      df <- read_csv(glue("{dir}/{fname}")) %>%
        mutate(runID = runID + max_runID)
      max_runID <<- max(df$runID)
      df
    })
  res_overall <- bind_rows(all_dfs_overall)
}



rm(all_files)
rm(all_dfs)
rm(all_dfs_overall)



if (F) {
  hits <- res
  
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
  write_csv(tidy_results, "results_sree/sree_sims.csv")
}