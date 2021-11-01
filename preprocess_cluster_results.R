
# get cluster results down into a single file

require(tidyverse)
require(glue)

dir <- "results2"
fhead <- "sim_results_edited"

all_files <- list.files(dir)

max_runID <- 0
all_dfs <- all_files[str_detect(all_files, glue("{fhead}-"))] %>%
  map(function(fname) {
    df <- read_csv(glue("{dir}/{fname}")) %>%
      mutate(runID = runID + max_runID)
    max_runID <<- max(df$runID)
    df
  })

max_runID <- 0
all_dfs_overall <- all_files[str_detect(all_files, glue("{fhead}_overall"))] %>%
  map(function(fname) {
    df <- read_csv(glue("{dir}/{fname}")) %>%
      mutate(runID = runID + max_runID)
    max_runID <<- max(df$runID)
    df
  })

res <- bind_rows(all_dfs)
res_overall <- bind_rows(all_dfs_overall)


rm(all_files)
rm(all_dfs)
rm(all_dfs_overall)
