# Create metaboanlyst mummichog-ready tables from SAP19_MET.csv using zea mays KEGG

library(magrittr)
library(data.table)
library(purrr)
library(furrr)
library(future)
library(MetaboAnalystR)
library(RJSONIO)
library(fitdistrplus)
library(parallel)
library(tidyverse)
library(memoise)

# Set working directory to your path
setwd("SAP_Metabolomics")

# Upload csv
met <- fread("./data/raw/untreated_CU_met/SAP19_MET.csv")

# Set up for parallel computation
plan(multisession, workers = (detectCores()-2)) # Use half of available cores 
furrr_options(scheduling = 2L) # dynamic distribution

# Format csv for mummichog analysis
# Create a list of categories to be evaluated with entries greater than 5
# change the group_by catagory to change the target catagories
sum <- met %>% group_by(RACE) %>% summarize(n()) %>%
  rename(ncount = 2) %>% 
  filter(ncount >= 5 & !is.na(RACE))

# Convert to a list
jobs <- as.list(sum$RACE)

# Replicate the met df for as many catagories as we have, lists must be equal to be evaluated
met_list <- replicate(length(jobs), met, simplify = FALSE)

# Create a function to calculate p values and perform mummichog analysis
met_ttest <- function(peak_list, job){
  
  df <- tibble()
  df2 <- peak_list[, GROUP := "average"][job, GROUP := job, on = "RACE"]
  
  df3 <- as_tibble(df2) %>%
    relocate(GROUP, .after = "SAMPLE")
  
  write_csv(df3, str_c("peak_list_", job, ".csv"))
  
  mSet<-InitDataObjects("pktable", "stat", FALSE)
  mSet<-Read.TextData(mSet, str_c("peak_list_", job, ".csv"), "rowu", "disc");
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet);
  mSet<-SanityCheckData(mSet)
  mSet<-FilterVariable(mSet, "none", -1, "F", 25, F)
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "QuantileNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
  mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
  mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
  mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "raw", FALSE)
  
  str_split_fixed(mSet[["dataSet"]][["prenorm.feat.nms"]], "/", n = 2)
  
  df <-  data.table(str_split_fixed(mSet[["dataSet"]][["prenorm.feat.nms"]], "/", n = 2), mSet[["analSet"]][["tt"]][["t.score"]],
                    mSet[["analSet"]][["tt"]][["p.value"]]) %>% 
    rename("m.z" = 1, "ret" = 2, "t.score" = 3, "p.value" = 4) %>% 
    arrange(p.value)
  
  write_csv(df, str_c("t_test_", job, ".csv"))
  
  rm(mSet)
  
  # Perform mummichog analysis
  mSet<-InitDataObjects("mass_all", "mummichog", FALSE)
  mSet<-SetPeakFormat(mSet, "rmp")
  mSet<-UpdateInstrumentParameters(mSet, 5.0, "negative", "yes", 0.02);
  mSet<-Read.PeakListData(mSet, str_c("t_test_", job, ".csv"));
  mSet<-SanityCheckMummichogData(mSet)
  mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
  mSet<-SetMummichogPval(mSet, 0.05)
  mSet<-PerformPSEA(mSet, "zma_kegg", "current", 3 , 100)
  df_res <- as_tibble(mSet[["path.nms"]]) %>% 
    bind_cols(as_tibble(mSet[["mummi.resmat"]])) %>% 
    mutate(RACE = job)
  return(df_res)
}

# Wrap function using safely so it will proceed if a pvalue results in NA
met_ttest_safe <- safely(met_ttest, quiet = TRUE)

# Use purrr to perform function across the lists that we created 
peak_lists <- map2(met_list, jobs, ~met_ttest_safe(.x, .y), .progress = TRUE)

# Save as R data
saveRDS(peak_lists$results, "peak_lists.RData")

# Write CSV
final_peak_list <- peak_lists %>% future_map("result") %>% compact() %>% rbindlist()
write_csv(final_peak_list, "final_peak_list_RACE.csv")
