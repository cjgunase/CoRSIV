#data-raw/process.R
# Data import and processing pipeline

library(readr)
library(readxl)

probe_id <- read_csv("data-raw/CoRSIV_ESS_SIV_CG_sites_clusters_hg38.csv")
corvis_data <- read_csv("data-raw/corsiv_table3.csv")
kidney_probe_data <- read_csv("data-raw/Kidney.csv")

# Data processing code here...

# This should be the last line
devtools::use_data(probe_id, corvis_data, kidney_probe_data,overwrite = T)

