## code to prepare `DATASET` dataset goes here

usethis::use_data("DATASET")
probe_id <-
  read.csv('data-raw/CoRSIV_ESS_SIV_CG_sites_clusters_hg38.csv') %>%
  mutate(probe_id = 1)
usethis::use_data(probe_id)

corsiv_data <-
  read.csv('data-raw/corsiv_table3.csv') %>%
  mutate(corsiv_data = 1)
usethis::use_data(corsiv_data)


kidney_corsiv_data <-
  read.csv('data-raw/Kidney.csv') %>%
  mutate(kidney_corsiv_data = 1)
usethis::use_data(kidney_corsiv_data)

