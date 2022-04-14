source("code/log_smk.R")
library(tidyverse)
library(mikropml)

doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])
options(future.globals.maxSize= +Inf)

data_raw <- readr::read_csv(snakemake@input[["csv"]])

data_raw <- data_raw %>%
  mutate(across(snakemake@params[['outcome_colname']], ~as.character(.x))) #only necessary when binary measure and 0/1 being used

dim(data_raw)
data_processed <- preprocess_data(data_raw,
                                  outcome_colname = snakemake@params[['outcome_colname']],
                                  group_neg_corr = TRUE,
                                  remove_var = 'nzv')

summary(data_processed)
full_data_processed <- data_processed

full <- readr::read_delim(file = snakemake@input[["full"]],
                          delim = "\t")

sub <- readr::read_delim(file = paste0("data/pheno/",
                                       snakemake@params[['outcome_colname']],
                                       "/",
                                       snakemake@wildcards[["group"]],
                                       ".tsv"),
                         delim = "\t")

pos_select <- sapply(sub$genome_id, function(x){grep(x, full$genome_id)})

if(length(pos_select) != nrow(sub)){
  stop("selection vector is incorrect length")
}

data_processed$dat_transformed <- data_processed$dat_transformed[pos_select,,drop=FALSE]

if(dim(data_processed$dat_transformed)[1] != length(pos_select)){
  stop("selection of transformed data is incorrect number of rows")
}

if(dim(data_processed$dat_transformed)[2] != dim(full_data_processed$dat_transformed)[2]){
  stop("selection of transformed data is incorrect number of columns")
}

summary(data_processed)
dim(data_processed$dat_transformed)

saveRDS(data_processed, file = snakemake@output[["rds"]])
