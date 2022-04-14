# Builds individual phenotype file and generates necessary directories for snakemake workflow
# This file would need to be updated based on the measures that the investigator is taking to adjust the continuous measure
# In this example, IDSA_severity is split into multiple separate data frames for analysis grouped by elix weighted score
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#FILES ----
pheno <- readr::read_csv(file = snakemake@input[["metadata"]])
phenotype_cols <- snakemake@wildcards[["phenotype"]]

paths <- paste0(c("data/pheno/",
                  "results/",
                  "benchmarks/",
                  "figures/",
                  "log/",
                  "data/mikropml/"),
                phenotype_cols)
paths <- c("data/pheno", "data/mikropml", paths)

for(i in 1:length(paths)){
  
  if(dir.exists(paths[i]) == FALSE){
    
    dir.create(paths[i])
    
  }else{
    
    print(paste0("directory ", paths[i], " already exists"))
    
  }
}

new_dir3 <- paste0("results/", phenotype_cols, "/runs")

if(dir.exists(new_dir3) == FALSE){
  
  dir.create(new_dir3)
    
  }else{
    
    print(paste0("directory ", new_dir3, " already exists"))
    
  }

if(!all(sapply(paths, dir.exists))){
  stop("Not all directories created")
}

#UNADJUSTED FILE ----
unadjusted_out <- pheno %>%
  select(genome_id,
         all_of(phenotype_cols))
  
write_tsv(unadjusted_out,
          file = paste0("data/pheno/", phenotype_cols, "/full.tsv"))

#GENERATE SPLIT SEVERITY FILES ----
quartiles <- pheno %>%
    select(elix_weighted_update) %>%
    drop_na() %>%
    deframe() %>%
    quantile()

q1 <- pheno %>%
  filter(elix_weighted_update <= quartiles[2]) %>%
  select(genome_id,
         all_of(phenotype_cols))

q2_3 <- pheno %>%
  filter(elix_weighted_update > quartiles[2] & elix_weighted_update <= quartiles[4]) %>%
  select(genome_id,
         all_of(phenotype_cols))

q4 <- pheno %>%
  filter(elix_weighted_update > quartiles[4]) %>%
  select(genome_id,
         all_of(phenotype_cols))

if(nrow(q1)+nrow(q2_3)+nrow(q4) != nrow(pheno)){
  stop("quartile splits are missing data")
}

write_tsv(q1,
          file = paste0("data/pheno/", phenotype_cols, "/q1.tsv"))
write_tsv(q2_3,
          file = paste0("data/pheno/", phenotype_cols, "/q2_3.tsv"))
write_tsv(q4,
          file = paste0("data/pheno/", phenotype_cols, "/q4.tsv"))
  
even_split <- pheno %>%
    select(elix_weighted_update) %>%
    drop_na() %>%
    deframe() %>%
    quantile(probs = seq(0, 1, 1/3))

t1 <- pheno %>%
  filter(elix_weighted_update <= even_split[2]) %>%
  select(genome_id,
         all_of(phenotype_cols))

t2 <- pheno %>%
  filter(elix_weighted_update > even_split[2] & elix_weighted_update <= even_split[3]) %>%
  select(genome_id,
         all_of(phenotype_cols))

t3 <- pheno %>%
  filter(elix_weighted_update > even_split[3]) %>%
  select(genome_id,
         all_of(phenotype_cols))

if(nrow(t1)+nrow(t2)+nrow(t3) != nrow(pheno)){
  stop("tertile splits are missing data")
}

write_tsv(t1,
          file = paste0("data/pheno/", phenotype_cols, "/t1.tsv"))
write_tsv(t2,
          file = paste0("data/pheno/", phenotype_cols, "/t2.tsv"))
write_tsv(t3,
          file = paste0("data/pheno/", phenotype_cols, "/t3.tsv"))
