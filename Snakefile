configfile: 'config/config.yml'

phenotype = config['phenotype']
tree = config['tree']
genome = config['genome']

ncores = config['ncores']
ml_methods = config['ml_methods']
kfold = config['kfold']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

# attempt limit is set to 5 in the file code/submit_slurm.sh under restart-times
def get_mem_mb_low(wildcards, attempt):
     mem = attempt*4
     return "%dGB" % (mem)

def get_mem_mb_med(wildcards, attempt):
     mem = attempt*8
     return "%dGB" % (mem)

def get_mem_mb_high(wildcards, attempt):
     mem = attempt*16
     return "%dGB" % (mem)

rule all:
    input:
        expand("aggregated/{phenotype}.runs.csv", phenotype=phenotype)

# This checkpoint is for when the continuous phenotype of interest may have a raw value and a value that is adjusted for a covariate
# The script prepro_overall.R needs to be edited by the user for their phenotype of interest and its relevant covariate(s)
# If there are no covariate(s) for adjusting the phenotype, within the prepro_overall.R script the adjustment workflow can be commented out
# The group wildcard will define the different phenotype files generated during this checkpoint, which can be the single raw file, or multiple covariate files depending on user input in the prepro_overall.R script
# Once this checkpoint is complete, snakemake will re-evaluate the jobs that are required to complete the necessary downstream file creation
checkpoint prepro_overall:
    input:
        R = "code/prepro_overall.R",
        metadata = config['metadata']
    output:
        dat_dir = directory("data/pheno/{phenotype}")
    log:
        "log/{phenotype}/prepro_overall.txt"
    resources:
        ncores = ncores,
        mem_mb = get_mem_mb_low
    script:
        "code/prepro_overall.R"

# This rule generates the phenotype:genotype data frame in the necessary input format for mikropml's preprocessing function (see rule in mikropml.smk, preprocess_data)
rule generate_mikropml_df:
    input:
        R = "code/generate_mikropml_df.R",
        pheno = "data/pheno/{phenotype}/{group}.tsv"
    output:
        file_name = "data/mikropml/{phenotype}/{group}.{genome}.csv"
    params:
        core_path = config['core'],
        pan_path = config['pan'],
        gene_path = config['gene']
    log:
        "log/{phenotype}/{group}.{genome}.generate_mikropml_df.txt"
    resources:
        ncores = ncores,
        mem_mb = get_mem_mb_med
    script:
        "code/generate_mikropml_df.R"

include: "mikropml.smk"

rule run_treeWAS:
    input:
        R = "code/run_treeWAS.R",
        pheno = "data/pheno/{phenotype}/{group}.tsv",
        rds = rules.preprocess_data.output.rds
    output:
        rdata = 'results/{phenotype}/{group}.{genome}.treeWAS.RData',
        plot = 'results/{phenotype}/{group}.{genome}.treeWAS.pdf'
    params:
        tree = tree
    log:
        "log/{phenotype}/{group}.{genome}.treeWAS.txt"
    resources:
        ncores = ncores,
        mem_mb = get_mem_mb_high
    script:
        "code/run_treeWAS.R"

rule run_hogwash_ungrouped:
    input:
        R = "code/run_hogwash_ungrouped.R",
        pheno = "data/pheno/{phenotype}/{group}.tsv",
        rds = rules.preprocess_data.output.rds
    output:
        rdata = "results/{phenotype}/hogwash_continuous_{group}.{genome}.ungrouped.rda",
        plot = "results/{phenotype}/hogwash_continuous_{group}.{genome}.ungrouped.pdf"
    params:
        tree = tree,
        file_name = '{group}.{genome}.ungrouped',
        dir = "results/{phenotype}"
    log:
        "log/{phenotype}/{group}.{genome}.hogwash.ungrouped.txt"
    resources:
        ncores = ncores
    script:
        "code/run_hogwash_ungrouped.R"

# hogwash grouped analysis and relevant functions are commmented out until optimization of the function is complete
# rule run_hogwash_grouped:
#     input:
#         R = "code/run_hogwash_grouped.R",
#         pheno = "data/pheno/{phenotype}/{group}.tsv",
#         rds = rules.preprocess_data.output.rds
#     output:
#         rdata = protected("results/{phenotype}/hogwash_continuous_{group}.{genome}.grouped.rda"),
#         plot = protected("results/{phenotype}/hogwash_continuous_{group}.{genome}.grouped.pdf")
#     params:
#         tree = tree,
#         file_name = '{group}.{genome}.grouped',
#         dir = "results/{phenotype}",
#         gene_key = config['gene_key']
#     wildcard_constraints:
#         genome = "core"
#     log:
#         "log/{phenotype}/{group}.{genome}.hogwash.grouped.txt"
#     resources:
#         ncores = ncores
#     script:
#         "code/run_hogwash_grouped.R"

def aggregate_input1(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{phenotype}/{group}.{genome}.treeWAS.RData',
        phenotype=wildcards.phenotype,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
        genome=genome)

def aggregate_input2(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{phenotype}/hogwash_continuous_{group}.{genome}.ungrouped.rda',
        phenotype=wildcards.phenotype,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
        genome=genome)

# def aggregate_input3(wildcards):
#     checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
#     return expand('results/{phenotype}/hogwash_continuous_{group}.{genome}.grouped.rda',
#         phenotype=wildcards.phenotype,
#         group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
#         genome = "core")

def aggregate_input4(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{phenotype}/{group}.{genome}.report.md',
        phenotype=wildcards.phenotype,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
        genome=genome)

# if 'core' in genome:
#     finish_list = [aggregate_input1, aggregate_input2, aggregate_input3, aggregate_input4]
# else:
    # finish_list = [aggregate_input1, aggregate_input2, aggregate_input4]

finish_list = [aggregate_input1, aggregate_input2, aggregate_input4]

rule finish_test:
    input:
        finish_list
    output:
        "aggregated/{phenotype}.runs.csv"
    log:
        "log/{phenotype}/finish.txt"
    script:
        "code/assemble_files.py"
