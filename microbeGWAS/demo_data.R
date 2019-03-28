library(microbeGWAS)
library(ape)       # ape::ace function (ancestral reconstruction)
library(phytools)  # phylogenetic tree function library
library(grid) # plots for continuous phenotypes
library(gridExtra) # plots for continuous phenotypes
library(phangorn)
library(pheatmap) # plots for continuous phenotypesq
library(ggplot2) # plots for continuous phenotypes

test_dir <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-03-21_figures_for_GTP_talk/data/"
# test_dir   <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-03-18_gwas_fix_geno_to_be_trans/data/"
dataset <- 0 # Discrete, gene test built from SNPS
# dataset <- 1 # Discrete, gene test not built from SNPs
#dataset <- 2 # Continuous, gene test built from SNPS
#dataset <- 3 # Continuous, gene test not built from SNPs
#dataset <- 4 # Continuous, STOP SNP

if (dataset == 0){
  # Discrete, gene test built from SNPS
  test_pheno <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/fqR_pheno.tsv"
  test_tree  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/fqR.tree"
  test_geno  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/fqR_snp_stop.tsv"
  test_annot <- NULL
  test_name  <- "fqR_gene_built_from_stop_snps"
  test_perm  <- "1000"
  test_alpha <- "0.15"
  test_bootstrap <- "0.7"
  test_gene_snp_lookup <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-03-05_stop_snp_gene_lookup.tsv"
  test_built_from_snps <- TRUE
} else if (dataset == 1){
  # Discrete, gene test not built from SNPs
  test_pheno <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/fqR_pheno.tsv"
  test_tree  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/fqR.tree"
  test_geno  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/fqR_gene_stop.tsv"
  test_annot <- NULL
  test_name  <- "fqR_gene_stop"
  test_perm  <- "1000"
  test_alpha <- "0.15"
  test_bootstrap <- "0.7"
  test_gene_snp_lookup <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-03-05_stop_snp_gene_lookup.tsv"
  test_built_from_snps <- FALSE
} else if (dataset == 2){
  # Continuous, gene test built from SNPS
  test_pheno <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/log_toxin_pheno.tsv"
  test_tree  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/log_toxin.tree"
  test_geno  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/log_toxin_snp_stop.tsv"
  test_annot <- NULL
  test_name  <- "log_toxin_gene_built_from_stop_snps"
  test_perm  <- "1000"
  test_alpha <- "0.15"
  test_bootstrap <- "0.7"
  test_gene_snp_lookup <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-03-05_stop_snp_gene_lookup.tsv"
  test_built_from_snps <- TRUE
} else if (dataset == 3){
  #Continuous, gene test not built from SNPs
  test_pheno <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/log_toxin_pheno.tsv"
  test_tree  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/log_toxin.tree"
  test_geno  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/log_toxin_gene_stop.tsv"
  test_annot <- NULL
  test_name  <- "log_toxin_gene_stop"
  test_perm  <- "1000"
  test_alpha <- "0.15"
  test_bootstrap <- "0.7"
  test_gene_snp_lookup <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-03-05_stop_snp_gene_lookup.tsv"
  test_built_from_snps <- FALSE
} else if (dataset == 4){
  #Continuous, STOP SNPs
  test_pheno <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/log_toxin_pheno.tsv"
  test_tree  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/log_toxin.tree"
  test_geno  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/log_toxin_snp_stop.tsv"
  test_annot <- NULL
  test_name  <- "log_toxin_snp_stop"
  test_perm  <- "10000"
  test_alpha <- "0.05"
  test_bootstrap <- "0.7"
  test_gene_snp_lookup <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-03-05_stop_snp_gene_lookup.tsv" # TODO update this to be null
  test_built_from_snps <- FALSE
}

args                        <- NULL
args$test                   <- FALSE
args$phenotype              <- read_in_tsv_matrix(test_pheno)
args$tree                   <- read.tree(test_tree)
args$genotype               <- read_in_tsv_matrix(test_geno)
args$output_name            <- test_name # Ex: log_toxin_snp_stop
args$output_dir             <- test_dir # Directory in which all output files will be saved
args$perm                   <- as.numeric(test_perm) #typically 10,000
args$alpha                  <- as.numeric(test_alpha)
args$annotation             <- NULL #read_in_tsv_matrix(test_annot)
args$discrete_or_continuous <- check_input_format(args$phenotype, args$tree, args$genotype, args$output_name, args$output_dir, args$perm, args$alpha, args$annot)
args$bootstrap_cutoff       <- test_bootstrap
args$gene_snp_lookup        <- read_in_tsv_matrix(test_gene_snp_lookup)
args$built_from_snps        <- test_built_from_snps

run_phyc(args)
