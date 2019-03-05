library(microbeGWAS)
library(ape)       # ape::ace function (ancestral reconstruction)
library(phytools)  # phylogenetic tree function library
library(ComplexHeatmap) # to make final plots for discrete phenotypes
library(phangorn)
library(pheatmap) # plots for continuous phenotypes
library(grid) # plots for continuous phenotypes
library(gridExtra) # plots for continuous phenotypes
library(ggplot2) # plots for continuous phenotypes

test_pheno <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/fqR_pheno.tsv"
test_tree  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/fqR.tree"
test_geno  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-04_format_data_for_gwas/data/2019-02-04_gwas_formatted_data/fqR_snp_stop.tsv"
test_annot <- NULL
test_dir   <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-02-25_simplify_phyc_output/data"
test_name  <- "fqR_snp_stop"
test_perm  <- "1000"
test_alpha <- "0.3"
test_bootstrap <- "0.7"
args       <- NULL
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


run_phyc(args)
