library(ape)
context("Run binary original") #-----------------------------------------------#

test_that("run_binary_original() doesn't give any errors when given a discrete phenotype and snps that group into genes" ,{
  test_dir <- "../../data/outputs/"
  test_pheno <- "../../data/discrete_phenotype.tsv"
  test_tree  <- "../../data/tree.tree"
  test_geno  <- "../../data/grouped_genotype.tsv"
  test_gene_snp_lookup <- "../../data/snp_gene_key.tsv"
  test_annot <- NULL
  test_name  <- "dummy_group_snps_into_genes"
  test_perm  <- "1000"
  test_fdr <- "0.15"
  test_bootstrap <- "0.7"
  test_group_genotype <- TRUE
  args                        <- NULL
  args$test                   <- FALSE
  args$phenotype              <- read_in_tsv_matrix(test_pheno)
  args$tree                   <- read.tree(test_tree)
  args$tree$node.label[1]     <- 0
  args$tree$node.label        <- as.numeric(args$tree$node.label)
  args$genotype               <- read_in_tsv_matrix(test_geno)
  args$output_name            <- test_name # Ex: log_toxin_snp_stop
  args$output_dir             <- test_dir # Directory in which all output files will be saved
  args$perm                   <- as.numeric(test_perm) #typically 10,000
  args$fdr                    <- as.numeric(test_fdr)
  args$annotation             <- NULL #read_in_tsv_matrix(test_annot)
  args$discrete_or_continuous <- check_input_format(args$phenotype, args$tree, args$genotype, args$output_name, args$output_dir, args$perm, args$fdr, args$annot)
  args$bootstrap_cutoff       <- as.numeric(test_bootstrap)
  args$group_genotype         <- test_group_genotype
  if (args$group_genotype ) {
    args$gene_snp_lookup      <- read_in_tsv_matrix(test_gene_snp_lookup)
  } else {
    args$gene_snp_lookup <- NULL
  }
  expect_error(run_binary_original(args), NA)
})

test_that("run_binary_original() doesn't give any errors when given a discrete phenotype and snps (do not group into genes)" ,{
  test_dir <- "../../data/outputs/"
  test_pheno <- "../../data/discrete_phenotype.tsv"
  test_tree  <- "../../data/tree.tree"
  test_geno  <- "../../data/grouped_genotype.tsv"
  test_gene_snp_lookup <- NULL
  test_annot <- NULL
  test_name  <- "dummy_group_snps_into_genes"
  test_perm  <- "1000"
  test_fdr <- "0.15"
  test_bootstrap <- "0.7"
  test_group_genotype <- FALSE
  args                        <- NULL
  args$test                   <- FALSE
  args$phenotype              <- read_in_tsv_matrix(test_pheno)
  args$tree                   <- read.tree(test_tree)
  args$tree$node.label[1]     <- 0
  args$tree$node.label        <- as.numeric(args$tree$node.label)
  args$genotype               <- read_in_tsv_matrix(test_geno)
  args$output_name            <- test_name # Ex: log_toxin_snp_stop
  args$output_dir             <- test_dir # Directory in which all output files will be saved
  args$perm                   <- as.numeric(test_perm) #typically 10,000
  args$fdr                    <- as.numeric(test_fdr)
  args$annotation             <- NULL #read_in_tsv_matrix(test_annot)
  args$discrete_or_continuous <- check_input_format(args$phenotype, args$tree, args$genotype, args$output_name, args$output_dir, args$perm, args$fdr, args$annot)
  args$bootstrap_cutoff       <- as.numeric(test_bootstrap)
  args$group_genotype         <- test_group_genotype
  if (args$group_genotype ) {
    args$gene_snp_lookup      <- read_in_tsv_matrix(test_gene_snp_lookup)
  } else {
    args$gene_snp_lookup <- NULL
  }
  expect_error(run_binary_original(args), NA)
})

