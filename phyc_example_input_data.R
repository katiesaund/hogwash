# 2018-12-14
# Katie Saund
# Paths to correctly formatted data

data_type <- "continuous" # or "discrete

if (data_type == "continuous"){
  # continuous
  test_pheno <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/log_toxin_pheno.tsv" 
  test_tree  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/log_toxin.tree"
  test_geno  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/log_toxin_gene_stop.tsv"
  test_annot <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/log_toxin_annotation.tsv"
  test_dir   <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-12-14_phyc_fix_delta_pheno/data/2018-12-18_morning/"
  test_name  <- "log_toxin_gene_stop"
  test_perm  <- "1000"
  test_alpha <- "0.1"
  args       <- NULL
  args$test                   <- FALSE
  args$phenotype              <- read_in_tsv_matrix(test_pheno)
  args$tree                   <- read.tree(test_tree)
  args$genotype               <- read_in_tsv_matrix(test_geno)
  args$output_name            <- test_name # Ex: log_toxin_snp_stop
  args$output_dir             <- test_dir # Directory in which all output files will be saved
  args$perm                   <- as.numeric(test_perm) #typically 10,000
  args$alpha                  <- as.numeric(test_alpha)
  args$annotation             <- read_in_tsv_matrix(test_annot)
  args$discrete_or_continuous <- check_input_format(args$phenotype, args$tree, args$genotype, args$output_name, args$output_dir, args$perm, args$alpha, args$annot)
} else if (data_type == "discrete"){ 

  # Discrete
  test_pheno <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/fqR_treewas_pheno.tsv"
  test_tree  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/fqR.tree"
  test_geno  <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/fqR_gene_stop.tsv"
  test_annot <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/fqR_annotation.tsv"
  test_dir   <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-12-14_phyc_fix_delta_pheno/data/2018-12-17_evening_results/"
  test_name  <- "fqR_gene_stop"
  test_perm  <- "1000"
  test_alpha <- "0.05"
  args       <- NULL
  args$test                   <- FALSE
  args$phenotype              <- read_in_tsv_matrix(test_pheno)
  args$tree                   <- read.tree(test_tree)
  args$genotype               <- read_in_tsv_matrix(test_geno)
  args$output_name            <- test_name # Ex: log_toxin_snp_stop
  args$output_dir             <- test_dir # Directory in which all output files will be saved
  args$perm                   <- as.numeric(test_perm) #typically 10,000
  args$alpha                  <- as.numeric(test_alpha)
  args$annotation             <- read_in_tsv_matrix(test_annot)
  args$discrete_or_continuous <- check_input_format(args$phenotype, args$tree, args$genotype, args$output_name, args$output_dir, args$perm, args$alpha, args$annot)
} else {
  print("data_type should be discrete or continuous")
}



