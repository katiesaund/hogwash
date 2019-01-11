# Libraries --------------------------------------------------------------------
source("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/lib/2018-09-06_plot_gwas_hits_lib.R")

get_gee_hits <- function(path, pheno, geno){
  file_name <- paste(path, "gee_", pheno, geno, "_hits_genotype_names.tsv", sep = "")
  hits <- read_in_tsv_matrix(file_name)
  names <-  as.vector(unlist(hits))
  return(names)
} # end get_gee_hits()

get_phyc_hits <- function(path, pheno, geno, type_of_phyc){
  if (pheno %in% c("fqR", "severity")){ # binary phenotypes
    file_name <- paste(path, "phyc_", pheno, geno, "_", type_of_phyc, "_pvals_for_sig_hits.tsv", sep = "")
    if (file.exists(file_name)){
      hits <- read_in_tsv_matrix(file_name)
      names <- as.vector(row.names(hits))
    } else {
      names <- "no_sig_hits"
    }
    
    
    # This commented out section below is for whenI want to combine all reconstruction and transition results together. 
    # For the LAMGC & U01 meeting I'm only planning to look at transition phyC results (so discrete and continuous results are more comparable)
    # file_name_recon <- paste(path, "phyc_", pheno, geno, "_reconstruction_pvals_for_sig_hits.tsv", sep = "")
    # file_name_trans <- paste(path, "phyc_", pheno, geno, "_transition_pvals_for_sig_hits.tsv", sep = "")
    # if (file.exists(file_name_recon)){
    #   hits <- read_in_tsv_matrix(file_name_recon)
    #   recon_temp <- as.vector(row.names(hits))
    #   
    # } else {
    #   recon_temp <- NULL
    # }
    # 
    # if (file.exists(file_name_trans)){
    #   hits <- read_in_tsv_matrix(file_name_trans)
    #   trans_temp <- as.vector(row.names(hits))
    # } else {
    #   trans_temp <- NULL
    # }
    # 
    # if (is.null(trans_temp) & is.null(recon_temp)){
    #   names <- as.vector("no_sig_hits")
    # } else if (is.null(trans_temp) & !is.null(recon_temp)){
    #   names <- recon_temp
    # } else if (is.null(recon_temp) & !is.null(trans_temp)){
    #   names <- trans_temp
    # } else {
    #   names <- c(trans_temp, recon_temp)  
    # }
    
    
  } else { #log_toxin, log_growth, etc... Continuous phenotypes
    file_name <- paste(path, "phyc_", pheno, geno, "_all_transitions_pvals_for_sig_hits.tsv", sep = "")
    if (file.exists(file_name)){
      hits <- read_in_tsv_matrix(file_name)
      names <- as.vector(row.names(hits))
    } else {
      names <- as.vector("no_sig_hits")
    }
  }
  
  return(names)
} # end get_phyc_hits()

get_treewas_hits <- function(path, pheno, geno, type_of_hit){
  file_name <- paste(path, "treewas_", pheno, geno, "_", type_of_hit, "_sig_hits.tsv", sep = "")
  hits <- read_in_tsv_matrix(file_name)
  names <-  row.names(hits)
  return(names)
} # end get_treewas_hits()


# ------------------------------------------------------------------------------

phenotypes <- c("log_cfe", "log_sporulation", "log_germ_tc", "log_germ_tc_and_gly", "log_growth",  "log_toxin", "fqR", "severity") 
genotypes  <- c("_gene_ns", "_gene_del", "_gene_high", "_gene_stop", "_snp_del", "_snp_high", "_snp_stop", "_roary_pan_genome", "_pilon_sv") #left out _snp_ns because as of 2018-09-06 9pm snp_ns is still not done running for wany GWAS method

# change path to getwd()
gwas_intersection_results <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/compare_methods/all_methods_intersection_per_test_sig_hits.rda"
gwas_two_method_results <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/compare_methods/two_plus_methods_intersection_per_test_sig_hits.rda"
gwas_any_method_results <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/compare_methods/any_per_test_sig_hits.rda"
path <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/"


# CHECK INPUTS -----------------------------------------------------------------
check_file_exists(gwas_intersection_results)
check_file_exists(gwas_two_method_results)
check_file_exists(gwas_any_method_results)

load(gwas_intersection_results)
load(gwas_two_method_results)
load(gwas_any_method_results)

# FUNCTION ---------------------------------------------------------------------
summary_count <- matrix(0, nrow = 3, ncol = length(genotypes) * length(phenotypes))
row.names(summary_count) <- c("any", "two_plus", "all_three")
colnames(summary_count) <- 1:ncol(summary_count)

count_per_method <- matrix(0, nrow = 7, ncol = length(genotypes) * length(phenotypes))
row.names(count_per_method) <- c("gee", "phyc_trans", "phyc_recon", "treewas_pooled", "treewas_simultaneous", "treewas_subsequent", "treewas_terminal")
colnames(count_per_method) <- 1:ncol(count_per_method)


counter <- 0

for (g in 1:length(genotypes)){
  for (p in 1:length(phenotypes)){
    counter <- counter + 1
    name <- paste(phenotypes[p], genotypes[g], sep = "")
    colnames(summary_count)[counter] <- name
    colnames(count_per_method)[counter] <- name
    
    for (i in 1:length(names(intersection_all_methods))){
      if(names(intersection_all_methods)[i] == name){
        summary_count[3, counter] <- length(unlist(intersection_all_methods[i]))
      }
    }
    for (i in 1:length(names(intersection_two_plus_methods))){
      if(names(intersection_two_plus_methods)[i] == name){
        summary_count[2, counter] <- length(unlist(intersection_two_plus_methods[i]))
      }
    }
    for (i in 1:length(names(any_hits))){
      if(names(any_hits)[i] == name){
        summary_count[1, counter] <- length(unlist(any_hits[i]))
      }
    }
    
    gee          <- get_gee_hits(path, phenotypes[p], genotypes[g])
    phyc_trans   <- get_phyc_hits(path, phenotypes[p], genotypes[g], "transition")
    phyc_recon   <- get_phyc_hits(path, phenotypes[p], genotypes[g], "reconstruction")
    treewas_comb <- get_treewas_hits(path, phenotypes[p], genotypes[g], "combined")
    treewas_sim  <- get_treewas_hits(path, phenotypes[p], genotypes[g], "simultaneous")
    treewas_sub  <- get_treewas_hits(path, phenotypes[p], genotypes[g], "subsequent")
    treewas_term <- get_treewas_hits(path, phenotypes[p], genotypes[g], "terminal")
    
    
    count_per_method[1, counter] <- length(gee)
    count_per_method[2, counter] <- length(phyc_trans)
    count_per_method[3, counter] <- length(phyc_recon)
    
    if (length(treewas_comb) == 1){
      if (treewas_comb == "1"){ # catch when no sign hits found
        count_per_method[4, counter] <- 0
      }
    } else {
      count_per_method[4, counter] <- length(treewas_comb)
    }
    
    if (length(treewas_sim) == 1){
      if (treewas_sim == "1"){ # catch when no sign hits found
        count_per_method[5, counter] <- 0
      }
    } else {
      count_per_method[5, counter] <- length(treewas_sim)
    }
    
    if (length(treewas_sub) == 1){
      if (treewas_sub == "1"){ # catch when no sign hits found
        count_per_method[6, counter] <- 0
      }
    } else {
      count_per_method[6, counter] <- length(treewas_sub)
    }
    
    if (length(treewas_term) == 1){
      if (treewas_term == "1"){ # catch when no sign hits found
        count_per_method[7, counter] <- 0
      }
    } else {
      count_per_method[7, counter] <- length(treewas_term)
    }
    
    
  }
}

print(summary_count)
print(count_per_method)


white_to_black = colorRampPalette(c("white", "black")) # Set up presence/absence colors
num_color = length(min(summary_count):max(summary_count)) # We want a unique shade of grayscale for each level (should be 2 for presence/absence)
htmp_col = white_to_black(num_color) # Generate a character vector of colors with each shade of grayscale
htmp <- ComplexHeatmap::Heatmap(
  matrix = summary_count, 
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  row_names_side = "left", 
  col = htmp_col, 
  show_column_names = TRUE, 
  show_heatmap_legend = TRUE, 
  row_title = "GWAS methods",
  column_title = "Number of significant hits", column_names_gp = gpar(cex = 0.6))

draw(htmp)
pdf("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/summary_hits_per_method_heatmap.pdf")
draw(htmp)
dev.off()

write.table(summary_count, file = "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/summary_hits_per_method.tsv", sep = "\t")



white_to_black = colorRampPalette(c("white", "black")) # Set up presence/absence colors
num_color = length(min(count_per_method):max(count_per_method)) # We want a unique shade of grayscale for each level (should be 2 for presence/absence)
htmp_col = white_to_black(num_color) # Generate a character vector of colors with each shade of grayscale
htmp <- ComplexHeatmap::Heatmap(
  matrix = count_per_method, 
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  row_names_side = "left", 
  col = htmp_col, 
  show_column_names = TRUE, 
  show_heatmap_legend = TRUE, 
  row_title = "GWAS methods", 
  column_title = "Number of significant hits", column_names_gp = gpar(cex = 0.6))

draw(htmp)
pdf("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/hits_per_method_heatmap.pdf")
draw(htmp)
dev.off()

write.table(count_per_method, file = "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/hits_per_method.tsv", sep = "\t")


# END OF SCRIPT ----------------------------------------------------------------
