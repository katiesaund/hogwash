# Katie Saund
# 2018-09-06

# Script will gather all significant hits from running GEE, treeWAS, and phyC on 
# a set of phenotypes and genotypes. 


# Intersection means hits common to the datasets (can be present in others).
# Unique means hits found exclusively in the single dataset (cannot be present in others).
# All means any hit in a single dataset (can be present in others).

# lib
check_is_string <- function(char){
  # Function description ------------------------------------------------------- 
  # Check that the input is a character string. 
  #
  # Input: 
  # char. Character.
  #
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (!is.character(char)){
    stop("Object must be a character string.")
  }
} # end check_is_string()


check_file_exists <- function(file_name){
  # Function description ------------------------------------------------------- 
  # Check that the file exists. 
  #
  # Input: 
  # file_name. Character.   
  #
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (!file.exists(file_name)){
    stop("File does not exist")
  }
} # end check_file_exists()

read_in_tsv_matrix <- function(mat){
  # Function description ------------------------------------------------------- 
  # Read in the standardized GWAS matrix format: rows correspond to samples, columns correspond to genotypes/phenotypes. 
  # Data are tab-separated. 
  #
  # Inputs: 
  # mat. Character. Path to matrix file.  
  # 
  # Output: 
  # temp. Matrix. 
  # 
  # Check inputs ---------------------------------------------------------------
  check_is_string(mat)
  check_file_exists(mat)
  
  # Function -------------------------------------------------------------------
  temp <- read.table(mat,
                     sep = "\t",
                     row.names = 1, 
                     header = TRUE,
                     stringsAsFactors = FALSE)
  temp <- as.matrix(temp)
  
  # Check and return output ----------------------------------------------------
  if (class(temp) != "matrix"){stop("Ouput is incorrectly formatted")}
  return(temp)
} # end read_in_tsv_matrix()

find_intersection <- function(list1, list2){
  intersection <- intersect(list1, list2)
  return(intersection)
} # end find_intersection()

get_gee_hits <- function(path, pheno, geno){
  file_name <- paste(path, "gee_", pheno, geno, "_hits_genotype_names.tsv", sep = "")
  hits <- read_in_tsv_matrix(file_name)
  names <-  as.vector(unlist(hits))
  return(names)
} # end get_gee_hits()

get_phyc_hits <- function(path, pheno, geno){
  if (pheno %in% c("fqR", "severity")){ # binary phenotypes
    file_name_trans <- paste(path, "phyc_", pheno, geno, "_transition_pvals_for_sig_hits.tsv", sep = "")
    if (file.exists(file_name_trans)){
      hits <- read_in_tsv_matrix(file_name_trans)
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

get_treewas_hits <- function(path, pheno, geno){
  file_name <- paste(path, "treewas_", pheno, geno, "_combined_sig_hits.tsv", sep = "")
  hits <- read_in_tsv_matrix(file_name)
  names <-  as.vector(unlist(hits))
  return(names)
} # end get_treewas_hits()

get_unique_hits <- function(method1, method2, method3){
  unique <- method1[!(method1 %in% c(method2, method3))]
  return(unique)
} # end get_unique_hits

save_results <- function(results, name){
  if (length(results) > 0){
    write.table(x = results, file = name, sep = "\t", quote = FALSE, eol = "\n", row.names = TRUE, col.names = TRUE)
  }
} # end save_results()

save_list <- function(results, name){
  if (length(results) > 0){
    con <- file(name, "w")
    writeLines(sapply(names(results),function(x) paste(x,paste(results[[x]],collapse=" "))), con = con)
    close(con)
  }
} # end save_list()

append_findings <- function(lst, title, addition){
  if (length(addition) > 0){
    lst[[title]] <- addition
  }
  return(lst)
} # end append_findings()

# FUNCTION ---------------------------------------------------------------------

phenotypes <- c("log_cfe", "log_germ_tc_and_gly", "log_growth", "log_toxin", "log_sporulation", "log_germ_tc", "fqR", "severity")
genotypes  <- c("_gene_ns", "_gene_stop", "_gene_high", "_gene_del", "_snp_stop","_snp_high", "_snp_del", "_roary_pan_genome", "_pilon_sv") # left out snp_ns which isn't done running as of 2018-09-06 1pm
# change path to getwd()
path <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/"
intersection_all_methods <- intersection_two_plus_methods <- any_hits <- list()

if (!dir.exists(paste(path, "/compare_methods/", sep = ""))){
  system(paste("mkdir ", paste(path, "/compare_methods/", sep = ""), sep = ""))
}

main_dir <- paste(path, "/compare_methods/", sep = "")

for (p in 1:length(phenotypes)){
  for (g in 1:length(genotypes)){
    print("loop start")
    print(paste(phenotypes[p], genotypes[g], sep = ""))
    subdir <- paste(path, "/compare_methods/", paste(phenotypes[p], genotypes[g], "/", sep = ""), sep = "")
    system(paste("mkdir ", subdir, sep = ""))
    
    # only reads in significant hits
    gee     <- get_gee_hits(path, phenotypes[p], genotypes[g])
    phyc    <- get_phyc_hits(path, phenotypes[p], genotypes[g])
    treewas <- get_treewas_hits(path, phenotypes[p], genotypes[g])
    
    gee_unique     <- get_unique_hits(gee, phyc, treewas)
    phyc_unique    <- get_unique_hits(phyc, gee, treewas)
    treewas_unique <- get_unique_hits(treewas, gee, phyc)
    
    gee_phyc_treewas_intersection <- find_intersection(gee, find_intersection(phyc, treewas))
    gee_phyc_intersection         <- find_intersection(gee, phyc)
    gee_treewas_intersection      <- find_intersection(gee, treewas)
    phyc_treewas_intersection     <- find_intersection(phyc, treewas)
    
    # Save results
    save_results(gee, paste(subdir, "/gee_all_sig_hits.tsv", sep = ""))
    save_results(gee_unique, paste(subdir, "/gee_unique_sig_hits.tsv", sep = ""))
    save_results(phyc, paste(subdir, "/phyc_all_sig_hits.tsv", sep = ""))
    save_results(phyc_unique, paste(subdir, "/phyc_unique_sig_hits.tsv", sep = ""))
    save_results(treewas, paste(subdir, "/treewas_all_sig_hits.tsv", sep = ""))
    save_results(treewas_unique, paste(subdir, "/treewas_unique_sig_hits.tsv", sep = ""))

    save_results(gee_phyc_treewas_intersection, paste(subdir, phenotypes[p], genotypes[g],"_gee_phyc_treewas_intersection_sig_hits_intersection.tsv", sep = ""))
    save_results(gee_phyc_intersection,         paste(subdir, phenotypes[p], genotypes[g],"_gee_phyc_intersection_sig_hits_intersection.tsv", sep = ""))
    save_results(gee_treewas_intersection,      paste(subdir, phenotypes[p], genotypes[g],"_gee_treewas_intersection_sig_hits_intersection.tsv", sep = ""))
    save_results(phyc_treewas_intersection,     paste(subdir, phenotypes[p], genotypes[g],"_phyc_treewas_intersection_sig_hits_intersection.tsv", sep = ""))
    
    # Update summary information
    name <- paste(phenotypes[p], genotypes[g], sep = "")
    intersection_all_methods <- append_findings(intersection_all_methods, name, gee_phyc_treewas_intersection)
    all_two_plus <- unique(c(gee_phyc_intersection, gee_treewas_intersection, phyc_treewas_intersection))
    intersection_two_plus_methods <- append_findings(intersection_two_plus_methods, name, all_two_plus)
    all_hits <- unique(c(gee, phyc, treewas))
    any_hits <- append_findings(any_hits, name, all_hits)
  }
}
# Save the big summary files
save_list(intersection_all_methods,      paste(main_dir, "all_methods_intersection_per_test_sig_hits.txt",      sep = ""))
save_list(intersection_two_plus_methods, paste(main_dir, "two_plus_methods_intersection_per_test_sig_hits.txt", sep = ""))
save_list(any_hits,                      paste(main_dir, "any_per_test_sig_hits.txt",                           sep = ""))

save(intersection_all_methods,       file = paste(main_dir, "all_methods_intersection_per_test_sig_hits.rda",      sep = ""))
save(intersection_two_plus_methods,  file = paste(main_dir, "two_plus_methods_intersection_per_test_sig_hits.rda", sep = ""))
save(any_hits,                       file = paste(main_dir, "any_per_test_sig_hits.rda",                           sep = ""))

# END SCRIPT -------------------------------------------------------------------