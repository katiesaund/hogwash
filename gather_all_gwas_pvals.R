# Katie Saund
# 2018-09-07

# Script will gather the pvalue for each genotype-phenotype results after 
# running GEE, treeWAS, and phyC. 

library(ComplexHeatmap)

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

load_rda <- function(file_name){
  #loads an RData file, and returns it
  load(file_name)
  get(ls()[ls() != "file_name"])
} # end load_rda


get_gee_pvals <- function(path, pheno, geno){
  file_name <- paste(path, "gee_", pheno, geno, "_all_pval.tsv", sep = "")
  pvals <- read_in_tsv_matrix(file_name)
  return(pvals)
} # end get_gee_pvals()

get_phyc_trans_pvals <- function(path, pheno, geno){
  if (pheno %in% c("fqR", "severity")){ # binary phenotypes
    file_name_trans <- paste(path, "phyc_", pheno, geno, "_transition_pvals_for_all_hits.tsv", sep = "")
    if (file.exists(file_name_trans)){
      pvals <- read_in_tsv_matrix(file_name_trans)
    } 
  } else { #log_toxin, log_growth, etc... Continuous phenotypes
    file_name_trans <- paste(path, "phyc_", pheno, geno, "_all_transitions_pvals_for_all_hits.tsv", sep = "")
    if (file.exists(file_name_trans)){
      pvals <- read_in_tsv_matrix(file_name_trans)
    }
  }
  return(pvals)
} # end get_phyc_trans_pvals()

get_treewas_pvals <- function(path, pheno, geno){
  treewas_data <- load_rda(paste(path, "treewas_", pheno, geno, "_results.rda", sep = ""))
  ter_pval <- treewas_data$terminal$p.vals
  sub_pval <- treewas_data$subsequent$p.vals
  sim_pval <- treewas_data$simultaneous$p.vals
  names(ter_pval) <- names(sub_pval) <- names(sim_pval) <- colnames(treewas_data$dat$snps)
  ter_pval <- as.matrix(ter_pval)
  sub_pval <- as.matrix(sub_pval)
  sim_pval <- as.matrix(sim_pval)
  return(list("ter_pval" = ter_pval, "sub_pval" = sub_pval, "sim_pval" = sim_pval))
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

get_expected_gene_subset <- function(pval_mat, gene_list){
  if (class(pval_mat) == "matrix"){
    subset <- rep(FALSE, nrow(pval_mat))
    for (i in 1:length(gene_list)){
      for (j in 1:nrow(pval_mat)){
        if (length(grep(gene_list[i], row.names(pval_mat)[j], ignore.case = TRUE)) != 0){
          subset[j] <- TRUE
        }
      }
    }
  } else {
    stop("input bad")
  }
  controls <- pval_mat[ subset, , drop = FALSE]
  return(controls)
}    




subset_all_gwas_types <- function(gee_results, phyc_results, ter_results, sim_results, sub_results, gene_list){
  gee_pos_ctrl  <- get_expected_gene_subset(gee_results, gene_list)
  phyc_pos_ctrl <- get_expected_gene_subset(phyc_results, gene_list)
  ter_pos_ctrl  <- get_expected_gene_subset(ter_results, gene_list)
  sim_pos_ctrl  <- get_expected_gene_subset(sim_results, gene_list)
  sub_pos_ctrl  <- get_expected_gene_subset(sub_results, gene_list)
  return(list("gee_pos_ctrl" = gee_pos_ctrl, "phyc_pos_ctrl" = phyc_pos_ctrl, 
              "ter_pos_ctrl" = ter_pos_ctrl, "sim_pos_ctrl" = sim_pos_ctrl, 
              "sub_pos_ctrl" = sub_pos_ctrl))
}


collect_pvals_for_found_genes <- function(gee_pc, phyc_pc, ter_pc, sim_pc, sub_pc){
  found_genes <- unique(c(row.names(gee_pc), 
                          row.names(phyc_pc), 
                          row.names(ter_pc), 
                          row.names(sim_pc), 
                          row.names(sub_pc)))
  
  found_genes_with_pvals <- matrix(NA, nrow = length(found_genes), ncol = 5)
  print(found_genes)
  if (length(found_genes) == 0){
    stop("no found genes")
  }
  colnames(found_genes_with_pvals) <- c("gee", "phyc", "ter", "sim", "sub")
  row.names(found_genes_with_pvals) <- found_genes
  
  for (i in 1:nrow(found_genes_with_pvals)){
    for (j in 1:nrow(gee_pc)){
      if (row.names(found_genes_with_pvals)[i] %in% row.names(gee_pc)[j]){
        found_genes_with_pvals[i, 1] <- gee_pc[j, 29]
      }
    }
    for (j in 1:nrow(phyc_pc)){
      if (row.names(found_genes_with_pvals)[i] %in% row.names(phyc_pc)[j]){
        found_genes_with_pvals[i, 2] <- phyc_pc[j, 1]
      }
    }
    for (j in 1:nrow(ter_pc)){
      if (row.names(found_genes_with_pvals)[i] %in% row.names(ter_pc)[j]){
        found_genes_with_pvals[i, 3] <- ter_pc[j, 1]
      }
    }
    for (j in 1:nrow(sim_pc)){
      if (row.names(found_genes_with_pvals)[i] %in% row.names(sim_pc)[j]){
        found_genes_with_pvals[i, 4] <- sim_pc[j, 1]
      }
    }
    for (j in 1:nrow(sub_pc)){
      if (row.names(found_genes_with_pvals)[i] %in% row.names(sub_pc)[j]){
        found_genes_with_pvals[i, 5] <- sub_pc[j, 1]
      }
    }
  }
  return(found_genes_with_pvals)
}

plot_found_genes <- function(pval_mat, name, path_name){
  temp <- pval_mat
  pval_mat[pval_mat >= 0.01] <- 1
  pval_mat[pval_mat < 0.01]  <- 0
  print(table(pval_mat))
  print(dim(table(pval_mat)))
  temp_name <- gsub("\n", "_", name)
  temp_name <- gsub(" ", "_", temp_name)
  if (dim(table(pval_mat)) > 1){
    white_to_black = colorRampPalette(c("black", "white")) # Set up presence/absence colors
    #num_color = length(min(pval_mat):max(pval_mat)) # We want a unique shade of grayscale for each level (should be 2 for presence/absence)
    htmp_col = white_to_black(2) # Generate a character vector of colors with each shade of grayscale
    htmp <- ComplexHeatmap::Heatmap(
      matrix = pval_mat,
      na_col = "grey",
      cluster_columns = FALSE, 
      row_names_side = "left", 
      col = htmp_col, 
      show_column_names = TRUE, 
      row_names_max_width = unit(50, "mm"),
      show_heatmap_legend = FALSE, 
      column_names_max_height = unit(50, "mm"),
      row_title = "Positive Control Genes",
      column_title = paste(name, "\nBlack=Significant White=Not Sign. Gray=NA", sep = " "),  
      cluster_rows = TRUE, 
      show_row_dend = TRUE,
      row_names_gp = gpar(cex = 1.0), 
      column_names_gp = gpar(cex = 1.0), width = unit(5 *ncol(pval_mat), "mm")
    )

    pdf(paste(path_name, "/", temp_name, ".pdf", sep = ""))
    draw(htmp)
    dev.off()
  }
  write.table(pval_mat, file = paste(path_name, "/", temp_name, ".tsv", sep = ""), sep = "\t")
  
}


all_subset_steps <- function(gee_pvals, phyc_pvals, ter_pvals, sim_pvals, sub_pvals, gene_list, plot_name, dir_name){
  sub_results <- subset_all_gwas_types(gee_pvals, phyc_pvals, ter_pvals, sim_pvals, sub_pvals, gene_list)
  pvals <- collect_pvals_for_found_genes(sub_results$gee_pos_ctrl, sub_results$phyc_pos_ctrl, sub_results$ter_pos_ctrl, sub_results$sim_pos_ctrl, sub_results$sub_pos_ctrl)
  plot_found_genes(pvals, plot_name, dir_name)
}

# FUNCTION ---------------------------------------------------------------------

expected_sporulation_genes <- c("CD1352", "CD1492", "CD1579", "CD1949", 
                                "CD2492", "Spo0A", "SpoIIAA", "SpoIIAB", 
                                "SpoIIE", "Sigma factor F", "Sigma factor G", 
                                "Sigma factor E", "Sigma factor K", "SpoIVA", 
                                "SpoVM", "SipL", "CdeC", "BclA1", "bclA2", 
                                "bclA3", "CD3563", "CspBA", "SleC", "CspC", 
                                "CotA", "CotB", "CotCB", "CotD", "CotE", "CotF", 
                                "CotJB2", "CotG", "SodA", "cheB", "sigE", 
                                "sigF", "sigG", "spoIIAA", "spoIIAB", "spoIVB", 
                                "gpr", "spoIIE", "spoIIR", "spoIVA", "spoIVB2")

expected_germination_genes <- c("gerP", "cwlJ", "gerKA", "gerKB", "gerKC", 
                                "gerLA", "gerLB", "gerLC", "spoVA", "gerS", 
                                "cspA", "cspC", "cspB", "cspBA", "sipL", "ssb", 
                                "CD3571", "asd", "CD2894", "spoIIIAA", "CD0972")

expected_toxin_genes       <- c("rho", "tpi", "CD2808", "tcdA", "tcdC", "codY", 
                                "ldh", "polA", "CD1263", "tcdE", "tcdB", "cdtR", 
                                "ermB", "gluD", "gyrA", "CD2894A")
expected_fqR_genes         <- c("gyrA", "gyrB")



path <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/"

# for running snp_high: Make sure to comment out germ section
phenotypes <- c("log_cfe", "log_toxin", "log_sporulation") # "log_growth" "severity" "fqR" "log_germ_tc_and_gly" "log_germ_tc"
genotypes  <- c( "_snp_high") #, "_gene_ns")     # "_gene_stop", ,"_snp_ns", "_gene_high",  "_gene_del", "_snp_stop",", "_snp_del", "_pilon_sv", "_roary_pan_genome"

# for running gene_ns - make sure to uncomment germ section
# phenotypes <- c("log_cfe", "log_germ_tc_and_gly", "log_toxin", "log_sporulation", "log_germ_tc", "fqR") # "log_growth" "severity"
# genotypes  <- c( "_gene_ns") #, "_snp_high")     # "_gene_stop", ,"_snp_ns", "_gene_high",  "_gene_del", "_snp_stop",", "_snp_del", "_pilon_sv", "_roary_pan_genome"


pos_ctrl_dir <- paste(path, "/positive_controls/", sep = "")
system(paste("mkdir ", pos_ctrl_dir, sep = ""))


for (p in 1:length(phenotypes)){
  for (g in 1:length(genotypes)){
    print("loop start")
    print(paste(phenotypes[p], genotypes[g], sep = ""))
    
    gee     <- t(get_gee_pvals(path, phenotypes[p], genotypes[g])) # ignores if QIC is low  # gee isa matrix. # colnames are genotype, transpose to orient the same as phyC
    phyc    <- get_phyc_trans_pvals(path, phenotypes[p], genotypes[g]) #phyc is a matrix # need to transpose to make colnames genotype
    treewas <- get_treewas_pvals(path, phenotypes[p], genotypes[g])#treewas$ter_pval # vector #treewas$sub_pval #treewas$sim_pval

    if (phenotypes[p] == "log_toxin"){
      print("toxin")
      name <- paste(phenotypes[p], genotypes[g], "\nexpected_toxin_genes", sep = "")
      all_subset_steps(gee, phyc, treewas$ter_pval, treewas$sim_pval, treewas$sub_pval, expected_toxin_genes, name, pos_ctrl_dir)
    }
    
    if (phenotypes[p] %in% c("log_sporulation", "log_cfe")){
      print("spore")
      name <- paste(phenotypes[p], genotypes[g], "\nexpected_sporulation_genes", sep = "")
      all_subset_steps(gee, phyc, treewas$ter_pval, treewas$sim_pval, treewas$sub_pval, expected_sporulation_genes, name, pos_ctrl_dir)
    }

    # if (phenotypes[p] %in% c("log_germ_tc_and_gly", "log_germ_tc")){
    #   print("Germ")
    #   name <- paste(phenotypes[p], genotypes[g], "\nexpected_germination_genes", sep = "")
    #   all_subset_steps(gee, phyc, treewas$ter_pval, treewas$sim_pval, treewas$sub_pval, expected_germination_genes, name, pos_ctrl_dir)
    # }

    if (phenotypes[p] == "fqR"){
      print("fq")
      name <- paste(phenotypes[p], genotypes[g], "\nexpected_fqR_genes", sep = "")
      all_subset_steps(gee, phyc, treewas$ter_pval, treewas$sim_pval, treewas$sub_pval, expected_fqR_genes, name, pos_ctrl_dir)
    }
  }
}

# END SCRIPT -------------------------------------------------------------------



