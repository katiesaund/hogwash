# Katie Saund
# 2018-09-05
# as of 2018-09-05 left out pilon and snp_del results b/c the genotypes take too long to read into the file?
# TO DO: I added 3 repetitions of the if() statement that make me open up the genotype 3x instead of 1 time and adds redudnancy, go back to fix this. 

# Libraries --------------------------------------------------------------------
source("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/lib/2018-09-06_plot_gwas_hits_lib.R")
# ------------------------------------------------------------------------------

phenotypes <- c("log_cfe", "log_growth", "log_toxin", "log_germ_tc_and_gly", "log_germ_tc", "log_sporulation", "fqR", "severity") 
genotypes  <- c("_gene_ns", "_gene_stop", "_gene_high", "_gene_del", "_snp_stop","_snp_high", "_roary_pan_genome", "_snp_del" , "_pilon_sv") #left out _snp_ns because as of 2018-09-06 9pm snp_ns is still not done running for wany GWAS method

# change path to getwd()
results_path <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/compare_methods/"
gwas_format_path <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/"
gwas_intersection_results <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/compare_methods/all_methods_intersection_per_test_sig_hits.rda"
gwas_two_method_results <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/compare_methods/two_plus_methods_intersection_per_test_sig_hits.rda"
gwas_any_method_results <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/compare_methods/any_per_test_sig_hits.rda"

# args <- commandArgs(trailingOnly = TRUE) # Grab arguments from the PBS script
# results_path <- args[1]
# gwas_format_path <- args[2]
# gwas_intersection_results <- args[3]

# CHECK INPUTS -----------------------------------------------------------------
check_if_dir_exists(results_path)
check_if_dir_exists(gwas_format_path)
check_file_exists(gwas_intersection_results)
load(gwas_intersection_results)
load(gwas_two_method_results)
load(gwas_any_method_results)

# FUNCTION ---------------------------------------------------------------------
for (p in 1:length(phenotypes)){
  phenotype <- read_in_tsv_matrix(paste(gwas_format_path, phenotypes[p], "_pheno.tsv", sep = ""))
  colnames(phenotype) <- phenotypes[p]
  tree    <- read.tree(paste(gwas_format_path, phenotypes[p], ".tree", sep = ""))
  htmp_tr <- create_heatmap_compatible_tree(tree)
  annotation <- read_in_tsv_matrix(paste(gwas_format_path, phenotypes[p], "_annotation.tsv", sep = ""))
  simple_annotation <- annotation[ , 1, drop = FALSE]
  colnames(simple_annotation) <- "Ribotype"
  discrete_or_continuous <- assign_pheno_type(phenotype)
  
  for (g in 1:length(genotypes)){
    print(paste(phenotypes[p], genotypes[g], sep = ""))
    name <- paste(phenotypes[p], genotypes[g], sep = "")
    for (i in 1:length(names(intersection_all_methods))){
      if(names(intersection_all_methods)[i] == name){
        print("start all methods")
        genotype <- load_rda(paste(gwas_format_path, phenotypes[p], genotypes[g], ".rda", sep = ""))
        genotype <- format_genotype_name(genotype)
        #genotype <- read_in_tsv_matrix(paste(gwas_format_path, phenotypes[p], genotypes[g], ".tsv", sep = ""))
        method_names <- gsub("^X", "", unlist(intersection_all_methods[i]))
        hits <- genotype[ , colnames(genotype) %in%  method_names, drop = FALSE]
        white_to_black = colorRampPalette(c("black", "white")) # Set up presence/absence colors
        num_color = length(min(hits):max(hits)) # We want a unique shade of grayscale for each level (should be 2 for presence/absence)
        htmp_col = white_to_black(num_color) # Generate a character vector of colors with each shade of grayscale
        
        if (discrete_or_continuous == "continuous"){
          ha_combined = rowAnnotation(df = simple_annotation, 
                                      boxplot = row_anno_barplot(phenotype, axis = TRUE, gp = gpar(fill = "black")), 
                                      col = list(Ribotype = c("027" = "blue", "014-020" = "purple", "053-163" = "orange", "078-126" = "green", "other" = "grey")), 
                                      width = unit(40, "mm"),
                                      annotation_width = c(.5,10))
          
        } else if (discrete_or_continuous == "discrete") {
          discrete_phenotype <- phenotype
          colnames(discrete_phenotype) <- "Presence"

          combined_annotation <- data.frame(Ribotype = simple_annotation, presence = discrete_phenotype)
          ha_combined <- rowAnnotation(df = combined_annotation, 
                                      col = list(Ribotype= c("027" = "blue", "014-020" = "purple", "053-163" = "orange", "078-126" = "green", "other" = "grey"), 
                                                 Presence= c("0" = "white", "1" = "black")), 
                                      width = unit(30, "mm"),
                                      annotation_width = c(1, 1))
        } else {
          stop("Phenotype not defined as either discrete nor continuous.")
        }
        htmp <- ComplexHeatmap::Heatmap(
          matrix = hits, 
          cluster_columns = FALSE, 
          row_names_side = "left", 
          col = htmp_col, 
          show_column_names = TRUE, 
          row_names_max_width = unit(20, "mm"),
          show_heatmap_legend = FALSE, 
          column_names_max_height = unit(50, "mm"),
          row_title = "Isolates",
          column_title = paste("Hits from all 3 GWAS methods: ", phenotypes[p], genotypes[g], sep = ""),  
          cluster_rows = htmp_tr, 
          row_dend_side = "left", 
          row_dend_width = unit(40, "mm"), 
          row_names_gp = gpar(cex = 0.4), 
          column_names_gp = gpar(cex = 0.75), width = unit(2 *ncol(hits), "mm")
        ) + ha_combined
        
        pdf(paste(results_path, "all_methods_intersection_heatmap_", phenotypes[p], genotypes[g], ".pdf", sep = ""))
        draw(htmp)
        dev.off()
        print("done all methods")
      }
      #
      #
      #
    }
    for (i in 1:length(names(intersection_two_plus_methods))){
      if(names(intersection_two_plus_methods)[i] == name){
        print("start two methods")
        genotype <- load_rda(paste(gwas_format_path, phenotypes[p], genotypes[g], ".rda", sep = ""))
        genotype <- format_genotype_name(genotype)
        #genotype <- read_in_tsv_matrix(paste(gwas_format_path, phenotypes[p], genotypes[g], ".tsv", sep = ""))
        method_names <- gsub("^X", "", unlist(intersection_two_plus_methods[i]))
        #hits <- genotype[ , colnames(genotype) %in%  unlist(intersection_two_plus_methods[i]), drop = FALSE]
        hits <- genotype[ , colnames(genotype) %in%  method_names, drop = FALSE]
        
        white_to_black = colorRampPalette(c("black", "white")) # Set up presence/absence colors
        num_color = length(min(hits):max(hits)) # We want a unique shade of grayscale for each level (should be 2 for presence/absence)
        htmp_col = white_to_black(num_color) # Generate a character vector of colors with each shade of grayscale
        
        if (discrete_or_continuous == "continuous"){
          ha_combined = rowAnnotation(df = simple_annotation, 
                                      boxplot = row_anno_barplot(phenotype, axis = TRUE, gp = gpar(fill = "black")), 
                                      col = list(Ribotype = c("027" = "blue", "014-020" = "purple", "053-163" = "orange", "078-126" = "green", "other" = "grey")), 
                                      width = unit(40, "mm"),
                                      annotation_width = c(.5,10))
          
        } else if (discrete_or_continuous == "discrete") {
          discrete_phenotype <- phenotype
          colnames(discrete_phenotype) <- "Presence"
          
          combined_annotation <- data.frame(Ribotype = simple_annotation, presence = discrete_phenotype)
          ha_combined <- rowAnnotation(df = combined_annotation, 
                                       col = list(Ribotype= c("027" = "blue", "014-020" = "purple", "053-163" = "orange", "078-126" = "green", "other" = "grey"), 
                                                  Presence= c("0" = "white", "1" = "black")), 
                                       width = unit(30, "mm"),
                                       annotation_width = c(1, 1))
        } else {
          stop("Phenotype not defined as either discrete nor continuous.")
        }
        htmp <- ComplexHeatmap::Heatmap(
          matrix = hits, 
          cluster_columns = FALSE, 
          row_names_side = "left", 
          col = htmp_col, 
          show_column_names = TRUE, 
          row_names_max_width = unit(20, "mm"),
          show_heatmap_legend = FALSE, 
          column_names_max_height = unit(50, "mm"),
          row_title = "Isolates",
          column_title = paste("Hits from at least 2 GWAS methods: ", phenotypes[p], genotypes[g], sep = ""),  
          cluster_rows = htmp_tr, 
          row_dend_side = "left", 
          row_dend_width = unit(40, "mm"), 
          row_names_gp = gpar(cex = 0.4), 
          column_names_gp = gpar(cex = 0.75), width = unit(2 *ncol(hits), "mm")
        ) + ha_combined
        
        pdf(paste(results_path, "two_methods_intersection_heatmap_", phenotypes[p], genotypes[g], ".pdf", sep = ""))
        draw(htmp)
        dev.off()
        print("done two methods")

      }
    }
    for (i in 1:length(names(any_hits))){
      # repeat for any method
      if(names(any_hits)[i] == name){
        print("any hits")
        genotype <- load_rda(paste(gwas_format_path, phenotypes[p], genotypes[g], ".rda", sep = ""))
        genotype <- format_genotype_name(genotype)
        #genotype <- read_in_tsv_matrix(paste(gwas_format_path, phenotypes[p], genotypes[g], ".tsv", sep = ""))
        
        method_names <- gsub("^X", "", unlist(any_hits[i]))
        hits <- genotype[ , colnames(genotype) %in%  method_names, drop = FALSE]
        white_to_black = colorRampPalette(c("black", "white")) # Set up presence/absence colors
        num_color = length(min(hits):max(hits)) # We want a unique shade of grayscale for each level (should be 2 for presence/absence)
        htmp_col = white_to_black(num_color) # Generate a character vector of colors with each shade of grayscale
        
        if (discrete_or_continuous == "continuous"){
          ha_combined = rowAnnotation(df = simple_annotation, 
                                      boxplot = row_anno_barplot(phenotype, axis = TRUE, gp = gpar(fill = "black")), 
                                      col = list(Ribotype = c("027" = "blue", "014-020" = "purple", "053-163" = "orange", "078-126" = "green", "other" = "grey")), 
                                      width = unit(40, "mm"),
                                      annotation_width = c(.5,10))
          
        } else if (discrete_or_continuous == "discrete") {
          discrete_phenotype <- phenotype
          colnames(discrete_phenotype) <- "Presence"
          
          combined_annotation <- data.frame(Ribotype = simple_annotation, presence = discrete_phenotype)
          ha_combined <- rowAnnotation(df = combined_annotation, 
                                       col = list(Ribotype= c("027" = "blue", "014-020" = "purple", "053-163" = "orange", "078-126" = "green", "other" = "grey"), 
                                                  Presence= c("0" = "white", "1" = "black")), 
                                       width = unit(30, "mm"),
                                       annotation_width = c(1, 1))
        } else {
          stop("Phenotype not defined as either discrete nor continuous.")
        }
        htmp <- ComplexHeatmap::Heatmap(
          matrix = hits, 
          cluster_columns = FALSE, 
          row_names_side = "left", 
          col = htmp_col, 
          show_column_names = TRUE, 
          row_names_max_width = unit(20, "mm"),
          show_heatmap_legend = FALSE, 
          column_names_max_height = unit(50, "mm"),
          row_title = "Isolates",
          column_title = paste("Hits from any GWAS method: ", phenotypes[p], genotypes[g], sep = ""),  
          cluster_rows = htmp_tr, 
          row_dend_side = "left", 
          row_dend_width = unit(40, "mm"), 
          row_names_gp = gpar(cex = 0.4), 
          column_names_gp = gpar(cex = 0.75), width = unit(2 *ncol(hits), "mm")
        ) + ha_combined
        
        pdf(paste(results_path, "any_method_intersection_heatmap_", phenotypes[p], genotypes[g], ".pdf", sep = ""))
        draw(htmp)
        dev.off()
        print("done with any")
      }
    }
  }
}

# END OF SCRIPT ----------------------------------------------------------------
