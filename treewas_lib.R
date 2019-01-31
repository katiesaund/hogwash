# Katie Saund

library(treeWAS)
library(ape)
library(phangorn)

run_treeWAS <- function(tree, variant_matrix, phenotype, reconstruction_type, multiple_test_correction, filename, alpha_value){
  # tree is a rooted phylogenetic object
  # variant_matrix is a matrix of variants (what is proper orientation?)
  # phenotype is a vector??
  # reconstruction_type is a string ("parsimony")
  # multiple_test_correction is a string ("bonf" or "fdr")
  if (class(reconstruction_type) != typeof(reconstruction_type) | class(reconstruction_type) != "character"){
    stop("Reconstruction type should be a string (eg \"parsimony\")")
  }
  if (class(multiple_test_correction) != typeof(multiple_test_correction) | class(multiple_test_correction) != "character"){
    stop("multiple_test_correction should be a string (eg \"bonf\" or \"fdr\")")
  }
  if (class(tree) != "phylo" | !(is.rooted(tree))){
    stop("tree should be a rooted phylo object")
  }
  
  
  results <- treeWAS(snps = variant_matrix, 
                     phen = phenotype, 
                     tree = tree, 
                     n.subs = NULL, 
                     n.snps.sim = ncol(variant_matrix)*10, 
                     chunk.size = ncol(variant_matrix), 
                     test = c("terminal", "simultaneous", "subsequent"), 
                     snps.reconstruction = reconstruction_type, 
                     snps.sim.reconstruction = reconstruction_type, 
                     phen.reconstruction = reconstruction_type, 
                     phen.type = NULL,
                     na.rm = TRUE, 
                     p.value = alpha_value, 
                     p.value.correct = multiple_test_correction, 
                     p.value.by = "count", 
                     dist.dna.model = "JC69", 
                     plot.tree = TRUE, 
                     plot.manhattan = TRUE, 
                     plot.null.dist = TRUE, 
                     plot.dist = TRUE, 
                     snps.assoc = NULL, 
                     filename.plot = filename, 
                     seed = 1)
  return(results)
}

save_treewas_hits <- function(result_list, type, fname){
  write.table(result_list, 
              file = paste(fname, "_", type, "_sig_hits.tsv", sep = ""),
              sep = "\t",
              col.names = TRUE, 
              row.names = TRUE,
              quote = FALSE)
}

remove_nonexistant_genotypes <- function(geno_mat, prefix){
  geno_to_drop <- rep(FALSE, ncol(geno_mat))
  geno_to_drop[colSums(geno_mat) < 1] <- TRUE
  dropped_genotype_names <- colnames(geno_mat)[geno_to_drop]
  filename <- paste(prefix, "_genotypes_dropped_because_vars_not_in_samples.txt", sep = "")
  writeLines(dropped_genotype_names, con = filename, sep = "\n")
  geno_mat <- geno_mat[ , !geno_to_drop, drop = FALSE]
  return(geno_mat)
} # end remove_nonexistant_genotypes