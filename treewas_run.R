# Katie Saund

# LOAD LIBRARIES --------------------------------------------------------------#
source("/nfs/esnitkin/bin_group/pipeline/Github/gwas/treewas_lib.R")

# READ IN ARGUMENTS -----------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE) #arguments from the PBS script
pheno     <- read.csv(args[1], row.names = 1, sep = "\t")
tree      <- read.tree(args[2])
geno      <- read.csv(args[3], header = TRUE, row.names = 1, sep = "\t")
recon     <- args[4]
test_corr <- args[5]
fname     <- args[6]
alpha     <- args[7]

# RUN TREEWAS -----------------------------------------------------------------#
temp <- unlist(pheno)
names(temp) <- row.names(pheno)
pheno <- temp

if (!is.rooted(tree)){
  tree <- midpoint(tree)
}

pdf_fname <- paste(fname, ".pdf", sep = "")
geno <- remove_nonexistant_genotypes(geno)
results <- run_treeWAS(tree, geno, pheno, recon, test_corr, pdf_fname, alpha)
save(results, file = paste(fname, "_results.rda", sep = ""))
save_treewas_hits(results$terminal$sig.snps, "terminal", fname)
save_treewas_hits(results$simultaneous$sig.snps, "simultaneous", fname)
save_treewas_hits(results$subsequent$sig.snps, "subsequent", fname)
save_treewas_hits(results$treeWAS.combined$treeWAS.combined, "combined", fname)

# SAVE SESSION INFO -----------------------------------------------------------#
writeLines(capture.output(sessionInfo()), paste(fname, "_sessionInfo_", format(Sys.Date(), "%Y-%m-%d"), ".txt", sep = ""))

# END OF SCRIPT ---------------------------------------------------------------#
