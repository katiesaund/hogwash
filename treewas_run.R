# 2018-08-24
# Katie Saund

# LOAD LIBRARIES -----------------------------------------------------------------------------------------
library(ape)
library(phangorn)
source("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-24_treewas/lib/2018-08-24_treewas_lib.R")

# READ IN ARGUMENTS --------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE) #arguments from the PBS script
pheno <- read.csv(args[1], row.names = 1, sep = "\t")
tree <- read.tree(args[2])
geno <- read.csv(args[3], header = TRUE, row.names = 1, sep = "\t")
recon <- args[4]
test_corr <- args[5]
fname <- args[6]

# RUN TREEWAS --------------------------------------------------------------------------------------------
temp <- unlist(pheno)
names(temp) <- row.names(pheno)
pheno <- temp

# added is.rooted if statement on 2018-09-25 to deal with odd midpoint rooting issue for PSM dataset. 
if (!is.rooted(tree)){
  tree <- midpoint(tree)
}
# End of section added on 2018-09-25

pdf_fname <- paste(fname, ".pdf", sep = "")

results <- run_treeWAS(tree, geno, pheno, recon, test_corr, pdf_fname)
save(results, file = paste(fname, "_results.rda", sep = ""))
save_treewas_hits(results$terminal$sig.snps, "terminal", fname)
save_treewas_hits(results$simultaneous$sig.snps, "simultaneous", fname)
save_treewas_hits(results$subsequent$sig.snps, "subsequent", fname)
save_treewas_hits(results$treeWAS.combined$treeWAS.combined, "combined", fname)

# END OF SCRIPT ------------------------------------------------------------------------------------------
