library(ape)
library(phytools)
set.seed(1)
tree <- rtree(7)
tree$edge.length <- rep(sum(tree$edge.length)/Nedge(tree), Nedge(tree))
tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

discrete_phenotype <- matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = Ntip(tree), ncol = 1)
row.names(discrete_phenotype) <- tree$tip.label
colnames(discrete_phenotype) <- "resistance"

set.seed(1)
continuous_phenotype <- as.matrix(fastBM(tree))
row.names(continuous_phenotype) <- tree$tip.label
colnames(continuous_phenotype) <- "growth"

genotype1 <- matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = Ntip(tree), ncol = 1)
genotype2 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = Ntip(tree), ncol = 1)
genotype3 <- matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = Ntip(tree), ncol = 1)
genotype4 <- matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = Ntip(tree), ncol = 1)
genotype5 <- matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = Ntip(tree), ncol = 1)
genotype6 <- matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = Ntip(tree), ncol = 1)
genotype7 <- matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = Ntip(tree), ncol = 1)
genotype8 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = Ntip(tree), ncol = 1)
genotype <- cbind(genotype1, genotype2, genotype3, genotype4, genotype5, genotype6, genotype7, genotype8)
row.names(genotype) <- tree$tip.label
colnames(genotype) <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

snp_gene_key <- matrix(NA, nrow = ncol(genotype), ncol = 2)
colnames(snp_gene_key) <- c("SNP", "GENE")
snp_gene_key[ , 1] <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
snp_gene_key[ , 2] <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")

write.table(genotype, file =  "/nfs/esnitkin/bin_group/pipeline/Github/gwas/hogwash/data/grouped_genotype.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(discrete_phenotype, file = "/nfs/esnitkin/bin_group/pipeline/Github/gwas/hogwash/data/discrete_phenotype.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(continuous_phenotype, file = "/nfs/esnitkin/bin_group/pipeline/Github/gwas/hogwash/data/continuous_phenotype.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(snp_gene_key, file = "/nfs/esnitkin/bin_group/pipeline/Github/gwas/hogwash/data/snp_gene_key.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.tree(tree, file = "/nfs/esnitkin/bin_group/pipeline/Github/gwas/hogwash/data/tree.tree")
