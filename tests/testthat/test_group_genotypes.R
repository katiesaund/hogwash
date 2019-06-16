context("Group genotypes") #---------------------------------------------------#

# test build_gene_trans_from_snp_trans -----------------------------------------
test_that("build_gene_trans_from_snp_trans does X given Y", {
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <- rep(sum(temp_tree$edge.length)/ape::Nedge(temp_tree), ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <- matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <- matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <- matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <- matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype <- cbind(genotype1, genotype2, genotype5, genotype6, genotype7, genotype8)
  row.names(genotype) <- temp_tree$tip.label
  colnames(genotype) <- c("SNP1", "SNP2", "SNP5", "SNP6", "SNP7", "SNP8")

  set.seed(1)
  continuous_phenotype <- as.matrix(phytools::fastBM(temp_tree))
  row.names(continuous_phenotype) <- temp_tree$tip.label
  colnames(continuous_phenotype) <- "growth"


  snp_gene_key <- matrix(NA, nrow = ncol(genotype), ncol = 2)
  colnames(snp_gene_key) <- c("SNP", "GENE")
  snp_gene_key[ , 1] <- c("SNP1", "SNP2", "SNP5", "SNP6", "SNP7", "SNP8")
  snp_gene_key[ , 2] <- c("GENE1", "GENE2", "GENE5", "GENE6", "GENE7", "GENE7")

  AR <- prepare_ancestral_reconstructions(temp_tree, continuous_phenotype, genotype, "continuous")

  temp_results <- build_gene_trans_from_snp_trans(temp_tree, genotype, AR$geno_trans, snp_gene_key)

  expect_true(length(temp_results) == length(unique(snp_gene_key[ , 2])))
  expect_equal(temp_results[[1]]$transition, AR$geno_trans[[1]]$transition)
  expect_equal(temp_results[[2]]$transition, AR$geno_trans[[2]]$transition)
  expect_equal(temp_results[[3]]$transition, AR$geno_trans[[3]]$transition)
  expect_equal(temp_results[[4]]$transition, AR$geno_trans[[4]]$transition)
  expect_equal(temp_results[[5]]$transition, AR$geno_trans[[5]]$transition + AR$geno_trans[[6]]$transition)

  expect_equal(length(temp_results[[1]]$transition), ape::Nedge(temp_tree))
  expect_equal(length(temp_results[[1]]$trans_dir), ape::Nedge(temp_tree))

})

# test prepare_genotype --------------------------------------------------------
test_that("prepare_genotype gives expected genotype for a grouped input", {
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <- rep(sum(temp_tree$edge.length)/ape::Nedge(temp_tree), ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <- matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <- matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <- matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <- matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <- matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <- matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1, genotype2, genotype3, genotype4, genotype5, genotype6, genotype7, genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[ , 1] <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[ , 2] <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")

  temp_logical <- TRUE

  temp_result <- prepare_genotype(temp_logical, temp_geno, temp_tree, temp_lookup)
  expect_equal(temp_result$genotype, temp_geno[ , c(1, 2, 5, 6, 7, 8)])
  expect_equal(temp_result$gene_snp_lookup, temp_lookup[c(1, 2, 5, 6, 7, 8), ])
  expect_equal(as.numeric(unname(temp_result$snps_per_gene)), c(1, 1, 1, 1, 2))
  expect_equal(temp_result$unique_genes, c("GENE1", "GENE2", "GENE5", "GENE6", "GENE7"))

})
test_that("prepare_genotype gives expected genotype for an not grouped input", {
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <- rep(sum(temp_tree$edge.length)/ape::Nedge(temp_tree), ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <- matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <- matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <- matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <- matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <- matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <- matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1, genotype2, genotype3, genotype4, genotype5, genotype6, genotype7, genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[ , 1] <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[ , 2] <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")
  temp_logical <- FALSE
  temp_result <- prepare_genotype(temp_logical, temp_geno, temp_tree, temp_lookup)

  expect_null(temp_result$snps_per_gene)
  expect_equal(temp_result$genotype, temp_geno[ , c(1, 5, 6)])
  expect_equal(temp_result$convergence_not_possible_genotypes, c("SNP2", "SNP3", "SNP4", "SNP7", "SNP8"))

})

# test prepare_ungrouped_genotype ----------------------------------------------
test_that("prepare_ungrouped_genotype", {
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <- rep(sum(temp_tree$edge.length)/ape::Nedge(temp_tree), ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <- matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <- matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <- matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <- matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <- matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <- matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1, genotype2, genotype3, genotype4, genotype5, genotype6, genotype7, genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_result <- prepare_ungrouped_genotype(temp_geno, temp_tree)
  expect_null(temp_result$snps_per_gene)
  expect_equal(temp_result$genotype, temp_geno[ , c(1, 5, 6)])
  expect_equal(temp_result$convergence_not_possible_genotypes, c("SNP2", "SNP3", "SNP4", "SNP7", "SNP8"))
})

# test prepare_grouped_genotype ------------------------------------------------
test_that("prepare_grouped_genotype", {
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <- rep(sum(temp_tree$edge.length)/ape::Nedge(temp_tree), ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <- matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <- matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <- matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <- matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <- matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <- matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1, genotype2, genotype3, genotype4, genotype5, genotype6, genotype7, genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[ , 1] <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[ , 2] <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")

  temp_result <- prepare_grouped_genotype(temp_geno, temp_lookup)

  expect_equal(temp_result$genotype, temp_geno[ , c(1, 2, 5, 6, 7, 8)])
  expect_equal(temp_result$gene_snp_lookup, temp_lookup[c(1, 2, 5, 6, 7, 8), ])
  expect_equal(as.numeric(unname(temp_result$snps_per_gene)), c(1, 1, 1, 1, 2))
  expect_equal(temp_result$unique_genes, c("GENE1", "GENE2", "GENE5", "GENE6", "GENE7"))
})

# test group_genotypes ---------------------------------------------------------
test_that("group_genotypes does X given Y", {
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <- rep(sum(temp_tree$edge.length)/ape::Nedge(temp_tree), ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  set.seed(1)
  temp_pheno <- as.matrix(phytools::fastBM(temp_tree))
  row.names(temp_pheno) <- temp_tree$tip.label
  colnames(temp_pheno) <- "growth"
  temp_continuous <- "continuous"

  genotype1 <- matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <- matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <- matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <- matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <- matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <- matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <- matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1, genotype2, genotype3, genotype4, genotype5, genotype6, genotype7, genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[ , 1] <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[ , 2] <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")

  temp_group_logical <- TRUE

  geno <- prepare_genotype(temp_group_logical, temp_geno, temp_tree, temp_lookup)
  genotype <- geno$genotype
  AR <- prepare_ancestral_reconstructions(temp_tree, temp_pheno, genotype, temp_continuous)
  geno_trans_concomitant <- AR$geno_trans # Include all transition edges (WT -> mutant and mutant -> WT). For discrete concomitant and continuous tests.
  geno_trans_original    <- prepare_genotype_transitions_for_original_discrete_test(genotype, geno_trans_concomitant) # Keep only WT -> mutant transitions.

  temp_results <- group_genotypes(temp_tree, genotype, AR$geno_recon_and_conf, geno_trans_concomitant, geno_trans_original, geno$gene_snp_lookup, geno$unique_genes)
  expect_equal(length(temp_results$geno_recon_ordered_by_edges), ncol(temp_results$genotype))
  expect_identical(row.names(temp_results$genotype), temp_tree$tip.label)
  expect_equal(length(temp_results$geno_trans_concomitant[[1]]$transition), ape::Nedge(temp_tree))
  expect_equal(length(temp_results$geno_trans_orig[[1]]$transition), ape::Nedge(temp_tree))
  expect_equal(length(temp_results$geno_trans_concomitant[[1]]$trans_dir), ape::Nedge(temp_tree))
  expect_equal(length(temp_results$geno_trans_orig[[1]]$trans_dir), ape::Nedge(temp_tree))
})
