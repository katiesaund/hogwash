context("Group genotypes") #---------------------------------------------------#

# test prepare_genotype --------------------------------------------------------
test_that("prepare_genotype gives expected genotype for a pre-AR grouped input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <-
    rep(sum(temp_tree$edge.length) / ape::Nedge(temp_tree),
        ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <-
    matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <-
    matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <-
    matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <-
    matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <-
    matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <-
    matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1,
                     genotype2,
                     genotype3,
                     genotype4,
                     genotype5,
                     genotype6,
                     genotype7,
                     genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[, 1] <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[, 2] <-
    c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")

  temp_logical <- TRUE
  temp_result <- prepare_genotype(temp_logical,
                                  temp_geno,
                                  temp_tree,
                                  temp_lookup,
                                  "pre-ar")

  # Test
  # Expected output matrix
  expected_geno <- temp_geno[, c(1, 5, 6, 7)]
  expected_geno[, 4] <- rowSums(temp_geno[, 7:8])
  colnames(expected_geno) <- paste0("GENE", c(1, 5, 6, 7))
  expect_equal(temp_result$genotype, expected_geno)
  expect_equal(temp_result$gene_snp_lookup, temp_lookup[c(1, 2, 5, 6, 7, 8), ])
  expect_equal(as.numeric(unname(temp_result$snps_per_gene)), c(1, 1, 1, 1, 2))
  expect_equal(temp_result$unique_genes,
               c("GENE1", "GENE2", "GENE5", "GENE6", "GENE7"))
})

test_that("prepare_genotype gives expected genotype for a not grouped input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <-
    rep(sum(temp_tree$edge.length) / ape::Nedge(temp_tree),
        ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <-
    matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <-
    matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <-
    matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <-
    matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <-
    matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <-
    matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1,
                     genotype2,
                     genotype3,
                     genotype4,
                     genotype5,
                     genotype6,
                     genotype7,
                     genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[, 1] <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[, 2] <-
    c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")
  temp_logical <- FALSE
  temp_result <- prepare_genotype(temp_logical,
                                  temp_geno,
                                  temp_tree,
                                  temp_lookup,
                                  "pre-ar")

  # Test
  expect_null(temp_result$snps_per_gene)
  expect_equal(temp_result$genotype, temp_geno[, c(1, 5, 6)])
  expect_equal(temp_result$no_convergence_genotypes,
               c("SNP2", "SNP3", "SNP4", "SNP7", "SNP8"))
})

# test prepare_ungrouped_genotype ----------------------------------------------
test_that("prepare_ungrouped_genotype works for valid input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <-
    rep(sum(temp_tree$edge.length) / ape::Nedge(temp_tree),
        ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <-
    matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <-
    matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <-
    matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <-
    matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <-
    matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <-
    matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1,
                     genotype2,
                     genotype3,
                     genotype4,
                     genotype5,
                     genotype6,
                     genotype7,
                     genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_result <- prepare_ungrouped_genotype(temp_geno, temp_tree)

  # Test
  expect_null(temp_result$snps_per_gene)
  expect_equal(temp_result$genotype, temp_geno[, c(1, 5, 6)])
  expect_equal(temp_result$no_convergence_genotypes,
               c("SNP2", "SNP3", "SNP4", "SNP7", "SNP8"))
})

# test prepare_grouped_genotype_pre_ar -----------------------------------------
test_that("prepare_grouped_genotype_pre_ar works for valid inputs", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <-
    rep(sum(temp_tree$edge.length) / ape::Nedge(temp_tree),
        ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <-
    matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <-
    matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <-
    matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <-
    matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <-
    matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <-
    matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1,
                     genotype2,
                     genotype3,
                     genotype4,
                     genotype5,
                     genotype6,
                     genotype7,
                     genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[, 1] <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[, 2] <-
    c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")

  temp_result <-
    prepare_grouped_genotype_pre_ar(temp_geno, temp_lookup, temp_tree)

  # Test
  expected_geno <- temp_geno[, c(1, 5, 6, 7)]
  expected_geno[, 4] <- rowSums(temp_geno[,7:8])
  colnames(expected_geno) <- paste0("GENE", c(1, 5, 6, 7))
  expect_equal(temp_result$genotype, expected_geno)
  expect_equal(temp_result$gene_snp_lookup, temp_lookup[c(1, 2, 5, 6, 7, 8), ])
  expect_equal(as.numeric(unname(temp_result$snps_per_gene)), c(1, 1, 1, 1, 2))
  expect_equal(temp_result$unique_genes,
               c("GENE1", "GENE2", "GENE5", "GENE6", "GENE7"))
})

# test build_gene_genotype_from_snps_pre_ar ------------------------------------
test_that("build_gene_genotype_from_snps_pre_ar works for valid inputs", {
  # Set up
  ntip <- 4
  temp_geno <- matrix(c(0, 1),
                      ncol = 4,
                      nrow = ntip)
  temp_geno[2, 1] <- temp_geno[2, 3] <- temp_geno[4, 2] <- 0
  temp_geno[1, 3] <- 1

  colnames(temp_geno) <- c("a", "b", "d", "h")
  row.names(temp_geno) <- c("sample1", "sample2", "sample3", "sample4")
  temp_key <- matrix(c(letters[1:8],
                     c(rep("One", 3), rep("Two", 3), rep("Three", 2))),
                     ncol = 2,
                     nrow = 8)
  colnames(temp_key) <- c("snp", "gene")
  temp_key <- temp_key[temp_key[, 1] %in% colnames(temp_geno),, drop = FALSE]
  temp_tree <- ape::rcoal(ntip)
  temp_tree$node.label <- c(100, 100, 100)
  temp_tree$tip.label <- row.names(temp_geno)
  results <- build_gene_genotype_from_snps_pre_ar(temp_geno,
                                                  temp_key,
                                                  temp_tree)

  expected_results <- temp_geno[, c(1, 3:4)]
  expected_results[2, 1] <- 1
  colnames(expected_results) <- c("One", "Two", "Three")

  # Test
  expect_identical(results, expected_results)
})

test_that("build_gene_genotype_from_snps_pre_ar gives error for invalid inputs", {
  # Set up
  temp_geno <- matrix(c(0, 1),
                      ncol = 4,
                      nrow = 3)
  colnames(temp_geno) <- c("a", "b", "d", "h")
  row.names(temp_geno) <- c("sample1", "sample2", "sample3")
  temp_key <- matrix(c(letters[1:8],
                       c(rep("One", 3), rep("Two", 3), rep("Three", 2))),
                     ncol = 2,
                     nrow = 8)
  colnames(temp_key) <- c("snp", "gene")

  # Test
  expect_error(build_gene_genotype_from_snps_pre_ar(temp_geno, temp_key))
})

# Below: post-AR unit tests

context("Group genotypes") #---------------------------------------------------#
# test build_gene_trans_from_snp_trans_post_ar
test_that("build_gene_trans_from_snp_trans_post_ar works for valid input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <-
    rep(sum(temp_tree$edge.length) / ape::Nedge(temp_tree),
        ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <-
    matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <-
    matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <-
    matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <-
    matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype <-
    cbind(genotype1, genotype2, genotype5, genotype6, genotype7, genotype8)
  row.names(genotype) <- temp_tree$tip.label
  colnames(genotype) <- c("SNP1", "SNP2", "SNP5", "SNP6", "SNP7", "SNP8")

  set.seed(1)
  continuous_phenotype <- as.matrix(phytools::fastBM(temp_tree))
  row.names(continuous_phenotype) <- temp_tree$tip.label
  colnames(continuous_phenotype) <- "growth"


  snp_gene_key <- matrix(NA, nrow = ncol(genotype), ncol = 2)
  colnames(snp_gene_key) <- c("SNP", "GENE")
  snp_gene_key[, 1] <- c("SNP1", "SNP2", "SNP5", "SNP6", "SNP7", "SNP8")
  snp_gene_key[, 2] <- c("GENE1", "GENE2", "GENE5", "GENE6", "GENE7", "GENE7")

  AR <- prepare_ancestral_reconstructions(temp_tree,
                                          continuous_phenotype,
                                          genotype,
                                          "continuous")

  temp_results <- build_gene_trans_from_snp_trans_post_ar(temp_tree,
                                                  genotype,
                                                  AR$geno_trans,
                                                  snp_gene_key)

  # Test
  expect_true(length(temp_results) == length(unique(snp_gene_key[, 2])))
  expect_equal(temp_results[[1]]$transition, AR$geno_trans[[1]]$transition)
  expect_equal(temp_results[[2]]$transition, AR$geno_trans[[2]]$transition)
  expect_equal(temp_results[[3]]$transition, AR$geno_trans[[3]]$transition)
  expect_equal(temp_results[[4]]$transition, AR$geno_trans[[4]]$transition)
  expect_equal(temp_results[[5]]$transition,
               AR$geno_trans[[5]]$transition + AR$geno_trans[[6]]$transition)
  expect_equal(length(temp_results[[1]]$transition), ape::Nedge(temp_tree))
  expect_equal(length(temp_results[[1]]$trans_dir), ape::Nedge(temp_tree))
})

# test prepare_genotype --------------------------------------------------------
test_that("prepare_genotype gives expected genotype for a post-ar grouped input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <-
    rep(sum(temp_tree$edge.length) / ape::Nedge(temp_tree),
        ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <-
    matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <-
    matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <-
    matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <-
    matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <-
    matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <-
    matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1,
                     genotype2,
                     genotype3,
                     genotype4,
                     genotype5,
                     genotype6,
                     genotype7,
                     genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[, 1] <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[, 2] <-
    c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")

  temp_logical <- TRUE
  temp_result <- prepare_genotype(temp_logical,
                                  temp_geno,
                                  temp_tree,
                                  temp_lookup,
                                  "post-ar")

  # Test
  expect_equal(temp_result$genotype, temp_geno[, c(1, 2, 5, 6, 7, 8)])
  expect_equal(temp_result$gene_snp_lookup, temp_lookup[c(1, 2, 5, 6, 7, 8), ])
  expect_equal(as.numeric(unname(temp_result$snps_per_gene)), c(1, 1, 1, 1, 2))
  expect_equal(temp_result$unique_genes,
               c("GENE1", "GENE2", "GENE5", "GENE6", "GENE7"))

})
test_that("prepare_genotype gives expected genotype for a not grouped input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <-
    rep(sum(temp_tree$edge.length) / ape::Nedge(temp_tree),
        ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <-
    matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <-
    matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <-
    matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <-
    matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <-
    matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <-
    matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1,
                     genotype2,
                     genotype3,
                     genotype4,
                     genotype5,
                     genotype6,
                     genotype7,
                     genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[, 1] <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[, 2] <-
    c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")
  temp_logical <- FALSE
  temp_result <- prepare_genotype(temp_logical,
                                  temp_geno,
                                  temp_tree,
                                  temp_lookup,
                                  "post-ar")

  # Test
  expect_null(temp_result$snps_per_gene)
  expect_equal(temp_result$genotype, temp_geno[, c(1, 5, 6)])
  expect_equal(temp_result$no_convergence_genotypes,
               c("SNP2", "SNP3", "SNP4", "SNP7", "SNP8"))
})

# test prepare_grouped_genotype_post_ar ----------------------------------------
test_that("prepare_grouped_genotype_post_ar works for valid inputs", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <-
    rep(sum(temp_tree$edge.length) / ape::Nedge(temp_tree),
        ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  genotype1 <-
    matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <-
    matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <-
    matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <-
    matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <-
    matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <-
    matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1,
                     genotype2,
                     genotype3,
                     genotype4,
                     genotype5,
                     genotype6,
                     genotype7,
                     genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[, 1] <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[, 2] <-
    c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")

  temp_result <- prepare_grouped_genotype_post_ar(temp_geno, temp_lookup)

  # Test
  expect_equal(temp_result$genotype, temp_geno[, c(1, 2, 5, 6, 7, 8)])
  expect_equal(temp_result$gene_snp_lookup, temp_lookup[c(1, 2, 5, 6, 7, 8), ])
  expect_equal(as.numeric(unname(temp_result$snps_per_gene)), c(1, 1, 1, 1, 2))
  expect_equal(temp_result$unique_genes,
               c("GENE1", "GENE2", "GENE5", "GENE6", "GENE7"))
})

# test group_genotypes_post_ar -------------------------------------------------
test_that("group_genotypes_post_ar works for valid inputs", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(7)
  temp_tree$edge.length <-
    rep(sum(temp_tree$edge.length) / ape::Nedge(temp_tree),
        ape::Nedge(temp_tree))
  temp_tree$node.label <- c(100, 100, 50, 100, 100, 100) # 1 low confidence edge

  set.seed(1)
  temp_pheno <- as.matrix(phytools::fastBM(temp_tree))
  row.names(temp_pheno) <- temp_tree$tip.label
  colnames(temp_pheno) <- "growth"
  temp_continuous <- "continuous"

  genotype1 <-
    matrix(c(0, 1, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype2 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype3 <-
    matrix(c(0, 0, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype4 <-
    matrix(c(1, 1, 1, 1, 1, 1, 1), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype5 <-
    matrix(c(0, 1, 1, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype6 <-
    matrix(c(0, 1, 1, 1, 1, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype7 <-
    matrix(c(0, 1, 0, 0, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  genotype8 <-
    matrix(c(0, 0, 0, 1, 0, 0, 0), nrow = ape::Ntip(temp_tree), ncol = 1)
  temp_geno <- cbind(genotype1,
                     genotype2,
                     genotype3,
                     genotype4,
                     genotype5,
                     genotype6,
                     genotype7,
                     genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

  temp_lookup <- matrix(NA, nrow = ncol(temp_geno), ncol = 2)
  colnames(temp_lookup) <- c("SNP", "GENE")
  temp_lookup[, 1] <-
    c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_lookup[, 2] <-
    c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE7")

  temp_group_logical <- TRUE

  geno <- prepare_genotype(temp_group_logical,
                           temp_geno,
                           temp_tree,
                           temp_lookup,
                           "post-ar")
  genotype <- geno$genotype
  AR <- prepare_ancestral_reconstructions(temp_tree,
                                          temp_pheno,
                                          genotype,
                                          temp_continuous)

  # Include all transition edges (WT -> mutant and mutant -> WT).
  #  For synchronous and continuous tests.
  geno_trans_synchronous <- AR$geno_trans

  # Keep only WT -> mutant transitions.
  geno_trans_phyc <- prep_geno_trans_for_phyc(genotype, geno_trans_synchronous)

  temp_results <- group_genotypes_post_ar(temp_tree,
                                          genotype,
                                          AR$geno_recon_and_conf,
                                          geno_trans_synchronous,
                                          geno_trans_phyc,
                                          geno$gene_snp_lookup,
                                          geno$unique_genes)

  # Test
  expect_equal(length(temp_results$geno_recon_ordered_by_edges),
               ncol(temp_results$genotype))
  expect_identical(row.names(temp_results$genotype), temp_tree$tip.label)
  expect_equal(length(temp_results$geno_trans_synchronous[[1]]$transition),
               ape::Nedge(temp_tree))
  expect_equal(length(temp_results$geno_trans_phyc[[1]]$transition),
               ape::Nedge(temp_tree))
  expect_equal(length(temp_results$geno_trans_synchronous[[1]]$trans_dir),
               ape::Nedge(temp_tree))
  expect_equal(length(temp_results$geno_trans_phyc[[1]]$trans_dir),
               ape::Nedge(temp_tree))
})

# test build_gene_genotype_from_snps_post_ar
test_that("build_gene_genotype_from_snps_post_ar works for valid inputs", {
  # Set up
  temp_geno <- matrix(c(0, 1),
                      ncol = 4,
                      nrow = 3)
  colnames(temp_geno) <- c("a", "b", "d", "h")
  row.names(temp_geno) <- c("sample1", "sample2", "sample3")
  temp_key <- matrix(c(letters[1:8],
                       c(rep("One", 3), rep("Two", 3), rep("Three", 2))),
                     ncol = 2,
                     nrow = 8)
  colnames(temp_key) <- c("snp", "gene")
  temp_key <- temp_key[temp_key[, 1] %in% colnames(temp_geno),, drop = FALSE]
  results <- build_gene_genotype_from_snps_post_ar(temp_geno, temp_key)

  expected_results <- matrix(c(1, 1, 1, 0, 1, 0, 1, 0, 1), nrow = 3, ncol = 3)
  colnames(expected_results) <- c("One", "Two", "Three")
  row.names(expected_results) <- c("sample1", "sample2", "sample3")

  # Test
  expect_identical(results, expected_results)
})

test_that("build_gene_genotype_from_snps_post_ar gives error for invalid inputs", {
  # Set up
  temp_geno <- matrix(c(0, 1),
                      ncol = 4,
                      nrow = 3)
  colnames(temp_geno) <- c("a", "b", "d", "h")
  row.names(temp_geno) <- c("sample1", "sample2", "sample3")
  temp_key <- matrix(c(letters[1:8],
                       c(rep("One", 3), rep("Two", 3), rep("Three", 2))),
                     ncol = 2,
                     nrow = 8)
  colnames(temp_key) <- c("snp", "gene")

  # Test
  expect_error(build_gene_genotype_from_snps_post_ar(temp_geno, temp_key))
})
