library(microbeGWAS)
context("Transition edges") #----------------------------------------------#
# test identify_transition_edges

test_that("identify_transition_edges returns correct transition vector and trans_dir vector for discrete phenotype", {
  set.seed(1)
  temp_tree <- rtree(5)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_pheno <- matrix(c(1, 0, 1, 0, 1), ncol = 1)
  row.names(temp_pheno) <- temp_tree$tip.label
  set.seed(1)
  temp_recon <- ancestral_reconstruction_by_ML(temp_tree, temp_pheno, 1, "discrete")
  temp_results <- identify_transition_edges(temp_tree, temp_pheno, 1, temp_recon$node_anc_rec, "discrete")
  expect_equivalent(temp_results$transition, c(1, 1, 0, 0, 0, 1, 0, 1))
  expect_equivalent(temp_results$trans_dir, c(-1, 1, 0, 0, 0, -1, 0, 1))
})

test_that("identify_transition_edges returns correct transition vector and trans_dir vector for discrete genotype", {
  set.seed(1)
  temp_tree <- rtree(5)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_geno <- cbind(c(1, 0, 1, 0, 1), c(1, 0, 0, 0, 0))
  row.names(temp_geno) <- temp_tree$tip.label
  temp_recon <- temp_results <- rep(list(0), ncol(temp_geno))
  set.seed(1)
  for (i in 1:2){
    temp_recon[[i]] <- ancestral_reconstruction_by_ML(temp_tree, temp_geno, i, "discrete")
    temp_results[[i]] <- identify_transition_edges(temp_tree, temp_geno, i, temp_recon[[i]]$node_anc_rec, "discrete")
  }
  expect_equivalent(temp_results[[1]]$transition, c(1,  1, 0, 0, 0,  1, 0, 1))
  expect_equivalent(temp_results[[1]]$trans_dir,  c(-1, 1, 0, 0, 0, -1, 0, 1))
  expect_equivalent(temp_results[[2]]$transition, c(0,  1, 0, 0, 0,  0, 0, 0))
  expect_equivalent(temp_results[[2]]$trans_dir,  c(0,  1, 0, 0, 0,  0, 0, 0))
})


test_that("identify_transition_edges returns correct transition vector and trans_dir vector for continuous phenotype", {
  set.seed(1)
  temp_tree <- rtree(5)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_pheno <- matrix(rnorm(5), ncol = 1)
  row.names(temp_pheno) <- temp_tree$tip.label
  set.seed(1)
  temp_recon <- ancestral_reconstruction_by_ML(temp_tree, temp_pheno, 1, "continuous")
  temp_results <- identify_transition_edges(temp_tree, temp_pheno, 1, temp_recon$node_anc_rec, "continuous")
  expect_equivalent(temp_results$transition, NA)
  expect_equivalent(temp_results$trans_dir, c(1, -1,  1, -1,  1, -1, -1, -1))
})

# test keep_at_least_two_high_conf_trans_edges
test_that("keep_at_least_two_high_conf_trans_edges returns a vector c(TRUE, FALSE) for this test genotype", {
  set.seed(1)
  temp_tree <- rtree(5)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_geno <- cbind(c(1, 0, 1, 0, 1), c(1, 0, 0, 0, 0))
  row.names(temp_geno) <- temp_tree$tip.label
  fake_confidence <- temp_recon <- temp_results <- rep(list(0), ncol(temp_geno))
  set.seed(1)
  for (i in 1:2){
    temp_recon[[i]] <- ancestral_reconstruction_by_ML(temp_tree, temp_geno, i, "discrete")
    temp_results[[i]] <- identify_transition_edges(temp_tree, temp_geno, i, temp_recon[[i]]$node_anc_rec, "discrete")
    fake_confidence[[i]] <- rep(1, Nedge(temp_tree))
  }
  expect_equivalent(keep_at_least_two_high_conf_trans_edges(temp_results, fake_confidence), c(TRUE, FALSE))
  for (i in 1:2){
    fake_confidence[[i]] <- rep(0, Nedge(temp_tree))
  }
  expect_equivalent(keep_at_least_two_high_conf_trans_edges(temp_results, fake_confidence), c(FALSE, FALSE))
})

# test keep_hits_with_more_change_on_trans_edges
test_that("keep_hits_with_more_change_on_trans_edges does X given Y", {
  num_genotypes <- 5
  temp_fdr <- 0.05
  temp_results <- temp_pvals <- NULL
  temp_results$trans_median <- rep(10, num_genotypes)
  temp_results$all_edges_median <- c(8:12)
  temp_pvals$hit_pvals <- as.data.frame(matrix(0.01, ncol = 1, nrow = num_genotypes))
  colnames(temp_pvals$hit_pvals) <- "fdr_corrected_pvals"
  row.names(temp_pvals$hit_pvals) <- letters[1:num_genotypes]
  results <- keep_hits_with_more_change_on_trans_edges(temp_results, temp_pvals, temp_fdr)
  expect_equivalent(row.names(results), c("a", "b"))
})
