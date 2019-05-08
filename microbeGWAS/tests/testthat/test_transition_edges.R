library(microbeGWAS)
context("Transition edges") #----------------------------------------------#
# TODO write tests for all functions

# TODO test identify_transition_edges

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

# keep_at_least_two_high_conf_trans_edges <- function(genotype_transition, genotype_confidence){

#   has_at_least_two_high_confidence_transition_edges <- rep(FALSE, length(genotype_transition))
#   for (p in 1:length(genotype_transition)){
#     if (sum(genotype_transition[[p]]$transition * genotype_confidence[[p]]) > 1){
#       has_at_least_two_high_confidence_transition_edges[p] <- TRUE
#     }
#   }
#
#   # Check and return output ----------------------------------------------------
#   return(has_at_least_two_high_confidence_transition_edges)
# } # end keep_at_least_two_high_conf_trans_edges()
#
# TODO test keep_hits_with_more_change_on_trans_edges
# keep_hits_with_more_change_on_trans_edges <- function(results, pvals, a){
#   # TODO
#   # 1) Update description of inputs
#   # 2) check inputs.
#   # 3) Update description of output.
#   # 4) Check output.
#   #
#   # Function description -------------------------------------------------------
#   # SUBSET SIGNIFICANT HITS WHERE THE MEDIAN(DELTA PHENOTYPE) ON TRANSITION EDGES IS > MEDIAN(DELTA PHENOTYPE) ON ALL EDGES
#   #
#   # Inputs:
#   # results.
#   # pvals.
#   # a.      Number. Alpha (significance threshold).
#   # Output:
#   # has_at_least_two_high_confidence_transition_edges.
#   #
#   # Check inputs ---------------------------------------------------------------
#   check_if_alpha_valid(a)
#
#   # Function -------------------------------------------------------------------
#   temp <- pvals$hit_pvals[(results$trans_median > results$all_edges_median), , drop = FALSE]
#   hits <- temp[temp[ , 1] < a, , drop = FALSE]
#
#   # Check and return output ----------------------------------------------------
#   return(hits)
# } # end keep_hits_with_more_change_on_trans_edges()
