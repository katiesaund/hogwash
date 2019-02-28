library(microbeGWAS)

context("first test")
test_that("check_tree_is_valid works", {
  for (i in 2:10){
    expect_true(check_tree_is_valid(rtree(i)))
  }
  invalid_tree <- rtree(10)
  for (j in 1:Nedge(invalid_tree)){
    if (invalid_tree$edge[j, 2] == 1){
      invalid_tree$edge[j, 2] <- 10000
      break
    }
  }

  expect_that(check_tree_is_valid(invalid_tree), throws_error())
})


test_that("ancestral reconstructions have correct length", {
  tree <- rtree(9) # TO DO ADD NODE VALS
  test_mat <- matrix(rep(1, 0), nrow = Ntip(tree), ncol = 9)
  dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "discrete", 0.70) # add a continuous one
  dummy_geno <-  ancestral_reconstruction_by_ML(tree, test_mat, 2, "discrete", 0.70)
  expect_identical(length(dummy_pheno$tip_and_node_recon), Nedge(tree))
  expect_identical(length(dummy_geno$tip_and_node_recon), Nedge(tree))
})
