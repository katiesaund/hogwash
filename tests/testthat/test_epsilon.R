test_that("calculate_continuous_gamma gives expected results for valid continuous input", {
  # Set up
  set.seed(1)
  num_tip <- 5
  tree <- ape::rcoal(num_tip)
  tree$node.label <- rep(100, ape::Nnode(tree))
  cont_pheno <- rnorm(num_tip, mean = 0, sd = 2)
  tree_recon <- ape::ace(x = cont_pheno, phy = tree, type = "continuous", method = "ML", model = "BM", marginal = FALSE)
  tree_recon_tip_and_node <- c(cont_pheno, tree_recon$ace)
  pheno_recon_mat <- convert_to_edge_mat(tr = tree, tree_recon_tip_and_node)
  high_conf <- NULL
  num_edge <- ape::Nedge(tree)
  num_geno <- 5
  high_conf$tr_and_pheno_hi_conf <- rep(TRUE, num_edge)
  high_conf$genotype_transition <- rep(list(list(0)), num_geno)
  high_conf$high_conf_ordered_by_edges <- rep(list(0), num_geno)

  for (i in 1:num_geno) {
    high_conf$genotype_transition[[i]]$transition <- rep(c(0, 1), num_edge / 2)
    high_conf$high_conf_ordered_by_edges[[i]] <- rep(1, num_edge)
  }

  output <- calculate_continuous_gamma(pheno_recon_mat, high_conf)
  # Test
  gamma <- 1
  beta_pheno <- 3
  beta_geno <- 4
  epsilon <- 2 * gamma / (beta_pheno + beta_geno)
  expect_equal(output$num_hi_conf_edges, rep(num_edge, num_geno))
  expect_equal(output$gamma_count, rep(gamma, num_geno))
  expect_equal(output$pheno_beta, rep(beta_pheno, num_geno))
  expect_equal(output$geno_beta, rep(beta_geno, num_geno))
  expect_equal(output$epsilon, rep(epsilon, num_geno))

  # Gives rather low epsilon because the two distributions overlap a lot and
  #   transition edges are left skewed.

  # Repeat for more right skewed genotype transition edges and left skewed
  #   non-transition edges
  for (i in 1:num_geno) {
    high_conf$genotype_transition[[i]]$transition <- c(1, 1, 0, 0, 0, 0, 1, 1)
  }
  new_output <- calculate_continuous_gamma(pheno_recon_mat, high_conf)
  gamma <- 4
  beta_pheno <- 6
  beta_geno <- 4
  epsilon <- 2 * gamma / (beta_pheno + beta_geno)
  expect_equal(new_output$num_hi_conf_edges, rep(num_edge, num_geno))
  expect_equal(new_output$gamma_count, rep(gamma, num_geno))
  expect_equal(new_output$pheno_beta, rep(beta_pheno, num_geno))
  expect_equal(new_output$geno_beta, rep(beta_geno, num_geno))
  expect_equal(new_output$epsilon, rep(epsilon, num_geno))
})
