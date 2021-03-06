context("Convergence") #---------------------------------------------------#
test_that("calculate_continuous_convergence gives expected results for valid continuous input", {
  # Set up
  set.seed(1)
  num_tip <- 5
  tree <- ape::rcoal(num_tip)
  tree$node.label <- rep(100, ape::Nnode(tree))
  cont_pheno <- rnorm(num_tip, mean = 0, sd = 2)
  tree_recon <-  suppressWarnings(ape::ace(x = cont_pheno,
                                           phy = tree,
                                           type = "continuous",
                                           method = "REML",
                                           model = "BM",
                                           marginal = FALSE))
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

  output <- calculate_continuous_convergence(pheno_recon_mat, high_conf)
  delta_pheno <- abs(pheno_recon_mat[, 1] - pheno_recon_mat[, 2])
  scaled_delta_pheno <- scales::rescale(delta_pheno)

  # What I expect
  beta_pheno <- sum(scaled_delta_pheno) # all are high confidence
  beta_geno <- 4
  intersection <-
    sum(scaled_delta_pheno * high_conf$genotype_transition[[1]]$transition)
  # all are high confidence

  union <- beta_pheno + beta_geno - intersection
  epsilon <- intersection / union

  # Test
  expect_equal(output$num_hi_conf_edges, rep(num_edge, num_geno))
  expect_equal(output$intersection, rep(intersection, num_geno))
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
  new_output <- calculate_continuous_convergence(pheno_recon_mat, high_conf)

  # What I expect
  beta_pheno <- sum(scaled_delta_pheno) # all are high confidence
  beta_geno <- 4
  intersection <-
    sum(scaled_delta_pheno * high_conf$genotype_transition[[1]]$transition)
  # all are high confidence

  union <- beta_pheno + beta_geno - intersection
  epsilon <- intersection / union

  # Test
  expect_equal(new_output$num_hi_conf_edges, rep(num_edge, num_geno))
  expect_equal(new_output$intersection, rep(intersection, num_geno))
  expect_equal(new_output$pheno_beta, rep(beta_pheno, num_geno))
  expect_equal(new_output$geno_beta, rep(beta_geno, num_geno))
  expect_equal(new_output$epsilon, rep(epsilon, num_geno))
})

test_that("calculate_phyc_convergence works as expected", {
  # Set up
  set.seed(1)
  num_tip <- 5
  tree <- ape::rcoal(num_tip)
  tree$node.label <- rep(100, ape::Nnode(tree))
  bin_pheno <- c(1, 0, 1, 0, 1)
  tree_recon_tip_and_node <- c(bin_pheno, 1, 0, 1, 1)
  pheno_recon_mat <- convert_to_edge_mat(tr = tree, tree_recon_tip_and_node)
  high_conf <- geno_transition_list <- NULL
  num_edge <- ape::Nedge(tree)
  num_geno <- 5
  high_conf$tr_and_pheno_hi_conf <- rep(TRUE, num_edge)
  high_conf$genotype_transition <- rep(list(list(0)), num_geno)
  high_conf$high_conf_ordered_by_edges <- rep(list(0), num_geno)

  for (i in 1:num_geno) {
    high_conf$genotype_transition[[i]]$transition <- rep(c(0, 1), num_edge / 2)
    high_conf$high_conf_ordered_by_edges[[i]] <- rep(1, num_edge)
    geno_transition_list[[i]] <- high_conf$genotype_transition[[i]]$transition
  }

  pheno_recon_vec <- c(1, 1, 0, 0, 0, 1, 1, 0)

  output <-
    calculate_phyc_convergence(geno_trans_edge_list = geno_transition_list,
                               pheno_recon_vec = pheno_recon_vec,
                               high_conf = high_conf)

  # What I expect
  beta_pheno <- sum(pheno_recon_vec) # all are high confidence
  beta_geno <- sum(geno_transition_list[[1]])
  intersection <- sum(geno_transition_list[[1]] == 1 & pheno_recon_vec == 1)
  # all are high confidence
  epsilon <- 2 * intersection / (beta_pheno + beta_geno)

  # Test
  expect_equal(output$num_hi_conf_edges, rep(num_edge, num_geno))
  expect_equal(output$intersection, rep(intersection, num_geno))
  expect_equal(output$pheno_beta, beta_pheno)
  expect_equal(output$geno_beta, rep(beta_geno, num_geno))
  expect_equal(output$epsilon, rep(epsilon, num_geno))
})


test_that("calculate_synchronous_convergence works as expected", {
  # Set up
  set.seed(1)
  num_tip <- 5
  tree <- ape::rcoal(num_tip)
  tree$node.label <- rep(100, ape::Nnode(tree))
  bin_pheno <- c(1, 0, 1, 0, 1)
  tree_recon_tip_and_node <- c(bin_pheno, 1, 0, 1, 1)
  pheno_recon_mat <- convert_to_edge_mat(tr = tree, tree_recon_tip_and_node)
  high_conf <- geno_transition_list <- NULL
  num_edge <- ape::Nedge(tree)
  num_geno <- 5
  high_conf$tr_and_pheno_hi_conf <- rep(TRUE, num_edge)
  high_conf$genotype_transition <- rep(list(list(0)), num_geno)
  high_conf$high_conf_ordered_by_edges <- rep(list(0), num_geno)

  for (i in 1:num_geno) {
    high_conf$genotype_transition[[i]]$transition <- rep(c(0, 1), num_edge / 2)
    high_conf$high_conf_ordered_by_edges[[i]] <- rep(1, num_edge)
    geno_transition_list[[i]] <- high_conf$genotype_transition[[i]]$transition
  }

  pheno_trans_vec <- NULL
  pheno_trans_vec$transition <- c(1, 1, 0, 0, 0, 1, 1, 0)

  output <-
    calculate_synchronous_convergence(
      geno_trans_edge_list = geno_transition_list,
      pheno_trans_vec = pheno_trans_vec,
      high_conf = high_conf)

  # What I expect
  beta_pheno <- sum(pheno_trans_vec$transition) # all are high confidence
  beta_geno <- sum(geno_transition_list[[1]])
  intersection <-
    sum(geno_transition_list[[1]] == 1 & pheno_trans_vec$transition == 1)
  # all are high confidence
  epsilon <- 2 * intersection / (beta_pheno + beta_geno)

  # Test
  expect_equal(output$num_hi_conf_edges, rep(num_edge, num_geno))
  expect_equal(output$intersection, rep(intersection, num_geno))
  expect_equal(output$pheno_beta, beta_pheno)
  expect_equal(output$geno_beta, rep(beta_geno, num_geno))
  expect_equal(output$epsilon, rep(epsilon, num_geno))
})
