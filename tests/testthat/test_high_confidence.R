context("High confidence") #----------------------------------------------#
# TODO write unit tests for assign_high_confidence_to_transition_edges

# test discretize_conf_with_cutoff
test_that("discretize_conf_with_cutoff should give this expected result", {
  temp_vector <- c(1:100)
  temp_threshold <- 50
  expect_equal(discretize_conf_with_cutoff(temp_vector,
                                                     temp_threshold),
               c(rep(0, 49), rep(1, 51)))
})

test_that("discretize_conf_with_cutoff should give this expected result -
          now with fractions", {
  temp_vector <- seq(from = 0, to = 1, by = 0.1)
  temp_threshold <- 0.5
  expect_equal(discretize_conf_with_cutoff(temp_vector,
                                                     temp_threshold),
               c(rep(0, 5), rep(1, 6)))
})

# test report_num_high_confidence_trans_edge
test_that("report_num_high_confidence_trans_edge returns expected outcome for
          this test set", {
  fake_geno_names <- letters[1:3]
  fake_trans <- fake_hi_conf_edges <- rep(list(NULL), length(fake_geno_names))
  for (i in 1:length(fake_geno_names)) {
    fake_trans[[i]]$transition <- c(1, 1, 1)
    fake_hi_conf_edges[[i]] <- c(1, 0, 1)
  }
  fake_trans[[1]]$transition <- c(1, 1, 1)
  fake_trans[[2]]$transition <- c(0, 0, 0)
  fake_trans[[3]]$transition <- c(1, 1, 0)
  expected_result <- c(2, 0, 1)
  names(expected_result) <- fake_geno_names
  expect_equal(report_num_high_confidence_trans_edge(fake_trans,
                                                     fake_hi_conf_edges,
                                                     fake_geno_names),
               expected_result)
})

# test assign_high_confidence_to_transition_edges
test_that("assign_high_confidence_to_transition_edges returns the edges that are
          high confidence transition edges for this tree", {
  set.seed(1)
  num_samples <- 5
  temp_tree <- ape::rtree(num_samples)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_trans <-
    temp_confidence <-
    expected_result <-
    rep(list(NULL), num_samples)
  for (i in 1:num_samples) {
    temp_trans[[i]]$transition <- c(0, 0, 0, 1, 0, 1, 1, 0)
    temp_confidence[[i]] <- c(0, 0, 0, 0, 0, 1, 1, 0)
    expected_result[[i]] <- c(0, 0, 0, 0, 0, 1, 1, 0)
  }

  # FYI this geno does not match up with the fake transitions I made up
  temp_geno <- matrix(1, ncol = num_samples, nrow = num_samples)

  temp_tree$tip.label <- row.names(temp_geno) <- letters[1:num_samples]
  foo <- assign_high_confidence_to_transition_edges(temp_tree,
                                                    temp_confidence,
                                                    temp_trans,
                                                    temp_geno)
  expect_equal(foo[[1]], expected_result[[1]])
})


# test prepare_high_confidence_objects()
test_that("prepare_high_confidence_objects returns objects of expected sizes for
          continuous data", {
  set.seed(1)
  num_samples <- 11
  temp_pheno <- as.matrix(c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1), ncol = 1)
  temp_tree <- ape::rtree(num_samples)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_trans <- expected_result <- rep(list(NULL), num_samples)
  for (i in 1:num_samples) {
    temp_trans[[i]]$transition <-
      c(0, 0, 0, 1, 0, 1, 1, 1,  1,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    temp_trans[[i]]$trans_dir <-
      c(0, 0, 0, 1, 0, 1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  }

  # FYI this geno does not match up with the fake transitions I made up
  temp_geno <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1),
                      ncol = num_samples,
                      nrow = num_samples)
  temp_tree$tip.label <- row.names(temp_geno) <- letters[1:num_samples]
  colnames(temp_geno) <- as.character(c(1:11))
  temp_discrete <- "discrete"
  temp_AR <- prepare_ancestral_reconstructions(temp_tree,
                                               temp_pheno,
                                               temp_geno,
                                               temp_discrete)

  for (i in 1:num_samples) {
    temp_AR$geno_recon_and_conf[[i]]$tip_and_node_rec_conf <-
      rep(1, ape::Ntip(temp_tree) + ape::Nnode(temp_tree))
  }

  # Keep only WT -> mutant transitions.
  temp_trans_phyc <- prep_geno_trans_for_phyc(temp_geno, temp_trans)

  temp_bootstrap <- 0.5
  temp_snps_per_gene <- NULL

  temp_geno_conf_order_by_edges <-
    temp_geno_recon_ord_by_edges <-
    rep(list(0), ncol(temp_geno))
  for (k in 1:ncol(temp_geno)) {
    temp_geno_conf_order_by_edges[[k]] <-
      reorder_tip_and_node_to_edge(
        temp_AR$geno_recon_and_conf[[k]]$tip_and_node_rec_conf,
        temp_tree)
    temp_geno_recon_ord_by_edges[[k]] <-
      reorder_tip_and_node_to_edge(
        temp_AR$geno_recon_and_conf[[k]]$tip_and_node_recon,
        temp_tree)
  }

  temp_AR$pheno_recon_and_conf$tip_and_node_rec_conf <-
    rep(1, ape::Nnode(temp_tree) + ape::Ntip(temp_tree))

  temp_hi_conf_orig <-
    prepare_high_confidence_objects(
      temp_trans_phyc,
      temp_tree,
      temp_AR$pheno_recon_and_conf$tip_and_node_rec_conf,
      temp_bootstrap,
      temp_geno,
      temp_geno_conf_order_by_edges,
      temp_geno_recon_ord_by_edges,
      temp_snps_per_gene)
  expect_equal(temp_hi_conf_orig$genotype_transition[[1]]$transition,
               temp_hi_conf_orig$genotype_transition[[1]]$trans_dir)
  expect_equal(temp_hi_conf_orig$genotype_transition[[1]]$transition,
               c(0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

})
