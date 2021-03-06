# Ancestral reconstruction ----------------------------------------------------#

# ancestral_reconstruction_by_ML -----------------------------------------------
test_that("ancestral_reconstruction_by_ML with discrete input works", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(9, rooted = TRUE)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  num_col <- 9
  num_cells <- num_col * ape::Ntip(temp_tree)
  test_mat <-
    matrix(rep(c(1, 0), num_cells), nrow = ape::Ntip(temp_tree), ncol = num_col)
  dummy_pheno <-
    ancestral_reconstruction_by_ML(temp_tree, test_mat, 1, "discrete")
  dummy_geno <-
    ancestral_reconstruction_by_ML(temp_tree, test_mat, 2, "discrete")

  # Test
  expected_length <- ape::Ntip(temp_tree) + ape::Nnode(temp_tree)
  expect_identical(length(dummy_pheno$tip_and_node_recon), expected_length)
  expect_identical(length(dummy_geno$tip_and_node_recon), expected_length)
  expect_identical(length(dummy_pheno$tip_and_node_rec_conf), expected_length)
  expect_identical(length(dummy_geno$tip_and_node_rec_conf), expected_length)

  expected_length <- ape::Nnode(temp_tree)
  expect_identical(length(dummy_pheno$node_anc_rec), expected_length)
  expect_identical(length(dummy_geno$node_anc_rec), expected_length)

  expected_rows <- nrow(temp_tree$edge)
  expected_columns <- ncol(temp_tree$edge)
  expect_identical(nrow(dummy_pheno$recon_edge_mat), expected_rows)
  expect_identical(ncol(dummy_pheno$recon_edge_mat), expected_columns)
  expect_identical(nrow(dummy_geno$recon_edge_mat), expected_rows)
  expect_identical(ncol(dummy_geno$recon_edge_mat), expected_columns)
})

test_that("ancestral_reconstruction_by_ML with continuous input works.", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(9, rooted = TRUE)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  num_col <- 9
  num_cells <- num_col * ape::Ntip(temp_tree)
  test_mat <- matrix(rnorm(num_cells, mean = 0, sd = 10),
                     nrow = ape::Ntip(temp_tree),
                     ncol = num_col)
  test_pheno_mat_1 <- test_mat[, 1, drop = FALSE]
  test_pheno_mat_2 <- test_mat[, 2, drop = FALSE]
  dummy_pheno1 <- ancestral_reconstruction_by_ML(temp_tree,
                                                 test_pheno_mat_1,
                                                 1,
                                                 "continuous")
  dummy_pheno2 <-  ancestral_reconstruction_by_ML(temp_tree,
                                                  test_pheno_mat_2,
                                                  1,
                                                  "continuous")

  # Test
  expected_length <- ape::Ntip(temp_tree) + ape::Nnode(temp_tree)
  expect_identical(length(dummy_pheno1$tip_and_node_recon), expected_length)
  expect_identical(length(dummy_pheno2$tip_and_node_recon), expected_length)
  expect_identical(length(dummy_pheno1$tip_and_node_rec_conf), expected_length)
  expect_identical(length(dummy_pheno2$tip_and_node_rec_conf), expected_length)

  expected_length <- ape::Nnode(temp_tree)
  expect_identical(length(dummy_pheno1$node_anc_rec), expected_length)
  expect_identical(length(dummy_pheno2$node_anc_rec), expected_length)

  expected_rows <- nrow(temp_tree$edge)
  expected_columns <- ncol(temp_tree$edge)
  expect_silent(check_class(dummy_pheno1$recon_edge_mat, "matrix"))
  expect_identical(nrow(dummy_pheno1$recon_edge_mat), expected_rows)
  expect_identical(ncol(dummy_pheno1$recon_edge_mat), expected_columns)
})

test_that("ancestral_reconstruction_by_ML with discrete input works.", {
  # Set up
  num_col <- 2
  num_tips <- 5
  set.seed(1)
  tree <- ape::rtree(num_tips, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- matrix(c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1),
                     nrow = ape::Ntip(tree),
                     ncol = num_col)
  dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "discrete")
  dummy_geno <- ancestral_reconstruction_by_ML(tree, test_mat, 2, "discrete")

  # Test
  expect_equivalent(dummy_pheno$tip_and_node_recon,
                    c(1, 1, 1, 1, 0, 1, 1, 1, 1))
  expect_equivalent(dummy_pheno$node_anc_rec,
                    c(1, 1, 1, 1))
  expect_equivalent(dummy_pheno$tip_and_node_rec_conf,
                    c(1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equivalent(dummy_pheno$recon_edge_mat[, 1, drop = TRUE],
                    c(1, 1, 1, 1, 1, 1, 1, 1))
  expect_equivalent(dummy_pheno$recon_edge_mat[, 2, drop = TRUE],
                    c(1, 1, 1, 1, 1, 1, 1, 0))
  expect_equivalent(dummy_geno$tip_and_node_recon,
                    c(0, 0, 0, 1, 1, 0, 0, 0, 1))
  expect_equivalent(dummy_geno$node_anc_rec,
                    c(0, 0, 0, 1))
  expect_equivalent(dummy_geno$tip_and_node_rec_conf,
                    c(1, 1, 1, 1, 1, 0, 0, 0, 0))
  expect_equivalent(dummy_geno$recon_edge_mat[, 1, drop = TRUE],
                    c(0, 0, 0, 0, 0, 1, 1, 0))
  expect_equivalent(dummy_geno$recon_edge_mat[, 2, drop = TRUE],
                    c(0, 0, 0, 0, 1, 0, 1, 1))
})


test_that("ancestral_reconstruction_by_ML produce errors given bogus input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(9, rooted = TRUE)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))

  # Test
  expect_error(ancestral_reconstruction_by_ML(temp_tree,
                                              "foobar",
                                              1,
                                              "discrete"))
})

# test continuous_ancestral_reconstruction -------------------------------------
test_that("continuous_ancestral_reconstruction works given valid inputs", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(10, rooted = TRUE)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  num_col <- 10
  num_cells <- num_col * ape::Ntip(temp_tree)
  test_mat <- matrix(rnorm(num_cells, mean = 0, sd = 10),
                     nrow = ape::Ntip(temp_tree),
                     ncol = num_col)
  test_pheno_mat_1 <- test_mat[, 1, drop = FALSE]
  test_pheno_mat_2 <- test_mat[, 2, drop = FALSE]
  reconstruction_method <- "REML"
  data_type <- "continuous"
  index <- 1
  dummy_pheno1 <- continuous_ancestral_reconstruction(temp_tree,
                                                      test_pheno_mat_1,
                                                      index,
                                                      data_type,
                                                      reconstruction_method)
  dummy_pheno2 <-  continuous_ancestral_reconstruction(temp_tree,
                                                       test_pheno_mat_2,
                                                       index,
                                                       data_type,
                                                       reconstruction_method)

  # Test
  expected_length <- ape::Ntip(temp_tree) + ape::Nnode(temp_tree)
  expect_identical(length(dummy_pheno1$tip_and_node_recon), expected_length)
  expect_identical(length(dummy_pheno2$tip_and_node_recon), expected_length)

  expected_length <- ape::Nnode(temp_tree)
  expect_identical(length(dummy_pheno1$ML_anc_rec), expected_length)
  expect_identical(length(dummy_pheno2$ML_anc_rec), expected_length)

  # tips should be included in tip and node reconstruction
  expect_identical(
    unname(dummy_pheno1$tip_and_node_recon[1:ape::Ntip(temp_tree)]),
    test_pheno_mat_1[, 1, drop = TRUE])

  # node reconstruction should be included in tip and node reconstruction
  expect_identical(
    dummy_pheno1$tip_and_node_recon[(
      ape::Ntip(temp_tree) + 1):(ape::Ntip(temp_tree) + ape::Nnode(temp_tree))],
    dummy_pheno1$ML_anc_rec)
})

# test continuous_get_recon_confidence -----------------------------------------
test_that("continuous_get_recon_confidence works for valid input", {
  # Set up
  temp_reconstruction_vector <- c(1:10)
  temp_conf <- continuous_get_recon_confidence(temp_reconstruction_vector)

  # Test
  expect_identical(temp_conf, rep(1, length(temp_reconstruction_vector)))
})

test_that("continuous_get_recon_confidence throws error with invalid input", {
  # Test
  expect_error(continuous_get_recon_confidence(matrix(1, 10, 10)))
  expect_error(continuous_get_recon_confidence("foobar"))
  expect_error(continuous_get_recon_confidence(NA))
})

# test convert_to_edge_mat -----------------------------------------------------
test_that("convert_to_edge_mat works for valid input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
  test_vector <- rep(1, ape::Ntip(temp_tree) + ape::Nnode(temp_tree))
  output_matrix <- matrix(1, ncol = 2, nrow = ape::Nedge(temp_tree))
  results <- convert_to_edge_mat(temp_tree, test_vector)

  # Test
  expect_error(results, NA)
  expect_equivalent(results, output_matrix)
  expect_equal(nrow(results), ape::Nedge(temp_tree))
  expect_equal(ncol(results), 2)

  # Set up
  temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
  test_vector <- c(1:sum(ape::Ntip(temp_tree) + ape::Nnode(temp_tree)))
  test_vector <- test_vector + 10
  output_matrix <- temp_tree$edge + 10

  # Test
  expect_error(convert_to_edge_mat(temp_tree, test_vector), NA)
  expect_equivalent(convert_to_edge_mat(temp_tree, test_vector), output_matrix)

  # Set up
  test_vector <- rep(1:9)
  expected_parent_nodes <- c(6, 6, 7, 8, 8, 9, 9, 7)
  expected_child_nodes <- c(1, 7, 8, 2, 9, 3, 4, 5)
  expected_output <- cbind(expected_parent_nodes, expected_child_nodes)
  results <- convert_to_edge_mat(temp_tree, test_vector)
  expect_equal(unname(expected_output), results)
})

test_that("convert_to_edge_mat gives an error when not given a tree", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
  test_vector <- rep(1, ape::Ntip(temp_tree) + ape::Nnode(temp_tree))
  output_matrix <- matrix(1, ncol = 2, nrow = ape::Nedge(temp_tree))

  # Test
  not_a_tree <- "foobar"
  expect_error(convert_to_edge_mat(not_a_tree, test_vector))
  not_a_tree <- matrix(1, 10, 10)
  expect_error(convert_to_edge_mat(not_a_tree, test_vector))
  not_a_tree <- NA
  expect_error(convert_to_edge_mat(not_a_tree, test_vector))
  not_a_tree <- NULL
  expect_error(convert_to_edge_mat(not_a_tree, test_vector))
})

test_that("convert_to_edge_mat gives an error given invalid input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))

  # Test
  not_the_reconstruction <- rep(1, ape::Ntip(temp_tree))
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
  not_the_reconstruction <- NA
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
  not_the_reconstruction <- NULL
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
  not_the_reconstruction <- matrix(1, 10, 10)
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
  not_the_reconstruction <- "foobar"
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
})

# test discrete_ancestral_reconstruction ---------------------------------------
test_that("discrete_ancestral_reconstruction gives results of expected size", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(9, rooted = TRUE)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  num_col <- 9
  num_cells <- num_col * ape::Ntip(temp_tree)
  test_mat <- matrix(rep(c(1, 0), num_cells),
                     nrow = ape::Ntip(temp_tree),
                     ncol = num_col)
  reconstruction_method <- "ML"
  dummy_pheno <- discrete_ancestral_reconstruction(temp_tree,
                                                   test_mat,
                                                   1,
                                                   "discrete",
                                                   reconstruction_method)
  dummy_geno <- discrete_ancestral_reconstruction(temp_tree,
                                                  test_mat,
                                                  2,
                                                  "discrete",
                                                  reconstruction_method)

  # Test
  expected_length <- ape::Ntip(temp_tree) + ape::Nnode(temp_tree)
  expect_identical(length(dummy_pheno$tip_and_node_recon), expected_length)
  expect_identical(length(dummy_geno$tip_and_node_recon), expected_length)

  expected_length <- ape::Nnode(temp_tree)
  expect_identical(length(dummy_pheno$ML_anc_rec), expected_length)
  expect_identical(length(dummy_geno$ML_anc_rec), expected_length)

  expect_type(dummy_pheno$reconstruction, "list")
  expect_silent(check_class(dummy_pheno$reconstruction, "ace"))

  expect_identical(nrow(dummy_pheno$reconstruction$lik.anc),
                   ape::Nnode(temp_tree))
  expect_equal(ncol(dummy_pheno$reconstruction$lik.anc), 2)
})

test_that("discrete_ancestral_reconstruction gives error for invalid input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(9, rooted = TRUE)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  num_col <- 9
  num_cells <- num_col * ape::Ntip(temp_tree)
  test_mat <- matrix(rep(c(1, 0), num_cells),
                     nrow = ape::Ntip(temp_tree),
                     ncol = num_col)
  reconstruction_method <- "ML"
  expect_error(discrete_ancestral_reconstruction(temp_tree,
                                                   test_mat,
                                                   1,
                                                   "foo",
                                                   reconstruction_method))
})

test_that("discrete_ancestral_reconstruction gives error for invalid input", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(9, rooted = TRUE)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  num_col <- 9
  num_cells <- num_col * ape::Ntip(temp_tree)
  test_mat <- matrix(rep(c(1, 0), num_cells),
                     nrow = ape::Ntip(temp_tree),
                     ncol = num_col)
  reconstruction_method <- "REML"
  expect_error(discrete_ancestral_reconstruction(temp_tree,
                                                 test_mat,
                                                 1,
                                                 "continuous",
                                                 reconstruction_method))
})

test_that("discrete_ancestral_reconstruction throws error given bogus inputs", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(9, rooted = TRUE)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  num_col <- 9
  num_cells <- num_col * ape::Ntip(temp_tree)
  test_mat <- matrix(rep(c(1, 0), num_cells),
                     nrow = ape::Ntip(temp_tree),
                     ncol = num_col)
  reconstruction_method <- "ML"

  # Test
  expect_error(discrete_ancestral_reconstruction(temp_tree,
                                                 test_mat,
                                                 1,
                                                 "foobar",
                                                 reconstruction_method))
  expect_error(discrete_ancestral_reconstruction(temp_tree,
                                                 "foobar",
                                                 2,
                                                 "discrete",
                                                 reconstruction_method))

})

# test build_better_reconstruction ---------------------------------------------
test_that("build_better_reconstruction gives an error given invalid input", {
  # Set up
  num_col <- 9
  tree <- ape::rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- matrix(1:num_cells, nrow = ape::Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"

  # Test
  expect_error(build_better_reconstruction(test_mat,
                                           tree,
                                           discrete_or_continuous,
                                           index,
                                           reconstruction_method))
})

test_that("build_better_reconstruction works given valid input", {
  # Set up
  num_col <- 9
  tree <- ape::rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- matrix(rep(c(1, 0), num_cells),
                     nrow = ape::Ntip(tree),
                     ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"

  # Test
  expect_error(build_better_reconstruction(test_mat,
                                           tree,
                                           discrete_or_continuous,
                                           index,
                                           reconstruction_method),
               NA)
})

test_that("build_better_reconstruction gives an error given invalid input", {
  # Set up
  num_col <- 9
  tree <- ape::rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- matrix(1:num_cells, nrow = ape::Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "continuous"
  index <- 1
  reconstruction_method <- "REML"

  # Test
  expect_error(build_better_reconstruction(test_mat,
                                           tree,
                                           discrete_or_continuous,
                                           index,
                                           reconstruction_method))
})

test_that("build_better_reconstruction gives an error given invalid input", {
  # Set up
  num_col <- 9
  tree <- ape::rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- matrix(0, nrow = ape::Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "continuous"
  index <- 1
  reconstruction_method <- "REML"

  # Test
  expect_error(build_better_reconstruction(test_mat,
                                           tree,
                                           discrete_or_continuous,
                                           index,
                                           reconstruction_method))
})


test_that("build_better_reconstruction gives an error for invalid input", {
  # Set up
  num_col <- 9
  tree <- ape::rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- matrix(rep(c(1, 0), num_cells),
                     nrow = ape::Ntip(tree),
                     ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "fake_method"

  # Test
  expect_error(build_better_reconstruction(test_mat,
                                           tree,
                                           discrete_or_continuous,
                                           index,
                                           reconstruction_method))
})

test_that("build_better_reconstruction gives an error when it given a bad index", {
  # Set up
  num_col <- 9
  tree <- ape::rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- matrix(rep(c(1, 0), num_cells),
                     nrow = ape::Ntip(tree),
                     ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- num_col + 10
  reconstruction_method <- "ML"

  # Test
  expect_error(build_better_reconstruction(test_mat,
                                           tree,
                                           discrete_or_continuous,
                                           index,
                                           reconstruction_method))
})

test_that("build_better_reconstruction gives an error with invalid input", {
  # Set up
  num_col <- 9
  tree <- ape::rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- rep(c(1, 0), num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"

  # Test
  expect_error(build_better_reconstruction(test_mat,
                                           tree,
                                           discrete_or_continuous,
                                           index,
                                           reconstruction_method))
})

test_that("build_better_reconstruction returns 'ER' for this test data", {
  # Set up
  num_col <- 9
  set.seed(1)
  tree <- ape::rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- matrix(rep(c(1, 0), num_cells),
                     nrow = ape::Ntip(tree),
                     ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"
  reconstruction <- build_better_reconstruction(test_mat,
                                                tree,
                                                discrete_or_continuous,
                                                index,
                                                reconstruction_method)
  # Test
  expect_identical(reconstruction$call$model, "ER")
})

test_that("build_better_reconstruction returns 'ARD' for this test particular data", {
  # Set up
  num_col <- 100
  set.seed(1)
  tree <- ape::rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, ape::Nnode(tree))
  num_cells <- num_col * ape::Ntip(tree)
  test_mat <- matrix(rep(c(1, 1, 1, 1, 1, 1, 1, 0), num_cells),
                     nrow = ape::Ntip(tree),
                     ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"
  reconstruction <- build_better_reconstruction(test_mat,
                                                tree,
                                                discrete_or_continuous,
                                                index,
                                                reconstruction_method)
  # Test
  expect_identical(reconstruction$call$model, "ARD")
})

# test prepare_ancestral_reconstructions ---------------------------------------
test_that("prepare_ancestral_reconstructions works for valid input", {
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
  temp_geno <-
    cbind(genotype1, genotype2, genotype5, genotype6, genotype7, genotype8)
  row.names(temp_geno) <- temp_tree$tip.label
  colnames(temp_geno) <- c("SNP1", "SNP2", "SNP5", "SNP6", "SNP7", "SNP8")
  temp_continuous <- "continuous"
  temp_AR <- prepare_ancestral_reconstructions(temp_tree,
                                               temp_pheno,
                                               temp_geno,
                                               temp_continuous)

  # Test
  expect_equal(length(temp_AR$geno_trans), ncol(temp_geno))
  expect_equal(length(temp_AR$geno_trans[[1]]$transition),
               ape::Nedge(temp_tree))
  expect_equal(length(temp_AR$geno_trans[[1]]$trans_dir),
               ape::Nedge(temp_tree))
  expect_equal(length(temp_AR$geno_recon_and_conf[[1]]$tip_and_node_rec_conf),
               (ape::Ntip(temp_tree) + ape::Nnode(temp_tree)))
  expect_equal(length(temp_AR$geno_recon_and_conf[[1]]$tip_and_node_recon),
               (ape::Ntip(temp_tree) + ape::Nnode(temp_tree)))
  expect_equal(length(temp_AR$geno_recon_and_conf[[1]]$node_anc_rec),
               ape::Nnode(temp_tree))
  expect_equal(nrow(temp_AR$geno_recon_and_conf[[1]]$recon_edge_mat),
               ape::Nedge(temp_tree))
  expect_equal(ncol(temp_AR$geno_recon_and_conf[[1]]$recon_edge_mat), 2)
  expect_equal(length(temp_AR$pheno_recon_and_conf$tip_and_node_rec_conf),
               (ape::Ntip(temp_tree) + ape::Nnode(temp_tree)))
  expect_equal(length(temp_AR$pheno_recon_and_conf$tip_and_node_recon),
               (ape::Ntip(temp_tree) + ape::Nnode(temp_tree)))
  expect_equal(length(temp_AR$pheno_recon_and_conf$node_anc_rec),
               ape::Nnode(temp_tree))
  expect_equal(nrow(temp_AR$pheno_recon_and_conf$recon_edge_mat),
               ape::Nedge(temp_tree))
  expect_equal(ncol(temp_AR$pheno_recon_and_conf$recon_edge_mat), 2)
  expect_equal(temp_AR$pheno_recon_and_conf$tip_and_node_rec_conf,
               rep(1, ape::Ntip(temp_tree) + ape::Nnode(temp_tree)))
  expect_equal(
    unname(
      temp_AR$pheno_recon_and_conf$tip_and_node_recon[1:ape::Ntip(temp_tree)]),
    unname(temp_pheno[, 1, drop = TRUE]))
})
