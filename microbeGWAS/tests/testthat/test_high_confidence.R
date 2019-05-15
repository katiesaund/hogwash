library(microbeGWAS)
context("High confidence") #----------------------------------------------#
# TODO write unit tests for assign_high_confidence_to_transition_edges

# test discretize_confidence_using_threshold
test_that("discretize_confidence_using_threshold should give this expected result", {
  temp_vector <- c(1:100)
  temp_threshold <- 50
  expect_equal(discretize_confidence_using_threshold(temp_vector, temp_threshold), c(rep(0, 49), rep(1, 51)))
})

test_that("discretize_confidence_using_threshold should give this expected result - now with fractions", {
  temp_vector <- seq(from = 0, to = 1, by = 0.1)
  temp_threshold <- 0.5
  expect_equal(discretize_confidence_using_threshold(temp_vector, temp_threshold), c(rep(0, 5), rep(1, 6)))
})

# test report_num_high_confidence_trans_edge
test_that("report_num_high_confidence_trans_edge returns expected outcome for this test set", {
  fake_geno_names <- letters[1:3]
  fake_trans <- fake_hi_conf_edges <- rep(list(NULL), length(fake_geno_names))
  for (i in 1:length(fake_geno_names)){
    fake_trans[[i]]$transition <- c(1,1,1)
    fake_hi_conf_edges[[i]] <- c(1, 0, 1)
  }
  fake_trans[[1]]$transition <- c(1, 1, 1)
  fake_trans[[2]]$transition <- c(0, 0, 0)
  fake_trans[[3]]$transition <- c(1, 1, 0)
  expected_result <- c(2, 0, 1)
  names(expected_result) <- fake_geno_names
  expect_equal(report_num_high_confidence_trans_edge(fake_trans, fake_hi_conf_edges,fake_geno_names), expected_result)
})

# test assign_high_confidence_to_transition_edges_including_parent_info
test_that("assign_high_confidence_to_transition_edges_including_parent_info returns the edges that are high confidence transition edges for this tree", {
  set.seed(1)
  num_samples <- 5
  temp_tree <- rtree(num_samples)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  plot(temp_tree)
  edgelabels()
  temp_trans <- temp_confidence <- expected_result <- rep(list(NULL), num_samples)
  for (i in 1:num_samples){
    temp_trans[[i]]$transition <- c(0, 0, 0, 1, 0, 1, 1, 0)
    temp_confidence[[i]] <- c(0, 0, 0, 0, 0, 1, 1, 0)
    expected_result[[i]] <- c(0, 0, 0, 0, 0, 0, 1, 0)
  }
  temp_geno <- matrix(1, ncol = num_samples, nrow = num_samples) # FYI this geno does not match up with the fake transitions I made up
  temp_tree$tip.label <- row.names(temp_geno) <- letters[1:num_samples]
  foo <- assign_high_confidence_to_transition_edges(temp_tree, temp_confidence, temp_trans, temp_geno)
  expect_equal(foo[[1]], expected_result[[1]])
})
