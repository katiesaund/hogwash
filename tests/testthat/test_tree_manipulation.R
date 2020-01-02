context("Tree manipulation functions") #---------------------------------------#

# test reorder_tip_and_node_to_edge
test_that("reorder_tip_and_node_to_edge works for valid inputs", {
  # Set up
  # Assume a tree with 5 tips, 8 edges, and 4 nodes
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_tip_node <- c(rep(1, ape::Ntip(temp_tree)),
                     rep(0, ape::Nnode(temp_tree)))

  results <- reorder_tip_and_node_to_edge(temp_tip_node, temp_tree)

  # Test
  expect_equal(results, c(1, 0, 0, 1, 0, 1, 1, 1))
})

test_that("reorder_tip_and_node_to_edge gives errors for invalid inputs", {
  # Set up
  # Assume a tree with 5 tips, 8 edges, and 4 nodes
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_tip_node <- c(rep(1, ape::Ntip(temp_tree)),
                     rep(0, ape::Nnode(temp_tree)))
  temp_tree <- "bad_tree"

  # Test
  expect_error(reorder_tip_and_node_to_edge(temp_tip_node, temp_tree))
})

test_that("reorder_tip_and_node_to_edge gives errors for invalid inputs", {
  # Set up
  # Assume a tree with 5 tips, 8 edges, and 4 nodes
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_tip_node <- letters[1:9]

  # Test
  expect_error(reorder_tip_and_node_to_edge(temp_tip_node, temp_tree))
})

# test identify_short_edges
test_that("identify_short_edges returns the only 1s in this test tree", {
  set.seed(1)
  temp_tree <- ape::rtree(11)
  expect_identical(identify_short_edges(temp_tree),
                   rep(1, ape::Nedge(temp_tree)))
})

test_that("identify_short_edges returns a zero for the branch made artificially
          long in this test tree", {
  set.seed(1)
  temp_tree <- ape::rtree(15)
  tree_with_long_edge <- temp_tree
  tree_with_long_edge$edge.length[9] <-
    .15 * sum(tree_with_long_edge$edge.length)
  short_edges <- identify_short_edges(tree_with_long_edge)
  expect_identical(short_edges[9], 0)
})

test_that("identify_short_edges returns an error on this tree, where iteratively
          each edge is too long.", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  expect_error(identify_short_edges(temp_tree))
})

# test get_bootstrap_confidence()
test_that("get_bootstrap_confidence returns correct numbers for this dummy
          tree", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  temp_tree$node.label <- c(91:99)
  temp_treshold <- .94
  expected_result <- c(rep(1, 10), rep(0, 3), rep(1, 6))
  expect_identical(get_bootstrap_confidence(temp_tree, temp_treshold),
                   expected_result)
})

test_that("get_bootstrap_confidence returns correct numbers for this dummy
          tree", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  temp_tree$node.label <- c(91:99)
  temp_treshold <- 1
  expected_result <- c(rep(1, 10), rep(0, 9))
  expect_identical(get_bootstrap_confidence(temp_tree, temp_treshold),
                   expected_result)
})

test_that("get_bootstrap_confidence returns correct numbers for this dummy
          tree", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  temp_tree$node.label <- c(91:99)
  temp_treshold <- .5
  expected_result <- c(rep(1, ape::Nnode(temp_tree) + ape::Ntip(temp_tree)))
  expect_identical(get_bootstrap_confidence(temp_tree, temp_treshold),
                   expected_result)
})

test_that("get_bootstrap_confidence returns error for invalid threshold", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  temp_tree$node.label <- c(91:99)
  temp_treshold <- 50
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})

test_that("get_bootstrap_confidence returns error for tree without node
          labels", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  temp_treshold <- 0.5
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})

test_that("get_bootstrap_confidence returns error for tree without numeric node
          labels", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  temp_tree$node.label <- letters[1:ape::Nnode(temp_tree)]
  temp_treshold <- 0.5
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})

test_that("get_bootstrap_confidence returns error for incomplete node labels", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  temp_tree$node.label[1:5] <- c(1:5)
  temp_treshold <- 0.5
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})


test_that("get_bootstrap_confidence returns error for invalid input", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  temp_tree$node.label[1:10] <- c(92:101)
  temp_treshold <- 0.5
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})

test_that("get_bootstrap_confidence returns error for invalid input", {
  set.seed(1)
  temp_tree <- ape::rtree(10)
  temp_tree$node.label[1:10] <- c(-5:4)
  temp_treshold <- 0.5
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})
