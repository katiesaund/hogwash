context("Tree manipulation functions") #---------------------------------------#

# TODO test reorder_tips_and_nodes_to_edges() -- see if this function ever gets used before writing unit tests
# reorder_tips_and_nodes_to_edges <- function(tips_and_node_vector, tr){
#   # Function description -------------------------------------------------------
#   # TODO ??
#   #
#   # Input:
#   # Tr. Phylo.
#   # tips_and_node_vector. ?
#   #
#   # Output:
#   # ordered_by_edges ?
#   #
#   # Check input ----------------------------------------------------------------
#   check_tree_is_valid(tr)
#   # TODO add check of length of edges vs tips_and_node_vector
#
#   # Function -------------------------------------------------------------------
#   ordered_by_edges <- rep(NA, Nedge(tr))
#   for (i in 1:Nedge(tr)){
#     ordered_by_edges[i] <- tips_and_node_vector[tr$edge[i, 2]]
#   }
#
#   # Return output --------------------------------------------------------------
#   return(ordered_by_edges)
# } # end reorder_tips_and_nodes_to_edges()

# test find_parent_edge
test_that("find_parent_edge returns edge 2 when given edge 3", {
  set.seed(1)
  temp_tree <- rtree(4)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_num <- 3
  expected_result <- 2
  expect_equal(find_parent_edge(temp_tree, temp_num), expected_result)
})

test_that("find_parent_edge returns NA when given edge 1", {
  set.seed(1)
  temp_tree <- rtree(4)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_num <- 1
  expected_result <- NA
  expect_equivalent(find_parent_edge(temp_tree, temp_num), expected_result)
})

test_that("find_parent_edge returns NA when given edge 2", {
  set.seed(1)
  temp_tree <- rtree(4)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_num <- 2
  expected_result <- NA
  expect_equal(find_parent_edge(temp_tree, temp_num), expected_result)
})

test_that("find_parent_edge returns an error when given edge 0", {
  set.seed(1)
  temp_tree <- rtree(4)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_num <- 0
  expect_error(find_parent_edge(temp_tree, temp_num))
})

test_that("find_parent_edge returns an error when an edge that doesn't exist", {
  set.seed(1)
  temp_tree <- rtree(4)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_num <- 100
  expect_error(find_parent_edge(temp_tree, temp_num))
})

# test identify_short_edges
test_that("identify_short_edges returns the only 1s in this test tree", {
  set.seed(1)
  temp_tree <- rtree(11)
  expect_identical(identify_short_edges(temp_tree), rep(1, Nedge(temp_tree)))
})

test_that("identify_short_edges returns a zero for the branch made artificially long in this test tree", {
  set.seed(1)
  temp_tree <- rtree(11)
  tree_with_long_edge <- temp_tree
  tree_with_long_edge$edge.length[9] <- .15 * sum(tree_with_long_edge$edge.length)
  short_edges <- identify_short_edges(tree_with_long_edge)
  expect_identical(short_edges[9], 0)
})

test_that("identify_short_edges returns an error on this tree, where iteratively each edge is too long.", {
  set.seed(1)
  temp_tree <- rtree(10)
  expect_error(identify_short_edges(temp_tree))
})

# test get_bootstrap_confidence()
test_that("get_bootstrap_confidence returns correct numbers for this dummy tree", {
  set.seed(1)
  temp_tree <- rtree(10)
  temp_tree$node.label <- c(91:99)
  temp_treshold <- .94
  expected_result <- c(rep(1, 10), rep(0, 3), rep(1, 6))
  expect_identical(get_bootstrap_confidence(temp_tree, temp_treshold), expected_result)
})

test_that("get_bootstrap_confidence returns correct numbers for this dummy tree", {
  set.seed(1)
  temp_tree <- rtree(10)
  temp_tree$node.label <- c(91:99)
  temp_treshold <- 1
  expected_result <- c(rep(1, 10), rep(0, 9))
  expect_identical(get_bootstrap_confidence(temp_tree, temp_treshold), expected_result)
})

test_that("get_bootstrap_confidence returns correct numbers for this dummy tree", {
  set.seed(1)
  temp_tree <- rtree(10)
  temp_tree$node.label <- c(91:99)
  temp_treshold <- .5
  expected_result <- c(rep(1, Nnode(temp_tree) + Ntip(temp_tree)))
  expect_identical(get_bootstrap_confidence(temp_tree, temp_treshold), expected_result)
})

test_that("get_bootstrap_confidence returns error for invalid threshold", {
  set.seed(1)
  temp_tree <- rtree(10)
  temp_tree$node.label <- c(91:99)
  temp_treshold <- 50
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})

test_that("get_bootstrap_confidence returns error for tree without node labels", {
  set.seed(1)
  temp_tree <- rtree(10)
  temp_treshold <- 0.5
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})

test_that("get_bootstrap_confidence returns error for tree without numeric node labels", {
  set.seed(1)
  temp_tree <- rtree(10)
  temp_tree$node.label <- letters[1:Nnode(temp_tree)]
  temp_treshold <- 0.5
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})

test_that("get_bootstrap_confidence returns error for incomplete node labels", {
  set.seed(1)
  temp_tree <- rtree(10)
  temp_tree$node.label[1:5] <- c(1:5)
  temp_treshold <- 0.5
  expect_error(get_bootstrap_confidence(temp_tree, temp_treshold))
})
# End of script ----------------------------------------------------------------