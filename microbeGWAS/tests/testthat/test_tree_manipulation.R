library(microbeGWAS)
context("Tree manipulation functions") #---------------------------------------#
# TODO add tests

# TODO find_parent_edge() -- see if this function ever gets used before writing unit tests


# find_parent_edge <- function(tr, edge_num){
#   # Function description -------------------------------------------------------
#   # Given the number of a tree edge, return the parent edge number.
#   #
#   # Inputs:
#   # tr. Phylo.
#   # edge_num. Number.
#   #
#   # Outputs:
#   # parent_edge. Number.
#   #
#   # Check input ----------------------------------------------------------------
#   check_tree_is_valid(tr)
#   check_for_root_and_bootstrap(tr)
#   check_is_number(edge_num)
#   if (edge_num %% 1 != 0 | edge_num < 1){
#     stop("Node number must be a positive integer")
#   }
#   if (edge_num > Nedge(tr)){
#     stop("Number must be an edge in the tree.")
#   }
#   # Function -------------------------------------------------------------------
#
#   # given an edge number, get edge number of parent edge
#   parent_node <- tr$edge[edge_num, 1]
#   parent_edge <- which(tr$edge[ , 2] == parent_node)
#
#   # Return output --------------------------------------------------------------
#   return(parent_edge)
#   # TODO it breaks on 1st edge, need to deal with that, but not sure what's best yet
# } # end find_parent_edge
#
#

# test identify_short_edges

test_that("identify_short_edges returns the longest branch in this test tree", {
  set.seed(1)
  temp_tree <- rtree(10)
  temp_tree$edge.length
  identify_short_edges(temp_tree)

  tree_with_long_edge <- temp_tree
  sum(tree_with_long_edge$edge.length)
  tree_with_long_edge$edge.length[1] <- 1
  identify_short_edges(tree_with_long_edge)
})


while(max(temp_tree$edge.length[as.logical(short_edges)]) >= (0.1 * sum(temp_tree$edge.length[as.logical(short_edges)]))){
  short_edges[temp_tree$edge.length == max(temp_tree$edge.length[as.logical(short_edges)])] <- FALSE
}
short_edges

# identify_short_edges <- function(tr){

#   # Function -------------------------------------------------------------------
#   short_edges <- rep(1, Nedge(tr))
#   while(max(tr$edge.length[as.logical(short_edges)]) >= (0.1 * sum(tr$edge.length[as.logical(short_edges)]))){
#     short_edges[tr$edge.length == max(tr$edge.length[as.logical(short_edges)])] <- 0
#   }
#
#   # Return output --------------------------------------------------------------
#   return(short_edges)
# } # end identify_short_edges()
#
# get_bootstrap_confidence <- function(tr, confidence_threshold){
#   # Account for confidence in RAxML phylogenetic tree.
#   node_confidence <- tr$node.label
#   node_confidence <- as.numeric(node_confidence)/100
#   node_confidence <- discretize_confidence_using_threshold(node_confidence, confidence_threshold)
#   tree_tip_and_node_confidence <- c(rep(1, Ntip(tr)), node_confidence)
#   if (length(tree_tip_and_node_confidence) != sum(1, Nedge(tr))){
#     stop("tree confidence made incorrectly")
#   }
#   return(tree_tip_and_node_confidence)
# } # end get_bootstrap_confidence()
#
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
