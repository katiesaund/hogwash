find_parent_edge <- function(tr, edge_num){
  # Function description -------------------------------------------------------
  # Given the number of a tree edge, return the parent edge number.
  #
  # Inputs:
  # tr. Phylo.
  # edge_num. Number.
  #
  # Outputs:
  # parent_edge. Number.
  #
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)
  check_is_number(edge_num)
  if (edge_num %% 1 != 0 | edge_num < 1){
    stop("Node number must be a positive integer")
  }
  if (edge_num > Nedge(tr)){
    stop("Number must be an edge in the tree.")
  }
  # Function -------------------------------------------------------------------

  # given an edge number, get edge number of parent edge
  parent_node <- tr$edge[edge_num, 1]
  parent_edge <- which(tr$edge[ , 2] == parent_node)

  # Return output --------------------------------------------------------------
  return(parent_edge)
  # TODO it breaks on 1st edge, need to deal with that, but not sure what's best yet
} # end find_parent_edge


identify_short_edges <- function(tr){
  # Function description -------------------------------------------------------
  # Removes any edges that make up a 10% or more of the total tr length.
  #
  # Input:
  # Tr. Phylo.
  #
  # Output:
  # short_edges. Vector of edge indices.
  #
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)

  # Function -------------------------------------------------------------------
  short_edges <- rep(1, Nedge(tr))
  while(max(tr$edge.length[as.logical(short_edges)]) >= (0.1 * sum(tr$edge.length[as.logical(short_edges)]))){
    short_edges[tr$edge.length == max(tr$edge.length[as.logical(short_edges)])] <- 0
  }

  # Return output --------------------------------------------------------------
  return(short_edges)
} # end identify_short_edges()

get_bootstrap_confidence <- function(tr, confidence_threshold){
  # Account for confidence in RAxML phylogenetic tree.
  node_confidence <- tr$node.label
  node_confidence <- as.numeric(node_confidence)/100
  node_confidence <- discretize_confidence_using_threshold(node_confidence, confidence_threshold)
  tree_tip_and_node_confidence <- c(rep(1, Ntip(tr)), node_confidence)
  if (length(tree_tip_and_node_confidence) != sum(1, Nedge(tr))){
    stop("tree confidence made incorrectly")
  }
  return(tree_tip_and_node_confidence)
} # end get_bootstrap_confidence()
