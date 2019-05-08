find_parent_edge <- function(tr, edge_num){
  # Function description -------------------------------------------------------
  # Given the number of a tree edge, return the parent edge number.
  #
  # Inputs:
  # tr. Phylo.
  # edge_num. Number.
  #
  # Outputs:
  # parent_edge. Number or NA. NA when given edge doesn't have a parent.
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
  edge_with_basal_node <- Ntip(tr) + 1
  if (tr$edge[edge_num, 1] == edge_with_basal_node){
    parent_edge <- NA
  } else {
    parent_node <- tr$edge[edge_num, 1]
    parent_edge <- which(tr$edge[ , 2] == parent_node)
  }

  # Return output --------------------------------------------------------------
  return(parent_edge)
} # end find_parent_edge


identify_short_edges <- function(tr){
  # Function description -------------------------------------------------------
  # Removes any edges that make up a 10% or more of the total tr length.
  #
  # Input:
  # Tr. Phylo.
  #
  # Output:
  # short_edges. Binary vector. Length = Nedge(tr). Each entry corresponds to
  #                             an edge. 0 == long edge (low confidence in
  #                             reconstruction value). 1 == short edge (high
  #                             confidence in reconstruction value).
  #
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)

  # Function -------------------------------------------------------------------
  short_edges <- rep(1, Nedge(tr))
  while(max(tr$edge.length[as.logical(short_edges)]) >= (0.1 * sum(tr$edge.length[as.logical(short_edges)]))){
    short_edges[tr$edge.length == max(tr$edge.length[as.logical(short_edges)])] <- 0
    if (sum(short_edges) == 0){
      stop("Tree edge lengths are unreasonably long compared to the other edges.")
    }
  }

  # Return output --------------------------------------------------------------
  check_if_binary_vector(short_edges)
  if (length(short_edges) != Nedge(tr)){stop("Short_edges should have length == Nedge(tr)")}
  return(short_edges)
} # end identify_short_edges()

get_bootstrap_confidence <- function(tr, confidence_threshold){
  # Function description -------------------------------------------------------
  # Extracts the bootstrap confidence stored in the node.labels of the tree,
  # convertns them to a fraction between 0 and 1, and then classifies each
  # node as high or low confidence based on the input confidence threshold.
  #
  # Notes:
  # This method requires that the bootstrap values be assigned to node
  # labels.
  #
  # It makes an assumption that the bootstrap values are stored as numbers
  # between 1 and 100, that need to then be converted to from this percentage
  # to a fraction.
  #
  # A better input would be Bayesian posterior probabilities than ML derived
  # bootstrap values...but most of what we're working with is boostrap values.
  #
  # Input:
  # Tr. Phylo.
  # confidence_threshold. Numeric. Between 0 and 1.
  #
  # Output:
  # tree_tip_and_node_confidence.  Binary vector. Length = Ntip(tr) + Nnode(tr).
  #                                The first section is the tip confidence,
  #                                which is always, by definition 1. Then the
  #                                remaining entries correspond to the nodes
  #                                and is either 0 or 1.
  #
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_is_number(confidence_threshold)
  check_if_alpha_valid(confidence_threshold)
  if (max(tr$node.label) > 100 | min(tr$node.label) < 0){
    stop("Tree$node.label are expected to be positive numbers between 0 and 100")
  }

  # Function -------------------------------------------------------------------
  node_confidence <- tr$node.label

  if (max(node_confidence) > 1){
    node_confidence <- as.numeric(node_confidence)/100
  }
  node_confidence <- discretize_confidence_using_threshold(node_confidence, confidence_threshold)
  tree_tip_and_node_confidence <- c(rep(1, Ntip(tr)), node_confidence)

  # Check and return output ----------------------------------------------------
  check_if_binary_vector(tree_tip_and_node_confidence)
  if (length(tree_tip_and_node_confidence) != sum(Ntip(tr) + Nnode(tr))){
    stop("tree confidence made incorrectly")
  }

  return(tree_tip_and_node_confidence)
} # end get_bootstrap_confidence()

reorder_tips_and_nodes_to_edges <- function(tips_and_node_vector, tr){
  # Function description -------------------------------------------------------
  # Convert
  #
  # Input:
  # Tr. Phylo.
  # tips_and_node_vector. ?
  #
  # Output:
  # ordered_by_edges ?
  #
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  # TODO add check of length of edges vs tips_and_node_vector

  # Function -------------------------------------------------------------------
  ordered_by_edges <- rep(NA, Nedge(tr))
  for (i in 1:Nedge(tr)){
    ordered_by_edges[i] <- tips_and_node_vector[tr$edge[i, 2]]
  }

  # Return output --------------------------------------------------------------
  return(ordered_by_edges)
} # end reorder_tips_and_nodes_to_edges()

