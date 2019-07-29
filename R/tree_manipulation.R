#' identify_short_edges
#'
#' @description Removes any edges that make up a 10% or more of the total tr
#'  length.
#' @param tr  Phylo.
#'
#' @return short_edges. Binary vector. Length = Nedge(tr). Each entry
#'  corresponds to an edge. 0 == long edge (low confidence in reconstruction
#'  value). 1 == short edge (high confidence in reconstruction value).
#' @noRd
identify_short_edges <- function(tr){
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)

  # Function -------------------------------------------------------------------
  short_edges <- rep(1, ape::Nedge(tr))
  while (max(tr$edge.length[as.logical(short_edges)]) >=
         (0.1 * sum(tr$edge.length[as.logical(short_edges)]))) {
    short_edges[tr$edge.length ==
                  max(tr$edge.length[as.logical(short_edges)])] <- 0
    if (sum(short_edges) == 0) {
      stop("Tree edge lengths are unreasonably long compared to the other
           edges.")
    }
  }

  # Return output --------------------------------------------------------------
  check_if_binary_vector(short_edges)
  check_equal(length(short_edges), ape::Nedge(tr))
  return(short_edges)
} # end identify_short_edges()

#' get_bootstrap_confidence
#'
#' @description Extracts the bootstrap confidence stored in the node.labels of
#'  the tree, convertns them to a fraction between 0 and 1, and then classifies
#'  each node as high or low confidence based on the input confidence threshold.
#'
#' @details This method requires that the bootstrap values be assigned to node
#'  labels. It makes an assumption that the bootstrap values are stored as
#'   numbers between 1 and 100, that need to then be converted to from this
#'   percentage to a fraction. A better input would be Bayesian posterior
#'   probabilities than ML derived bootstrap values...but most of what we're
#'   working with is boostrap values.
#'
#' @param tr Phylo.
#' @param confidence_threshold Numeric. Between 0 and 1.
#'
#' @return tree_tip_and_node_confidence. Binary vector. Length = Ntip(tr) +
#'  Nnode(tr). The first section is the tip confidence, which is always, by
#'  definition 1. Then the remaining entries correspond to the nodes and is
#'  either 0 or 1.
#'
#'  @noRd
#'
get_bootstrap_confidence <- function(tr, confidence_threshold){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_is_number(confidence_threshold)
  check_num_between_0_and_1(confidence_threshold)
  if (max(tr$node.label) > 100 | min(tr$node.label) < 0) {
    stop("Tree$node.label are expected to be positive numbers between 0-100")
  }

  # Function -------------------------------------------------------------------
  node_confidence <- tr$node.label

  if (max(node_confidence) > 1) {
    node_confidence <- as.numeric(node_confidence) / 100
  }
  node_confidence <-
    discretize_conf_with_cutoff(node_confidence, confidence_threshold)
  tree_tip_and_node_confidence <- c(rep(1, ape::Ntip(tr)), node_confidence)

  # Check and return output ----------------------------------------------------
  check_if_binary_vector(tree_tip_and_node_confidence)
  check_equal(length(tree_tip_and_node_confidence),
              sum(ape::Ntip(tr) + ape::Nnode(tr)))
  return(tree_tip_and_node_confidence)
} # end get_bootstrap_confidence()

#' reorder_tip_and_node_to_edge
#'
#' @description Reorder a vector organzied by tips then nodes into a vector
#'  organized by tree edges.
#'
#'  @details This function grabs the value of each child from every node and
#'   stores them in the order of the tree's edge matrix. This effectively
#'   drops the value of the tree's root because the root is never a child and
#'   therefore not in the child column of the edge matrix. Each node/tip that
#'   is not the root is a child exactly once.
#'
#' @param tips_and_node_vector Numeric vector. Length = Ntip(tr) + Nnode(tr).
#' @param tr Phylo.
#'
#' @return ordered_by_edges. Numeric vector. Length = Nedge(tr).
#' @noRd
reorder_tip_and_node_to_edge <- function(tips_and_node_vector, tr){
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)
  check_equal(length(tips_and_node_vector), sum(ape::Ntip(tr), ape::Nnode(tr)))
  check_is_number(tips_and_node_vector[1])

  # Function -------------------------------------------------------------------
  ordered_by_edges <- rep(NA, ape::Nedge(tr))
  for (i in 1:ape::Nedge(tr)) {
    ordered_by_edges[i] <- tips_and_node_vector[tr$edge[i, 2]]
  }

  # Return output --------------------------------------------------------------
  return(ordered_by_edges)
}# end reorder_tip_and_node_to_edge()
