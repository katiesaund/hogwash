discretize_confidence_using_threshold <- function(confidence_vector, threshold){
  # Function description -------------------------------------------------------
  # Given a vector with values that describe confidence, binarize the vector a
  # accoriding to a cutoff value.
  #
  # Input:
  # Confidence vector. Numeric vector.
  # Threshold. Number.
  #
  # Output:
  # Confidence vector. Binary vector.
  #
  # Check inputs ---------------------------------------------------------------
  check_is_number(threshold)
  if (!is.vector(confidence_vector)){stop("input must be a numeric vector")}
  check_is_number(confidence_vector[1])

  # Function -------------------------------------------------------------------
  confidence_vector[confidence_vector  < threshold] <- 0
  confidence_vector[confidence_vector >= threshold] <- 1

  # Check and return output ----------------------------------------------------
  check_if_binary_vector(confidence_vector)
  return(confidence_vector)
} # end discretize_confidence_using_threshold()


report_num_high_confidence_trans_edge <- function(genotype_transition,
                                                  high_conf_edges, geno_names){
  # Function description -------------------------------------------------------
  # Given a genotype for which you have: a list of vectors that stores if there
  # is a genotype transition or not for each edge (genotype_transition), a list
  # of vectors that stores if that edge is high confidence or not
  # (high_conf_edges), and a character vector of the genotype names -- create an
  # object that stores the number of high confidence transition edges per
  # genotype.
  #
  # Inputs:
  # genotype_transition. List of numeric vectors. Number of lists == number of
  #                      genotypes. Length of vector == Nedge(tr)
  # high_conf_edges. Binary vector. List of numeric vectors. Number of lists ==
  #                  number of genotypes. Length of vector == Nedge(tr)
  # geno_names. Character vector. Length == ncol(genotype_matrix)
  #
  # Outputs:
  # num_high_confidence_transition_edges. Numeric vector. Count of number of
  #                                       high confidence transitions per
  #                                       genotype. Vector is named with
  #                                       genotype names.
  #
  # Check input ----------------------------------------------------------------
  if (length(genotype_transition) != length(geno_names)){
    stop("genotype_transition should have one vector for each genotype")
  }
  if (length(high_conf_edges) != length(geno_names)){
    stop("high_conf_edges should have one vector for each genotype")
  }
  if(!is.vector(genotype_transition[[1]]$transition)){stop("Input must have a vector.")}
  if(!is.vector(high_conf_edges[[1]]$transition)){stop("Input must have a vector.")}
  if(!is.vector(geno_names)){stop("Input must be a vector.")}
  if(!is.character(geno_names[1])){stop("Input must be a character vector.")}

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  num_high_confidence_transition_edges <- rep(0, length(high_conf_edges))
  for (p in 1:length(high_conf_edges)){
    num_high_confidence_transition_edges[p] <- sum(genotype_transition[[p]]$transition * high_conf_edges[[p]])
  }

  names(num_high_confidence_transition_edges) <- geno_names
  return(num_high_confidence_transition_edges)
} # end report_num_high_confidence_trans_edge

assign_high_confidence_to_transition_edges <- function(tr,
                                                       all_confidence_by_edge,
                                                       genotype_transition_by_edges,
                                                       geno){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # tr. Phylo.
  # all_confidence_by_edge. ?
  # genotype_transition_by_edges. ?
  # geno. Matrix. Binary.
  #
  # Outputs:
  # edge_and_parent_confident_and_trans_edge. ?
  #
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  if (length(genotype_transition_by_edges[[1]]$transition) != Nedge(tr)){
    stop("Dimension mismatch")
  }
  check_for_root_and_bootstrap(tr)
  if (length(all_confidence_by_edge[[1]]) != Nedge(tr)){
    stop("Dimension mismatch")
  }
  # Function -------------------------------------------------------------------
  # Identify all edges for which the edge and the parent edge are both high confidence
  edge_and_parent_both_confident <- edge_and_parent_confident_and_trans_edge <- rep(list(rep(0, Nedge(tr))), ncol(geno))
  for (ge in 1:ncol(geno)){
    for (ed in 2:(Nedge(tr) - 1)){ # start at 2 because there isn't a parent edge to edge 1, end at Nedge- 1 because no parent to the last edge either
      parent_edge <- find_parent_edge(tr, ed)
      if (all_confidence_by_edge[[ge]][ed] == 1 & all_confidence_by_edge[[ge]][parent_edge] == 1){
        edge_and_parent_both_confident[[ge]][ed] <- 1
      }
    }
    edge_and_parent_both_confident[[ge]][1] <- all_confidence_by_edge[[ge]][1] # have to account for the fact that there isn't a parent edge to edge 1
    edge_and_parent_both_confident[[ge]][Nedge(tr)] <- all_confidence_by_edge[[ge]][Nedge(tr)] # have to account for the fact that there isn't a parent edge to last edge
  }

  # Identify high confidence transition edges by overlapping the above and transitions
  for (k in 1:ncol(geno)){
    edge_and_parent_confident_and_trans_edge[[k]] <- as.numeric((edge_and_parent_both_confident[[k]] + genotype_transition_by_edges[[k]]$transition) == 2)
  }

  # Return output --------------------------------------------------------------
  # Return that overlap as high confidence transitions
  return(edge_and_parent_confident_and_trans_edge)
} # end assign_high_confidence_to_transition_edges()
