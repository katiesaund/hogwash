# Note on phylogenetic tree structure
# tr$edge:
#       [,1] [,2]
# [1,]  111  112
# [2,]  112  113
# [3,]  113  114
# [4,]  114  115
# Dim = 2*(Ntip(tr)) x 2
# Each row is an edge.
# The first colum is the older node and the second column is the younger node.


#' Identify transition edges on a tree
#'
#' @description Given a reconstruction identify which edges on tree are
#'   transitions. A transition edge is one in which the parent node and child
#'   node differ.
#'
#' @param tr Phylo.
#' @param mat Matrix. Phenotype (either continuous or binary) or a binary
#'   genotype. Dim: nrow = Ntrip(tr) x ncol = {1 if phenotype or number of
#'   genotypes}.
#' @param num Integer. Index of current genotype (column number in genotype
#'   matrix).
#' @param node_recon Numeric vector. Either pheno_recon_and_conf$node_anc_rec or
#'   geno_recon_and_conf[[k]]$node_anc_rec. Length = Nnode(tr). Ancestral
#'   reconstruction value for each node.
#' @param disc_cont Character string. Either "discrete" or "continuous".
#'
#' @return List.
#'   \describe{
#'    \item{transition}{Continuous: NA, because this is meaningless for
#'    continuous data as all edges will be transitions. Discrete: Numeric vector
#'    of 0 or 1. 0 indicates parent and child node are identical. 1 indicates
#'    parent and child node differ.}
#'    \item{trans_dir}{Numeric Vector. Continuous and discrete: gives the
#'    direction of the transition. When parent < child or parent_0_child_1 value
#'    is +1. When parent > child or parent_1_child_0 value is -1.}
#'   }
#' @noRd
identify_transition_edges <- function(tr, mat, num, node_recon, disc_cont){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  check_is_number(num)
  check_str_is_discrete_or_continuous(disc_cont)

  # FUNCTION -------------------------------------------------------------------
  transition <- transition_direction <-
    parent_node <- child_node <- integer(ape::Nedge(tr))
  older <- 1 # older node is 1st column in tr$edge
  younger <- 2 # younger node is 2nd column in tr$edge
  parent_0_child_1 <- 1
  parent_1_child_0 <- -1
  parent_equals_child <- 0
  both_parent_and_child_are_one <- 2


  for (i in 1:ape::Nedge(tr)) {
    if (is_tip(tr$edge[i, older], tr)) {
      stop("tree invalid")
    }
    parent_node[i] <- node_recon[tr$edge[i, older] - ape::Ntip(tr)]
    if (is_tip(tr$edge[i, younger], tr)) {
      # child is a tip
      child_node[i]  <- mat[, num][tr$edge[i, younger]]
    } else {
      # child is internal nodes
      child_node[i]  <- node_recon[tr$edge[i, younger] - ape::Ntip(tr)]
    }

    transition[i] <- sum(parent_node[i] + child_node[i])
    # transition[i] is either 0, 1, or 2 for discrete traits
    # transition[i] is not to be used when the trait is continuous because all,
    # or very nearly all edges are transition edges.

    if (parent_node[i] > child_node[i]) {
      transition_direction[i] <- parent_1_child_0
    } else if (parent_node[i] < child_node[i]) {
      transition_direction[i] <- parent_0_child_1
    }

    if (disc_cont == "discrete") {
      transition[transition == both_parent_and_child_are_one] <-
        parent_equals_child # parent_node == child_node, then no transition.
    } else {
      transition <- NA # Not used when the trait is continuous.
    }
  }
  # Check and return output ----------------------------------------------------
  results <- list("transition" = transition, "trans_dir" = transition_direction)
  return(results)
}

#' Keep only genotypes with two or more high confidence transition edges
#'
#' @description Since we're looking for convergence of transitions we need a
#'   second quality control step where we remove genotypes that have only 1
#'   transition-edge or where the transition edges are identical!
#'
#' @param genotype_transition List of multiple vectors ($transition and
#'   $trans_dir). Length(list) = number of genotypes. Length(vector) =
#'   Nedge(tr).
#' @param genotype_confidence List of vectors. Length(list) = number of
#'   genotypes. Length(vector) = Nedge(tr). Binary.
#'
#' @return at_least_two_hi_conf_trans_ed Logical vector. Length ==
#'   length(genotype_transition) == length(genotype_confidence).
#' @noRd
keep_two_plus_hi_conf_tran_ed <- function(genotype_transition,
                                          genotype_confidence){
  # Check inputs ---------------------------------------------------------------
  check_equal(length(genotype_transition), length(genotype_confidence))
  if (!is.vector(genotype_transition[[1]]$transition)) {
    stop("Input must be a numeric vector")
  }
  check_if_binary_vector(genotype_confidence[[1]])

  # Function -------------------------------------------------------------------
  at_least_two_hi_conf_trans_ed <-
    rep(FALSE, length(genotype_transition))
  for (p in 1:length(genotype_transition)) {
    if (sum(genotype_transition[[p]]$transition *
            genotype_confidence[[p]]) > 1) {
      at_least_two_hi_conf_trans_ed[p] <- TRUE
    }
  }

  # Check and return output ----------------------------------------------------
  return(at_least_two_hi_conf_trans_ed)
}

#' Prepare genotype transition object for PhyC Test
#'
#' @description Discrete testing requires two different definitions of genotype
#'   transition, one for PhyC and one for Synchronous Test. This function
#'   converts geno_trans$transition from the object created for synchronous test
#'   to the version required for the PhyC test.
#'
#' @details For PhyC prepare genotype transition as below: keep only WT ->
#'   mutant transitions (0 -> 1)
#'
#' @param geno Matrix. Nrow = Ntip(Tree). Ncol = number of genotypes. Binary.
#' @param genotype_transition List of lists. Each sublist has two vectors:
#'   \describe{
#'     \item{transition}{Length == Nedge(tree). 0/1}
#'     \item{trans_dir}{-1/0/1. Length == Nedge(tree).}
#'   }
#'
#' @return genotype_transition. List with $transition and $trans_dir.
#' @noRd
prep_geno_trans_for_phyc <- function(geno, genotype_transition){
  # Check inputs ---------------------------------------------------------------
  check_if_binary_matrix(geno)
  check_equal(length(genotype_transition), ncol(geno))
  check_if_binary_vector(genotype_transition[[1]]$transition)

  # Function -------------------------------------------------------------------
  for (k in 1:ncol(geno)) {
    parent_WT_child_mutant <- 1 # 1 implies parent < child
    parent_mutant_child_WT <- -1 # -1 implies parent > child
    no_transition <- 0 # 0 implies parent == child
    genotype_transition[[k]]$transition <-
      as.numeric(genotype_transition[[k]]$trans_dir == parent_WT_child_mutant)
    genotype_transition[[k]]$trans_dir[genotype_transition[[k]]$trans_dir
                                       == parent_mutant_child_WT] <-
      no_transition # erase transitions in the opposite direction
  }

  # Return output --------------------------------------------------------------
  return(genotype_transition)
}
