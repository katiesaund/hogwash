identify_transition_edges <- function(tr, mat, num, node_recon, disc_cont){
  # TODO-- it's very imporant that I come back to this and not rely solely on this function for the three distinct tests!
  # TODO for example with original phyc test -- I need to ignore any transitions with transition direction -1 and only look at transition direction +1.
  # Function description -------------------------------------------------------
  # Given a reconstruction identify which edges on tree are transitions.
  # A transition edge is one in which the parent node and child node differ.
  #
  # Inputs:
  # tr: phylogenetic tree.
  # mat: Matrix. Phenotype (either continuous or binary) and a binary genotype.
  # num: numeric. Column index (genotype index).
  # node_recon. Either pheno_recon_and_conf$node_anc_rec or geno_recon_and_conf[[k]]$node_anc_rec.
  # disc_cont: string. Either "discrete" or "continuous".
  #
  # Outputs:
  # List.
  #   $transition:
  #     Continuous: NA, because this is meaningless for continuous data (all edges will be transitions)
  #     Discrete: Numeric vector of 0 or 1. 0 indicates parent and child node are identical. 1 indicates parent and child node differ.
  #   $trans_dir: Numeric Vector.
  #     Continuous and discrete: gives the direction of the transition.
  #     When parent < child or parent_0_child_1 value is +1.
  #     When parent > child or parent_1_child_0 value is -1.

  # Check input ---------------------------------------------------------------#
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  check_is_number(num)
  check_str_is_discrete_or_continuous(disc_cont)

  # FUNCTION ------------------------------------------------------------------#
  transition <- transition_direction <- parent_node <- child_node <- integer(Nedge(tr)) # initialize all as zeroes
  older <- 1 # older node is 1st column in tr$edge
  younger <- 2 # younger node is 2nd column in tr$edge
  parent_0_child_1 <- 1
  parent_1_child_0 <- -1
  parent_equals_child <- 0
  both_parent_and_child_are_one <- 2


  for (i in 1:Nedge(tr)){
    if (is_tip(tr$edge[i, older], tr)){
      stop("tree invalid")
    }
    parent_node[i] <- node_recon[tr$edge[i, older] - Ntip(tr)]   # Assign node value
    if (is_tip(tr$edge[i, younger], tr)){ # child is a tip
      child_node[i]  <- mat[ , num][tr$edge[i, younger]]           # Assign tip value
    } else {  # child is internal nodes
      child_node[i]  <- node_recon[tr$edge[i, younger] - Ntip(tr)] # Assign node value
    }

    transition[i] <- sum(parent_node[i] + child_node[i])
    # transition[i] is either 0, 1, or 2 for discrete traits
    # transition[i] is not to be used when the trait is continuous because all, or very nearly all edges are transition edges.

    if (parent_node[i] > child_node[i]){
      transition_direction[i] <- parent_1_child_0
    } else if (parent_node[i] < child_node[i]){
      transition_direction[i] <- parent_0_child_1
    }

    if (disc_cont == "discrete"){
      transition[transition == both_parent_and_child_are_one] <- parent_equals_child # If parent_node and child_node had same value (1) then no transition occured.
    } else {
      transition <- NA # transition[i] is not to be used when the trait is continuous because all, or very nearly all edges are transition edges.
    }
  }
  # Check and return output ----------------------------------------------------
  results <- list("transition" = transition, "trans_dir" = transition_direction)
  return(results)
} # end identify_transition_edges()

# TODO the accuracy of this function assumes that the genotype_confidence input is already the transition confidence as calculated by: assign_high_confidence_to_transition_edges
# TODO how to I add a check to make sure that's true? ^^^^
keep_at_least_two_high_conf_trans_edges <- function(genotype_transition, genotype_confidence){
  # Function description -------------------------------------------------------
  # Since we're looking for convergence of transitions we need a second quality
  # control step where we remove genotypes that have only 1 transition-edge or
  # where the transition edges are identical!
  #
  # Inputs:
  # genotype_transition. List of multiple vectors ($transition and $trans_dir).
  #                      Length(list) = number of genotypes. Length(vector) =
  #                      Nedge(tr).
  # genotype_confidence. List of vectors. Length(list) = number of genotypes.
  #                      Length(vector) = Nedge(tr). Binary.
  #
  # Output:
  # has_at_least_two_high_confidence_transition_edges. Logical vector. Length ==
  #                                                    length(genotype_transition)
  #                                                    == lenght(genotype_confidence)
  #
  # Check inputs ---------------------------------------------------------------
  if (length(genotype_transition) != length(genotype_confidence)){
    stop("Both transition and confidence should have length corresponding to number of genotypes.")
  }
  if(!is.vector(genotype_transition[[1]]$transition)){
    stop("Genotype transition should have a vector called 'transition'.")
  }
  check_if_binary_vector(genotype_confidence[[1]])

  # Function -------------------------------------------------------------------
  has_at_least_two_high_confidence_transition_edges <- rep(FALSE, length(genotype_transition))
  for (p in 1:length(genotype_transition)){
    if (sum(genotype_transition[[p]]$transition * genotype_confidence[[p]]) > 1){
      has_at_least_two_high_confidence_transition_edges[p] <- TRUE
    }
  }

  # Check and return output ----------------------------------------------------
  return(has_at_least_two_high_confidence_transition_edges)
} # end keep_at_least_two_high_conf_trans_edges()

keep_hits_with_more_change_on_trans_edges <- function(results, pvals, fdr){
  # Function description -------------------------------------------------------
  # Of all of the significant hits, keep only those where the
  # median(delta phenotype) on transition edges is > median(delta phenotype) on
  # all edges. Do this so that we're only selecting for genotypes that are
  # plausibly have large effect on phenotype, rather than miniscule effects on
  # phenotype. Essentially, we're enforcing a one tailed test instead of using
  # the results of a two-tailed test.
  #
  # Inputs:
  # results. List of 8.
  #          $pvals. Named numeric vector. Length == number of genotypes. Values between 1 and 0. Names are genotype names.
  #          $ks_statistics. List of numeric vectors. Length of list == number of genotypes. Each vector has length == number of permutations. Values between 1 and 0.
  #          $observed_pheno_trans_delta. List of numeric vectors. Length of list == number of genotypes. Vectors are of variable length because length is the number of transition edges for that particular genotype. Vectors are numeric.
  #          $observed_pheno_non_trans_delta. List of numeric vectors. Length of list == number of genotypes. Vectors are of variable length because length is the number of non-transition edges for that particular genotype. Vectors are numeric.
  #          $trans_median. Numberic. Vector. Length = number of genotypes. Describes median delta phenotype on all transition edges.
  #          $all_edges_median. Numeric vector. Length = number of genotypes. Describes median delta phenotype on all edges.
  #          $num_genotypes. Integer. The number of genotypes.
  #          $observed_ks_stat. Numeric Vector. Length = number of genotypes. Values between 1 and 0.
  # pvals.  List of 2.
  #         $hit_pvals. Dataframe. 1 column. Nrow = number of genotypes. Row.names = genotypes. Column name = "fdr_corrected_pvals". Values between 1 and 0.
  #         $sig_pvals. Dataframe. 1 column. Nrow = number of genotypes that are significant after FDR correction. Column name = "fdr_corrected_pvals[fdr_corrected_pvals < alpha]". Row.names = genotypes. Nrow = is variable-- could be between 0 and max number of genotypes. It will only have rows if the corrected p-value is less than the alpha value.
  # fdr.    Number. False discovery rate (significance threshold). Between 0 and 1.
  # Output:
  # hits. Data.frame. 1 column. Colnum names = "fdr_corrected_pvals". Nrow = variable. Number of genotypes that are (1) significant after multiple test correction and (2) have higher median delta phenotype on transition edges than on all edges. Values are between 1 and 0. Rownames are genotypes.
  #
  # Check inputs ---------------------------------------------------------------
  check_num_between_0_and_1(fdr)

  # Function -------------------------------------------------------------------
  # Keep only hits with transition edge median delta phenotype higher than all edge delta phenotype median.
  temp <- pvals$hit_pvals[(results$trans_median > results$all_edges_median), , drop = FALSE]

  # Keep only those that are also significant after FDR correction.
  hits <- temp[temp[ , 1] < fdr, , drop = FALSE]

  # Check and return output ----------------------------------------------------
  return(hits)
} # end keep_hits_with_more_change_on_trans_edges()

prepare_genotype_transitions_for_original_discrete_test <- function(geno, genotype_transition){
  # Function description -------------------------------------------------------
  # Discrete testing requires two different definitions of genotype
  # transition:
  # 1) Version based on the original phyC
  # 2) Version for requiring concomitant transition state change phenotype
  #    and genotype.
  # This function converstions geno_trans$transition from the object created
  # for concomitant test to the version required for the original test.
  #
  # Note: for original phyC prepare genotype transition as below: keep only
  #       WT -> mutant transitions (0 -> 1).
  #
  # Inputs ---------------------------------------------------------------------
  # geno. Matrix. Nrow = Ntip(Tree). Ncol = number of genotypes. Binary.
  # genotype_transition. List of lists. Each sublist has two vectors.
  #       1) $transition. Length  == Nedge(tree). 0/1
  #       2) $trans_dir. -1/0/1. Length == Nedge(tree).
  #
  # Check inputs ---------------------------------------------------------------
  check_if_binary_matrix(geno)
  if (length(genotype_transition) != ncol(geno)){
    stop("Must have transition information for each genotype")
  }
  # Function -------------------------------------------------------------------
  for (k in 1:ncol(geno)){
    parent_WT_child_mutant <- 1 # 1 implies parent < child, -1 implies parent > child, 0 implies parent == child
    genotype_transition[[k]]$transition <- as.numeric(genotype_transition[[k]]$trans_dir == parent_WT_child_mutant)
  }

  # Return output --------------------------------------------------------------
  return(genotype_transition)
} # prepare_genotype_transitions_for_original_discrete_test
