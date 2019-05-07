# Katie Saund

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

# FUNCTIONS FOR PHYC ----------------------------------------------------------#
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

identify_transition_edges <- function(tr, mat, num, node_recon, disc_cont){
  # Function description -------------------------------------------------------
  # Given a reconstruction identify which edges on tree are transitions.
  # A transition edge is one in which the parent node and child node differ.
  #
  # Inputs:
  # tr: phylogenetic tree.
  # mat: Matrix.
  # num: numeric.
  # node_recon.
  # disc_cont: string. Either "discrete" or "continuous".
  #
  # Outputs:
  # List.
  #   $transition:
  #     Continuous: NA, because this is meaningless for continuous data (all edges will be transitions)
  #     Discrete: Numeric vector of 0 or 1. 0 indicates parent and child node are identical. 1 indicates parent and child node differ.
  #   $trans_dir: Numeric Vector.
  #     Continuous and discrete: gives the direction of the transition.

  # Check input ---------------------------------------------------------------#
  # TODO ADD CHECKS


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
    } else if (disc_cont == "continuous"){
      transition <- NA # transition[i] is not to be used when the trait is continuous because all, or very nearly all edges are transition edges.
    } else {
      stop("disc_cont must be discrete or continuous")
    }

  }
  # Check and return output ----------------------------------------------------
  results <- list("transition" = transition, "trans_dir" = transition_direction)
  return(results)
} # end identify_transition_edges()

assign_pheno_type <- function(mat){
  # Function description -------------------------------------------------------
  # Determine if the matrix is discrete or continuous.
  #
  # Input:
  # mat: Matrix. Should be a phenotype with 1 column.
  #
  # Outputs:
  # type. Character. Either "discrete" or "continuous".

  # Check input ----------------------------------------------------------------
  check_dimensions(mat, NULL, 1, 1, 1)

  # Function -------------------------------------------------------------------
  type <- "discrete"
  if (sum(!(mat %in% c(0, 1))) > 0){
    type <- "continuous"
  }

  # Check and return output ----------------------------------------------------
  check_is_string(type)
  return(type)
} # end assign_pheno_type()

calculate_phenotype_change_on_edge <- function(edge_list, phenotype_by_edges){
  # Function description -------------------------------------------------------
  # Quantify absoluate value of phenotype change on each tree edge.
  #
  # Input:
  # edge_list.          Numeric vector. Each number is the index of the tree edge to be used
  # phenotype_by_edges. Mat.            Dimensions: Nedge x 2 matrix. Entries are the phenotype value at the node, where row is the edge, 1st column is the parent node and 2nd column is the child node.
  #
  # Outputs:
  # delta.              Numeric vector.   Length = length(edge_list).

  # Check input ----------------------------------------------------------------
  if (max(edge_list) > nrow(phenotype_by_edges)){
    stop("Cannot calculate phenotype change on edges.")
  }
  check_dimensions(phenotype_by_edges, exact_rows = NULL, min_rows = max(edge_list), exact_cols = NULL, min_cols = 2)

  # Function -------------------------------------------------------------------
  delta <- rep(NA, length(unique(edge_list)))
  for (j in 1:length(edge_list)){
    delta[j] <- abs(phenotype_by_edges[edge_list[j], 1] - phenotype_by_edges[edge_list[j], 2])
  }

  # Check and return output ----------------------------------------------------
  return(delta)
} # end calculate_phenotype_change_on_edge()

run_ks_test <- function(t_index, non_t_index, phenotype_by_edges){
  # Function description -------------------------------------------------------
  # Run a Kolmogorov-Smirnov test on the continuous phenotype. The phenotype is
  # divided into two groups: transition edges and non-transition edges.
  #
  # Input:
  # t_index.     Vector. Each number is the index phenotype transition edges on the tree.
  # non_t_index. Vector.  Each number is the index phenotype transition edges on the tree.
  # phenotype_by_edges. Vector(?). Continuous phenotype.
  #
  # Outputs:
  # A list of 4 elements:
  # $pval. Numeric. KS Test p-value.
  # $statistic. Numeric. KS Test statistic.
  # $pheno_trans_delta. Numeric vector. The value of the delta of the phenotype on transition edges.
  # $pheno_non_trans_delta. Numeric vector. The value of the delta of the phenotype on all non-transition edges.

  # Check input ----------------------------------------------------------------


  # TODO add unit tests
  # TODO deal with cases when there isn't enough data (aren't at least 1 of least t_index and non_t_index)
  # t_index = transition edges
  # not_t_index = non transition edges
  # phenotype_by_edges = change in phenotype on each edge

  # Check input ---------------------------------------------------------------#
  if (length(t_index) < 1 | length(non_t_index) < 1){
    stop("Not enough high confidence transition edges to use for KS test.")
  }

  # Function -------------------------------------------------------------------

  # subset index to only high confidence edges:
  p_trans_delta     <- calculate_phenotype_change_on_edge(t_index,     phenotype_by_edges)
  p_non_trans_delta <- calculate_phenotype_change_on_edge(non_t_index, phenotype_by_edges)

  ks_results        <- ks.test(p_trans_delta, p_non_trans_delta)

  # Return output --------------------------------------------------------------

  results <- list("pval"      = round(ks_results$p.value, digits = 20),
                  "statistic" = round(ks_results$statistic, digits = 20),
                  "pheno_trans_delta"     = p_trans_delta,
                  "pheno_non_trans_delta" = p_non_trans_delta)
  return(results)
} # end run_ks_test()


calculate_genotype_significance <- function(mat, permutations, genotype_transition_list, tr, pheno_recon_ordered_by_edges, genotype_confidence, genotype_reconstruction){
  # Function description -------------------------------------------------------
  #
  # Input:
  # mat: Matrix. Should be a phenotype with 1 column.
  # permutations,
  # genotype_transition_list
  # tr
  # pheno_recon_ordered_by_edges
  # genotype_confidence
  # genotype_reconstruction
  #
  # Outputs:
  # pvals
  # ks_statistics" = empirical_ks_pval_list
  # "observed_pheno_trans_delta" = observed_pheno_trans_delta
  # "observed_pheno_non_trans_delta" = observed_pheno_non_trans_delta
  # "trans_median" = trans_median
  # "all_edges_median" = all_edges_median,
  # "num_genotypes" = ncol(mat))
  #
  # mat is a genotype_matrix. Each row is a genotype. Each column is a 0/1 value on an edge.
  # permutations is a numeric integer.
  # genotype_transition_list is a list of length = ncol(mat) and each list has length = Nedge(tr). All entries are 0 or 1.
  # tr object of phylo.
  # phenotype_recon_ordered_by_edges is a matrix. Nrow = Nedge(tr). Ncol =2. where the phenotype reconstruction is ordered by edges. IN THE SAME FORMAT AS tr$EDGE. SO NODE VALUES WILL APPEAR MULTIPLE TIMES IN THE tr.
  #            THIS FORMAT WILL MAKE IT MUCH EASIER TO CALCULATE PHENOTYPE CHANGE ON EDGES.
  # genotype_confidence is a list of length = ncol(mat) and each list has length = Nedge(tr). All entries are 0 (low confidence) or 1 (high confidence).
  #            NOTE: genotype_confidence lists the confidence in each edge. High confidence means the edge is high confidence by genotype reconstruction, phenotype reconstruction, bootstrap value, and edge length.

  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_permutation_num_valid(permutations)

  # Function -------------------------------------------------------------------
  num_genotypes <- ncol(mat)
  pvals <- observed_ks_pval <- trans_median <- all_edges_median <- observed_ks_stat <- rep(NA, num_genotypes) # 2018-11-12 added observed_ks_stat
  names(observed_ks_pval) <- names(pvals) <- colnames(mat)
  empirical_ks_pval_list <- empirical_ks_stat_list <- observed_pheno_trans_delta <- observed_pheno_non_trans_delta <- rep(list(0), num_genotypes) # 2018-11-12 added empirical_ks_stat_list

  for (i in 1:num_genotypes){
    # GRAB THE IDS OF THE TRANSITION EDGES:
    trans_index     <- c(1:Nedge(tr))[as.logical(genotype_transition_list[[i]]$transition)]
    non_trans_index <- c(1:Nedge(tr))[!genotype_transition_list[[i]]$transition]
    # [1] EX:  8 12 13 16 19 26 27 31 37 44 52 56 64 67 68 76 77 80 89 92 97 98
    # THESE EDGES ARE DEFINED BY THE NODES IN THE CORRESPONDING ROWS OF tr$EDGE

    # SUBSET TRANSITION / NON-TRANSITION EDGES TO ONLY HIGH CONFIDENCE ONES
    hi_conf_trans_index     <- trans_index[    as.logical(genotype_confidence[[i]][trans_index    ])]
    hi_conf_non_trans_index <- non_trans_index[as.logical(genotype_confidence[[i]][non_trans_index])]

    # Run KS test to find out if the phenotype change on transition edges is significantly different from phenotype change on non-transition edges
    observed_results <- run_ks_test(hi_conf_trans_index, hi_conf_non_trans_index, pheno_recon_ordered_by_edges)
    observed_ks_pval[i] <- observed_results$pval
    observed_ks_stat[i] <- observed_results$statistic

    # save these for reporting / plots
    observed_pheno_trans_delta[[i]]     <- observed_results$pheno_trans_delta
    observed_pheno_non_trans_delta[[i]] <- observed_results$pheno_non_trans_delta
    trans_median[i]     <- median(observed_results$pheno_trans_delta)
    all_edges_median[i] <- median(c(observed_results$pheno_trans_delta, observed_results$pheno_non_trans_delta))
    #

    # do the permutation part
    num_sample           <- length(hi_conf_trans_index)
    all_edges            <- c(1:Nedge(tr))
    which_branches       <- all_edges[as.logical(genotype_confidence[[i]])]
    all_sampled_branches <- matrix(   nrow = permutations, ncol = num_sample)
    redistributed_hits   <- matrix(0, nrow = permutations, ncol = length(which_branches))

    # ii. Create a matrix where each row is a random sampling of the high
    #     confidence edges of the tr where probability of choosing edges is
    #     proportional to length of edge. Number of edges selected for the
    #     permuted data set is the number of times the empirical genotype
    #     appears.

    set.seed(1) # for reproducability of the sample() function
    for (j in 1:permutations){  # create a random sample of the tr
      curr_num_branch <- length(which_branches)
      all_sampled_branches[j, ] <- sample(1:curr_num_branch,
                                          size = num_sample,
                                          replace = TRUE,
                                          prob = tr$edge.length[which_branches]/sum(tr$edge.length[which_branches]))
    } # end for (j)

    # all_sampled_branches is my new, permuted "hi_conf_trans_index" where each row is a new list of transition genotype branches
    # BUT CAVEAT: these are just fake/null transitions and some of them are probably actually touching! If I wanted to be
    # Super legit I would recreate as many hits, calculate new transitions, and then use those in my permutation test, somehow
    # controlling for variable numbers of transitions. But not doing that for now.

    # calculate permuted pheno deltas
    empirical_ks_pval <- empirical_ks_stat <- rep(NA, permutations)
    for (k in 1:nrow(all_sampled_branches)){
      permuted_trans_index     <- unique(all_sampled_branches[k, ])
      permuted_non_trans_index <- c(1:length(which_branches))[!(c(1:length(which_branches)) %in% unique(all_sampled_branches[k, ]))]
      empirical_results        <- run_ks_test(permuted_trans_index, permuted_non_trans_index, pheno_recon_ordered_by_edges)
      empirical_ks_pval[k]     <- empirical_results$pval
      empirical_ks_stat[k]     <- empirical_results$statistic
    } # end for (k)

    # the observed ks.test statistic: observed_ks_stat[i] (fixed this 2018-11-12; before I had it wrong with observed_ks_pval)
    # empirical p value caluclation here: (1 + more extreme observations) / (1 + permutations)
    empirical_ks_pval_list[[i]] <- empirical_ks_pval
    empirical_ks_stat_list[[i]] <- empirical_ks_stat
    pvals[i] <- (sum(1 + sum(empirical_ks_stat > observed_ks_stat[i]))/(permutations + 1))
  } # end for (i)

  # Return output --------------------------------------------------------------
  results <- list("pvals" = pvals, "ks_statistics" = empirical_ks_stat_list,
                  "observed_pheno_trans_delta" = observed_pheno_trans_delta,
                  "observed_pheno_non_trans_delta" = observed_pheno_non_trans_delta,
                  "trans_median" = trans_median, "all_edges_median" = all_edges_median,
                  "num_genotypes" = ncol(mat), "observed_ks_stat" = observed_ks_stat) # 2018-11-28
  return(results)
} # end calculate_genotype_significance()

convert_matrix_to_vector <- function(mat){
  # Function description -------------------------------------------------------
  # Convert a single column matrix into a vector, retain row names as names of vector.
  #
  # Input:
  # mat. Matrix. Matrix should have only 1 column.
  #
  # Output:
  # vec. Vector.
  #
  # Check input ----------------------------------------------------------------
  check_dimensions(mat, NULL, 1, 1, 1)

  # Function -------------------------------------------------------------------
  vec <- as.vector(unlist(mat[ , 1]))
  names(vec) <- row.names(mat)

  # Check and return output ----------------------------------------------------
  check_if_vector(vec)
  return(vec)
} # end convert_matrix_to_vector()

create_file_name <- function(output_dir, output_name, other_info){
  # Function description -------------------------------------------------------
  # Create a file name string that includes the path to the output directory.
  #
  # Input:
  # output_dir.  Character.
  # output_name. Character.
  # other_info.  Character.
  #
  # Output:
  # file_name. Character.
  #
  # Check input ----------------------------------------------------------------
  check_if_dir_exists(output_dir)
  check_is_string(output_name)
  check_is_string(other_info)

  # Function -------------------------------------------------------------------
  file_name <- paste(output_dir, "/", "phyc_", output_name, "_", other_info, sep = "")

  # Check and return output ----------------------------------------------------
  check_is_string(file_name)
  return(file_name)
} # end create_file_name()

create_test_data <- function(){
  # Function description -------------------------------------------------------
  # Create a set of reproducible test data.
  #
  # Input:
  # None.
  #
  # Output:
  # tree.             Phylo.
  # phenotype_matrix. Matrix.
  # genotype_matrix.  Matrix.
  #
  # Function -------------------------------------------------------------------
  # Create tree
  set.seed(1)
  tips            <- 50
  tree            <- rtree(n = tips, rooted = TRUE)
  tree$node.label <- rtruncnorm(n = Nnode(tree), sd = 10, mean = 85, a = 0, b = 100) # dummy tree bootstrap values

  # Create continous phenotype
  phenotype_matrix <- as.matrix(fastBM(tree))

  # Create genotypes
  genotype_matrix <- matrix(NA, nrow = tips, ncol = 100)
  for (i in 1:ncol(genotype_matrix)){
    genotype_matrix[ , i] <- rbinom(tips, 1, 0.5)
  }
  row.names(genotype_matrix)  <- tree$tip.label
  colnames(genotype_matrix) <- paste("snp", c(1:100), sep = "_")

  # Check and return output ----------------------------------------------------
  check_for_root_and_bootstrap(tree)
  check_dimensions(phenotype_matrix, Ntip(tree), 2, 1, 1)
  check_dimensions(genotype_matrix, Ntip(tree), 2, NULL, 1)
  check_if_binary_matrix(genotype_matrix)
  check_rownames(phenotype_matrix, tree)
  check_rownames(genotype_matrix, tree)

  results <- list("tree" = tree, "phenotype" = phenotype_matrix, "genotype" = genotype_matrix)
  return(results)
} # end create_test_data()

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
  # Function -------------------------------------------------------------------
  confidence_vector[confidence_vector  < threshold] <- 0
  confidence_vector[confidence_vector >= threshold] <- 1

  # Check and return output ----------------------------------------------------
  check_if_binary_vector(confidence_vector)
  return(confidence_vector)
} # end discretize_confidence_using_threshold()

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

reorder_tips_and_nodes_to_edges <- function(tips_and_node_vector, tr){
  # Function description -------------------------------------------------------
  # TODO ??
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

get_sig_hits_while_correcting_for_multiple_testing <- function(hit_values, alpha){
  # Function description -------------------------------------------------------
  # get_sig_hits_while_correcting_for_multiple_testing creates 2 vectors:
  # 1) vector of all FDR corrected pvalues, &
  # 2) vector of just significant hits.
  #
  # Inputs:
  # hit_values: Character. Vector of the empirical p-value for each genotype. Each entry corresponds to respective column in genotype_matrix. Should be between 0 and 1.
  # alpha: significance threshold as input by the user.
  #
  # Output:
  # fdr_corrected_pvals. Vector. all FDR corrected pvalues.
  # sig_pvals.           Vector. Only the significant FDR corrected pvalues.
  #
  # Check inputs ---------------------------------------------------------------
  check_if_alpha_valid(alpha)

  # Function -------------------------------------------------------------------
  fdr_corrected_pvals            <- p.adjust(hit_values, method = "fdr")
  sig_pvals                      <- as.data.frame(fdr_corrected_pvals[fdr_corrected_pvals < alpha])
  sig_pvals                      <- as.data.frame(sig_pvals)
  row.names(sig_pvals)           <- names(hit_values)[fdr_corrected_pvals < alpha]
  fdr_corrected_pvals            <- as.data.frame(fdr_corrected_pvals)
  row.names(fdr_corrected_pvals) <- names(hit_values)

  # Check and return output ----------------------------------------------------
  results                        <- list("hit_pvals" = fdr_corrected_pvals, "sig_pvals" = sig_pvals)
  return(results)
} # end get_sig_hits_while_correcting_for_multiple_testing

keep_at_least_two_high_conf_trans_edges <- function(genotype_transition, genotype_confidence){
  # TODO
  # 1) Update description of inputs
  # 2) check inputs.
  # 3) Update description of output.
  # 4) Check output.
  #
  # Function description -------------------------------------------------------
  # Since we're looking for convergence of transitions we need a second quality
  # control step where we remove genotypes that have only 1 transition-edge or
  # where the transition edges are identical!
  #
  # Inputs:
  # genotype_transition.
  # genotype_confidence.
  #
  # Output:
  # has_at_least_two_high_confidence_transition_edges.
  #
  # Check inputs ---------------------------------------------------------------

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

keep_hits_with_more_change_on_trans_edges <- function(results, pvals, a){
  # TODO
  # 1) Update description of inputs
  # 2) check inputs.
  # 3) Update description of output.
  # 4) Check output.
  #
  # Function description -------------------------------------------------------
  # SUBSET SIGNIFICANT HITS WHERE THE MEDIAN(DELTA PHENOTYPE) ON TRANSITION EDGES IS > MEDIAN(DELTA PHENOTYPE) ON ALL EDGES
  #
  # Inputs:
  # results.
  # pvals.
  # a.      Number. Alpha (significance threshold).
  # Output:
  # has_at_least_two_high_confidence_transition_edges.
  #
  # Check inputs ---------------------------------------------------------------
  check_if_alpha_valid(a)

  # Function -------------------------------------------------------------------
  temp <- pvals$hit_pvals[(results$trans_median > results$all_edges_median), , drop = FALSE]
  hits <- temp[temp[ , 1] < a, , drop = FALSE]

  # Check and return output ----------------------------------------------------
  return(hits)
} # end keep_hits_with_more_change_on_trans_edges()

calc_raw_diff <- function(edge_list, ph_edges){
  # Function description ------------------------------------------------------#
  # Subtract child node phenotype value from parent node phenotype value.
  #
  # Inputs:
  # TODO
  # edge_list. ?
  # ph_edges. ?
  #
  # Output:
  # delta. ?
  #
  # Check inputs ---------------------------------------------------------------
  # TODO
  # Function -------------------------------------------------------------------
  delta <- rep(NA, length(unique(edge_list)))
  for (j in 1:length(edge_list)){
    delta[j] <- ph_edges[edge_list[j], 1] - ph_edges[edge_list[j], 2]
  }
  # Return outputs -------------------------------------------------------------
  return(delta)
} # end calc_raw_diff()

reduce_redundancy <- function(mat, tr, dir, name){
  # Function description -------------------------------------------------------
  # This function REMOVES GENOTYPES THAT ARE: RARE, TOO COMMON, OR IDENTICAL
  #
  # Inputs:
  # mat.  Character. Path to matrix file.
  # tr.   Phylo.
  # dir.  Character.
  # name. Character.
  #
  # Output:
  # mat. Matrix.
  #
  # Check inputs ---------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_if_binary_matrix(mat)

  # Function -------------------------------------------------------------------
  geno_to_drop <- rep(FALSE, ncol(mat))
  geno_to_drop[colSums(mat) <= 1] <- TRUE
  geno_to_drop[colSums(mat) == (Ntip(tr) - 1)] <- TRUE
  geno_to_drop[colSums(mat) == Ntip(tr)] <- TRUE
  dropped_genotype_names <- colnames(mat)[geno_to_drop]
  mat <- mat[ , !geno_to_drop, drop = FALSE]

  # Check and return output ----------------------------------------------------
  check_if_binary_matrix(mat)
  check_rownames(mat, tr)
  results <- list("mat" = mat, "dropped_genotype_names" = dropped_genotype_names)
  return(results)
} # end reduce_redundancy()



# DISCRETE PHYC LIBRARY -------------------------------------------------------#

count_hits_on_edges <- function(genotype_reconstruction, phenotype_reconstruction, high_confidence_edges, phenotype_confidence){
  # TODO: update description
  # 1) We're now workingwith ordered by edges rather than ordered by tips then nodes
  # 2) rename combined/genotype confidence. it's confusing.
  # 3) stop returning combined_confidence- don't need it.

  # genotype_reconstruction is a list of vectors. Each vector corresponds with 1 genotype from the geno_mat.
  # Each entry in the vector corresponds to edge
  # phenotype_reconstruction is a vector where each entry corresponds to a node or tip.
  # genotype_confidence is a list of vectors. Each vector corresponds with 1 genotype from the geno_mat.
  # Each entry in the vector corresponds to a node or tip. 1 means high confidence in the node, 0 means low confidence.
  # phenotype_confidence is a vector where each entry corresponds to a node of tip. Same encoding as above.

  # TODO add checks / validations
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------


  # FUNCTION ------------------------------------------------------------------#

  both_present <- sapply(1:length(high_confidence_edges), function(x) {
    sum(phenotype_reconstruction[as.logical(high_confidence_edges[[x]])] + genotype_reconstruction[[x]][as.logical(high_confidence_edges[[x]])] == 2)
  })

  only_geno_present <- sapply(1:length(high_confidence_edges), function(x) {
    sum(genotype_reconstruction[[x]][as.logical(high_confidence_edges[[x]])]) - both_present[x]
  })

  hit_counts <- list("both_present" = both_present, "only_geno_present" = only_geno_present)
  return(hit_counts)
} #end count_hits_on_edges()

calculate_hit_pvals_corrected <- function(hit_counts, phenotype_reconstruction, tr, mat, permutations, alpha, high_confidence_edges){
  # calculate_genotype_pvals is the "meat" of the phyC algorithm.
  # calculate_genotype_pvals returns the empirical p-value for each observed genotype.
  # Algorithm overview:
  # 1. Subset tree edges to those with high confidence (as determined by phenotype & genotype reconstructions as well as tree bootstrap values).
  # 2. For each genotype from the genotype_matrix:
  #    i.   Exclude any edges on which the genotype is not present
  #    ii.  Create a matrix where each row is a random sampling of the high confidence edges of the tree where probability of choosing edges is proportional to length of edge. Number of
  #         edges selected for the permuted data set is the number of times the empirical genotype appears.
  #    iii. Calculate when the randomly permuted genotype overlap with the high confidence phenotype (these values create the null distribution for the permuted genotypes)
  #    iv.  Calculate empirical p-values.
  #    v.   Plot observed values on the null distribution.

  # Inputs:
  # hit_counts.               List.    2 part list containing results of count_hits_on_edges
  # phenotype_reconstruction. Numeric. 0/1. Length of Nedge(tree)
  # tree:                     Phylo.   Rooted phylogenetic tree.
  # mat:                      Matrix.  The genotype matrix. All entries are either 0 or 1. No NAs.
  # permutations:             Number.  The number of permutations.
  # alpha.                    Number.  Significance threshold.

  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------

  # read in variables
  both_present      <- hit_counts$both_present
  only_geno_present <- hit_counts$only_geno_present

  # initialize some values
  num_sample                            <- both_present + only_geno_present
  distribution_of_hits_from_permutation <- rep(0, length(both_present))
  hit_pvals                             <- rep(NA, ncol(mat))
  num_branch                            <- sapply(high_confidence_edges, function(x) sum(x))
  all_edges                             <- c(1:Nedge(tr))
  which_branches                        <- list(rep(0, length(high_confidence_edges)))

  # 1. Subset tr edges to those with high confidence (as determined by phenotype & genotype reconstructions as well as tr bootstrap values).
  for (i in 1:length(high_confidence_edges)){
    if (length(all_edges) == length(high_confidence_edges[[i]])){
      which_branches[[i]] <- all_edges[as.logical(high_confidence_edges[[i]])]
    } else {
      stop("which_branches incorrect")
    }
  }

  # 2. For each genotype from the genotype_matrix:
  record_of_redistributed_both_present <- rep(list(0), length(both_present))
  for (i in 1:length(both_present)){ # looping over each genotype in the genotype_matrix
    if (num_sample[i] > num_branch[i]){
      stop("Too many hits on the branches")
    }
    if((both_present[[i]] + only_geno_present[[i]]) < 2){
      # If there are 1 or 0 high confidence edges with the genotype present then the p-value should be reported as 1.0;
      # both present and just hit present are made up of only high confidence branches as defined in count_hits_on_edges
      # which isn't quite true but will indicate that we cannot calculate a p-value because we cannot detect any
      # convergence with one or fewer affected branches. And that we should skip to the next genotype to run the permutation test.
      hit_pvals[i] <- 1.0
    } else {
      # initialize some counters/variables for this loop
      counter              <- 0
      resampling_record    <- numeric(0)
      all_sampled_branches <- matrix(nrow = permutations, ncol = num_sample[i])
      redistributed_hits   <- matrix(0, nrow = permutations, ncol = length(which_branches[[i]]))
      # create a random sample of the tr
      set.seed(1)
      for(j in 1:permutations){
        curr_num_branch <- num_branch[i]
        all_sampled_branches[j, ] <- sample(1:curr_num_branch,
                                            size = num_sample[i],
                                            replace = TRUE,
                                            prob = tr$edge.length[which_branches[[i]]]/sum(tr$edge.length[which_branches[[i]]]))
      }

      #this if statement deals with situation when num_sampled = 1
      if(nrow(all_sampled_branches) != permutations){
        all_sampled_branches <- t(all_sampled_branches)
      }

      for (m in 1:nrow(all_sampled_branches)){
        redistributed_hits[m, ][all_sampled_branches[m, ]] <- 1
      }
      # right now redistributed hits is simply a matrix marking which edges we hit during
      # the permutation step above- because the point of the permutation is to pretend as if
      # we're assigning the presence of the genotype to random edges in the tree.

      empirical_both_present <- sapply(1:nrow(redistributed_hits), function(x) {
        sum(phenotype_reconstruction[as.logical(high_confidence_edges[[i]])] + redistributed_hits[x, ] == 2)
      })

      empirical_only_geno_present <- sapply(1:nrow(redistributed_hits), function(x) {
        sum(redistributed_hits[x, ]) - empirical_both_present[x]
      })

      # now beginning calculating when the randomly permuted "genotypes" overlap with the hi confidence phenotype
      temp <- matrix(0, nrow = permutations, ncol = length(which_branches[[i]])) #used to be length(Nedge(tree))
      redistributed_both_present <- rep(0, permutations)
      two <- rep(2, length(which_branches[[i]]))

      for (p in 1:permutations){
        if (length(redistributed_hits[p, ]) == length(phenotype_reconstruction[as.logical(high_confidence_edges[[i]])])){
          temp[p, ] <- redistributed_hits[p, ] + phenotype_reconstruction[as.logical(high_confidence_edges[[i]])]
        } else {
          stop("Dimension mismatch redistributed_hits and phenotype_reconstruction.")
        }
      }

      # temp ranges from 0 to 2.
      # A value of 0 means that the edge (as indicated by column) was not selected in that permutation (as indicated by row) & it's a low confidence pheno edge
      # A value of 1 means that either the edge was not selected in that permutation but is a low confidence edge or vice versa
      # A value of 2 means that the edge was selected in that permutation and it's a high confidence edge.
      for (q in 1:permutations){
        redistributed_both_present[q] <- sum(temp[q, ] == two)
      }

      # x_on_r <- sum(empirical_both_present >= both_present[i])
      # y_on_s <- sum(empirical_only_geno_present <= only_geno_present[i])

      # count only times when permuted (empirical) overlap of genotype and phenotype is more common than obsered (both_present[i])
      # and when the permuted (empirical) genotype does not overlap with phenotype is less common than observed (only_geno_present[i])
      new_counter <- sum((empirical_both_present >= both_present[i]) * (empirical_only_geno_present <= only_geno_present[i]))
      temp_pval <- ((new_counter + 1)/(permutations + 1))
      # Katie TODO 2019-02-26: save the new counter for each permutation, then plot that with the real as a red bar, add p val on top.

      if (sort(redistributed_both_present, decreasing = FALSE)[(alpha * permutations)] == 0 & both_present[i] == 0){
        pval <- 1
      } else if (temp_pval == 0 | temp_pval == 1){
        pval <- 2/(permutations + 1)
      } else if (temp_pval > 0.5){
        pval <- ((1 - temp_pval)  * 2)
      } else if (temp_pval <= 0.5){
        pval <- (temp_pval * 2)
      }

      record_of_redistributed_both_present[[i]] <- redistributed_both_present
      hit_pvals[i] <- format(round(pval, 20), nsmall = 20)
    }
  }
  names(hit_pvals) <- colnames(mat)
  results <- list("hit_pvals" = hit_pvals, "permuted_count" = record_of_redistributed_both_present, "observed_overlap" = both_present)
  # return(hit_pvals)
  return(results)
} # end calculate_hit_pvals_corrected

get_dropped_genotypes <- function(geno, keepers){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  dropped_genotype_names <- colnames(geno)[!keepers]
  return(dropped_genotype_names)
} # end get_dropped_genotypes

report_num_high_confidence_trans_edge <- function(genotype_transition, high_conf_edges, geno_names){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  num_high_confidence_transition_edges <- rep(0, length(high_conf_edges))
  for (p in 1:length(high_conf_edges)){
    num_high_confidence_transition_edges[p] <- sum(genotype_transition[[p]]$transition * high_conf_edges[[p]])
  }

  names(num_high_confidence_transition_edges) <- geno_names
  return(num_high_confidence_transition_edges)
} # end report_num_high_confidence_trans_edge

gene_test_from_snps <- function(){
# goal of this is to rebuild the gene test from snps and keep snp annotation
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
} # end gene_test_from_snps

build_gene_anc_recon_and_conf_from_snp <- function(tr, geno, g_reconstruction_and_confidence, gene_to_snp_lookup_table){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  # VALIDATE INPUTS -----------------------------------------------------------#
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  check_dimensions(geno, Ntip(tr), 2, NULL, 2)
  if (length(g_reconstruction_and_confidence) != ncol(geno)){
    stop("wrong input")
  }
  check_dimensions(gene_to_snp_lookup_table, NULL, 1, 2, 2)
  check_if_binary_vector_numeric(g_reconstruction_and_confidence[[1]]$tip_and_node_recon)
  check_if_binary_vector_numeric(g_reconstruction_and_confidence[[1]]$tip_and_node_rec_conf)

  # FUNCTION ------------------------------------------------------------------#
  tip_nodes_by_snp_mat_recon <- tip_nodes_by_snp_mat_confi <- matrix(0, nrow = (Nnode(tr) + Ntip(tr)), ncol = ncol(geno))
  if (nrow(tip_nodes_by_snp_mat_recon) != length(g_reconstruction_and_confidence[[1]]$tip_and_node_recon)){
    stop("mismatch in size")
  }
  for (k in 1:ncol(geno)){
    tip_nodes_by_snp_mat_recon[ , k] <- g_reconstruction_and_confidence[[k]]$tip_and_node_recon
    tip_nodes_by_snp_mat_confi[ , k] <- g_reconstruction_and_confidence[[k]]$tip_and_node_rec_conf
  }

  row.names(tip_nodes_by_snp_mat_recon) <- row.names(tip_nodes_by_snp_mat_confi) <- c(1:nrow(tip_nodes_by_snp_mat_recon))
  colnames(tip_nodes_by_snp_mat_recon) <- colnames(tip_nodes_by_snp_mat_confi) <- colnames(geno)

  if (nrow(gene_to_snp_lookup_table) != ncol(tip_nodes_by_snp_mat_recon)){
    stop("mismatch")
  }

  recon_times_confi <- tip_nodes_by_snp_mat_recon * tip_nodes_by_snp_mat_confi
  tip_nodes_by_snp_mat_recon_with_gene_id <- rbind(tip_nodes_by_snp_mat_recon, unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))
  tip_nodes_by_snp_mat_confi_with_gene_id <- rbind(tip_nodes_by_snp_mat_confi, unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))
  recon_times_confi_with_gene_id          <- rbind(recon_times_confi,          unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))

  if (nrow(tip_nodes_by_snp_mat_recon_with_gene_id) != (nrow(tip_nodes_by_snp_mat_recon) + 1)){
    stop("rbind didn't work")
  }
  unique_genes <- unique(gene_to_snp_lookup_table[ , 2])

  gene_mat_built_from_snps <- gene_presence_confidence <- gene_absence_confidence <- matrix(0, nrow = nrow(tip_nodes_by_snp_mat_recon), ncol = length(unique_genes))
  for (j in 1:length(unique_genes)){

    # Matrix of just the SNPs found in gene "j"
    temp_mat               <- tip_nodes_by_snp_mat_recon_with_gene_id[1:(nrow(tip_nodes_by_snp_mat_recon_with_gene_id) - 1), tip_nodes_by_snp_mat_recon_with_gene_id[nrow(tip_nodes_by_snp_mat_recon_with_gene_id), ] == unique_genes[j], drop = FALSE]
    temp_recon_times_confi <-          recon_times_confi_with_gene_id[1:(nrow(recon_times_confi_with_gene_id)          - 1),           recon_times_confi_with_gene_id[nrow(recon_times_confi_with_gene_id)        , ] == unique_genes[j], drop = FALSE]
    temp_conf              <- tip_nodes_by_snp_mat_confi_with_gene_id[1:(nrow(tip_nodes_by_snp_mat_confi_with_gene_id) - 1), tip_nodes_by_snp_mat_confi_with_gene_id[nrow(tip_nodes_by_snp_mat_confi_with_gene_id), ] == unique_genes[j], drop = FALSE]
    class(temp_mat) <- class(temp_recon_times_confi) <- class(temp_conf) <- "numeric"
    for (r in 1:nrow(gene_presence_confidence)){
      if (rowSums(temp_mat)[r] == 0 & rowSums(temp_conf)[r] > 0){
        gene_absence_confidence[r, j] <- 1
      }
    }
    gene_mat_built_from_snps[ , j] <- rowSums(temp_mat)
    gene_presence_confidence[ , j] <- rowSums(temp_recon_times_confi)
  }

  gene_mat_built_from_snps <- gene_mat_built_from_snps > 0
  gene_presence_confidence <- gene_presence_confidence > 0
  class(gene_mat_built_from_snps) <- class(gene_presence_confidence) <- "numeric"

  gene_all_confidence <- gene_presence_confidence + gene_absence_confidence

  colnames(gene_mat_built_from_snps)  <- colnames(gene_all_confidence)  <- unique_genes
  row.names(gene_mat_built_from_snps) <- row.names(gene_all_confidence) <- c(1:nrow(gene_mat_built_from_snps))

  gene_list_built_from_snps <- gene_conf_list_built_from_snps<- rep(list(0), length(unique_genes))
  for (m in 1:length(unique_genes)){
    gene_list_built_from_snps[[m]]      <- gene_mat_built_from_snps[ , m, drop = TRUE]
    gene_conf_list_built_from_snps[[m]] <- gene_all_confidence[      , m, drop = TRUE]
  }
  names(gene_list_built_from_snps) <- names(gene_conf_list_built_from_snps) <- unique_genes
  return(list("tip_node_recon" = gene_list_built_from_snps, "tip_node_conf" = gene_conf_list_built_from_snps ))
} # end build_gene_anc_recon_from_snp()


build_node_anc_recon_from_gene_list <- function(gene_list, tr){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  # TODO test this function works.
  gene_list_built_from_snps_just_node_anc_rec <- rep(list(0), length(gene_list$tip_and_node_recon))
  for (m in 1:length(gene_list$tip_and_node_recon)){
    gene_list_built_from_snps_just_node_anc_rec[[m]] <- gene_list[[m]]$tip_and_node_recon[(Ntip(tr) + 1):(Ntip(tr) + Nedge(tr))]
  }
  return(gene_list_built_from_snps_just_node_anc_rec)
} # end build_node_anc_recon_from_gene_list()


build_gene_trans_from_snp_trans <- function(tr, geno, geno_transition, gene_to_snp_lookup_table){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------

  # VALIDATE INPUTS -----------------------------------------------------------#
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  check_dimensions(geno, Ntip(tr), 2, NULL, 2)
  if (length(geno_transition) != ncol(geno)){
    stop("wrong input")
  }
  check_dimensions(gene_to_snp_lookup_table, NULL, 1, 2, 2)
  check_if_binary_vector_numeric(geno_transition[[1]]$transition)

  # FUNCTION ------------------------------------------------------------------#

  edges_by_snp_mat <- matrix(0, nrow = Nedge(tr), ncol = ncol(geno))
  if (nrow(edges_by_snp_mat) != length(geno_transition[[1]]$transition)){
    stop("mismatch in size")
  }
  for (k in 1:ncol(geno)){
    edges_by_snp_mat[ , k] <- geno_transition[[k]]$transition
  }
  row.names(edges_by_snp_mat) <- c(1:nrow(edges_by_snp_mat))
  colnames(edges_by_snp_mat) <- colnames(geno)

  if (nrow(gene_to_snp_lookup_table) != ncol(edges_by_snp_mat)){
    stop("mismatch")
  }

  edges_by_snp_mat_with_gene_id <- rbind(edges_by_snp_mat, unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))
  if (nrow(edges_by_snp_mat_with_gene_id) != (nrow(edges_by_snp_mat) + 1)){
    stop("rbind didn't work")
  }

  unique_genes <- unique(gene_to_snp_lookup_table[ , 2])
  gene_mat_built_from_snps <- matrix(0, nrow = nrow(edges_by_snp_mat), ncol = length(unique_genes))
  for (j in 1:length(unique_genes)){
    temp_mat <- edges_by_snp_mat_with_gene_id[1:(nrow(edges_by_snp_mat_with_gene_id) - 1) , edges_by_snp_mat_with_gene_id[nrow(edges_by_snp_mat_with_gene_id), ] == unique_genes[j], drop = FALSE]
    class(temp_mat) <- "numeric"
    temp_column <- rowSums(temp_mat)
    gene_mat_built_from_snps[ , j] <- temp_column
  }

  gene_mat_built_from_snps <- gene_mat_built_from_snps > 0
  class(gene_mat_built_from_snps) <- "numeric"

  colnames(gene_mat_built_from_snps) <- unique_genes
  row.names(gene_mat_built_from_snps) <- c(1:nrow(gene_mat_built_from_snps))

  gene_list_built_from_snps <- rep(list(0), length(unique_genes))
  for (m in 1:length(unique_genes)){
    gene_list_built_from_snps[[m]] <- gene_mat_built_from_snps[ , m, drop = TRUE]
  }
  names(gene_list_built_from_snps) <- unique_genes

  return(gene_list_built_from_snps)
} # end build_gene_trans_from_snp_trans()


build_gene_genotype_from_snps <- function(geno, gene_to_snp_lookup_table){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  unique_genes <- unique(gene_to_snp_lookup_table[ , 2])
  samples_by_genes <- matrix(0, nrow = nrow(geno), ncol = length(unique_genes))
  colnames(samples_by_genes) <- unique_genes
  row.names(samples_by_genes) <- row.names(geno)


  snp_geno_with_gene_id <- rbind(geno, unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))
  if (nrow(snp_geno_with_gene_id) != (nrow(geno) + 1)){
    stop("rbind didn't work")
  }

  for (j in 1:length(unique_genes)){
    temp_mat <- snp_geno_with_gene_id[1:(nrow(snp_geno_with_gene_id) - 1) , snp_geno_with_gene_id[nrow(snp_geno_with_gene_id), ] == unique_genes[j], drop = FALSE]
    class(temp_mat) <- "numeric"
    temp_column <- rowSums(temp_mat)
    samples_by_genes[ , j] <- temp_column
  }

  samples_by_genes <- samples_by_genes > 0
  class(samples_by_genes) <- "numeric"
  return(samples_by_genes)
}

create_contingency_table <- function(genotype_by_edges, phenotype_by_edges, geno){
  # TODO
  # add names to contingency table results
  # add validations
  # add unit tests
  # add function description
  # genotype_by_edges should be a list
  # phenotype_by_edges should be a numeric vector
  # length of genotype_by_edges should be ncol of geno
  # length of _phenotype_by_edges should be same length as every entry in genotype_by_edges
  # all genotype_by_edges and phenotype_by_edges should be 0s or 1
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------

  all_tables <- rep(list(NULL), length(genotype_by_edges))
  contingency_table <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
  row.names(contingency_table) <- c("geno_present", "geno_absent")
  colnames(contingency_table) <- c("pheno_present", "pheno_absent")
  for (i in 1:length(genotype_by_edges)){
    temp_table <- contingency_table
    temp_table[1, 1] <- sum(genotype_by_edges[[i]]  + phenotype_by_edges == 2,  na.rm = TRUE)
    temp_table[2, 2] <- sum(genotype_by_edges[[i]]  + phenotype_by_edges == 0,  na.rm = TRUE)
    temp_table[1, 2] <- sum(-genotype_by_edges[[i]] + phenotype_by_edges == -1, na.rm = TRUE)
    temp_table[2, 1] <- sum(-genotype_by_edges[[i]] + phenotype_by_edges == 1,  na.rm = TRUE)
    all_tables[[i]] <- temp_table
  }

  names(all_tables) <- colnames(geno)
  return(all_tables)
} # end create_contigency_table()




find_parent_edge <- function(tr, edge_num){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  # given an edge number, get edge number of parent edge
  parent_node <- tr$edge[edge_num, 1]
  parent_edge <- which(tr$edge[ , 2] == parent_node)

  return(parent_edge)
  # TODO it breaks on 1st edge, need to deal with that, but not sure what's best yet
} # end find_parent_edge

check_if_phenotype_normal <- function(pheno, continuous_or_discrete){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  if (continuous_or_discrete == "continuous"){
    result <- shapiro.test(unlist(pheno))
    alpha <- 0.05
    if (result$p < alpha){
      print("Consider normalizing your phenotype")
    }
  }
} # end check_if_phenotype_normal

check_if_convergence_occurs <- function(pheno, tree, continuous_or_discrete){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  if (continuous_or_discrete == "continuous"){
    set.seed(1)
    geiger_BM <- fitContinuous(tree, pheno, model = "BM")
    geiger_white <- fitContinuous(tree, pheno, model = "white")

    if (geiger_white$opt$aicc < geiger_BM$opt$aicc){
      print("WN better than BM")
    }

    # TODO Add this as a plot to output?
    #pdf(paste(test_dir, "/evolutionary_model_barplot.pdf", sep =""))
    #par(mfrow = c(1,1))
    #barplot(c(geiger_BM$opt$aicc, geiger_OU$opt$aicc, geiger_white$opt$aicc),
    #        names.arg = c("Brownian Motion", "Ornstein-Uhlenbeck", "White Noise"),
    #        col = c("black", "black", "black"),ylim = c(0, 600),
    #        ylab = "AICc", xlab = "Evolutionary model",
    #        main = "Model values for toxin")
    #dev.off()
  }
} # end check_if_convergence_occurs()


assign_high_confidence_to_transition_edges <- function(tr, all_confidence_by_edge, genotype_transition_by_edges, geno){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  # VALIDATION
  if (length(genotype_transition_by_edges[[1]]$transition) != Nedge(tr)){
    stop("Dimension mismatch")
  }
  check_for_root_and_bootstrap(tr)
  if (length(all_confidence_by_edge[[1]]) != Nedge(tr)){
    stop("Dimension mismatch")
  }

  # FUNCTION ----------------------------------------------------------------#

  # Identify all edges for which the edge and the parent edge are both high confidence
  edge_and_parent_both_confident <- edge_and_parent_confident_and_trans_edge <- rep(list(rep(0, Nedge(tr))), ncol(geno))
  for (ge in 1:ncol(geno)){
    for (ed in 2:(Nedge(tr) - 1)){ # start at 2 because there isn't a parent edge to edge 1, end at Nedge- 1 because no parent to the last edge either
      parent_edge <- find_parent_edge(tr, ed)
      if (all_confidence_by_edge[[ge]][ed] == 1 & all_confidence_by_edge[[ge]][parent_edge] == 1){
        edge_and_parent_both_confident[[ge]][ed] <- 1
      }
    }
    edge_and_parent_both_confident[[ge]][1] <- all_confidence_by_edge[[ge]][1] # have to accoutn for the fact that there isn't a parent edge to edge 1
    edge_and_parent_both_confident[[ge]][Nedge(tr)] <- all_confidence_by_edge[[ge]][Nedge(tr)] # have to accoutn for the fact that there isn't a parent edge to last edge
  }


  # Identify high confidence transition edges by overlapping the above and transitions
  for (k in 1:ncol(geno)){
    edge_and_parent_confident_and_trans_edge[[k]] <- as.numeric((edge_and_parent_both_confident[[k]] + genotype_transition_by_edges[[k]]$transition) == 2)
  }

  # Return that overlap as high confidence transitions
  return(edge_and_parent_confident_and_trans_edge)
} # end assign_high_confidence_to_transition_edges()

# END OF SCRIPT ---------------------------------------------------------------#
