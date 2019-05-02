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

# LIBRARIES -------------------------------------------------------------------#
# library(ape)       # ape::ace function (ancestral reconstruction)
# library(phytools)  # phylogenetic tree function library
# library(ComplexHeatmap) # to make final plots for discrete phenotypes
# library(phangorn)
# library(pheatmap) # plots for continuous phenotypes
# library(grid) # plots for continuous phenotypes
# library(gridExtra) # plots for continuous phenotypes
# library(ggplot2) # plots for continuous phenotypes

# library(truncnorm) # truncnorm::rtruncnorm function (creates truncated normal distribution) # 2018-08-29 started getting a bug about truncnorm loading....will haveto figure that out later.

# FUNCTIONS FOR PHYC ----------------------------------------------------------#
ancestral_reconstruction_by_ML <- function(tr, mat, num, disc_cont){
  # TODO:
  # 1) Break this function into two subfuctions: one for continuous, one for discrete
  # 2) For those two subfunctions:
  #    a) subfunctionalize
  #    b) comment
  #    c) add output descriptions
  #    d) finish checking all inputs are correctly formatted.

  # Function description -------------------------------------------------------
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # tree:      Phylo.     Rooted phylogenetic tree.
  # mat:       Matrix.    Either the phenotype matrix or the genotype matrix.
  # num:       Numeric.   Indicating the row of the matrix from which the ancestral reconstruction is to be built.
  # disc_cont: Character. Either "discrete" or "continuous".
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec
  # "anc_rec_confidence" = tip_and_node_anc_rec_confidence
  #
  # Check input ----------------------------------------------------------------
  check_is_string(disc_cont)
  check_for_root_and_boostrap(tr)

  # Function -------------------------------------------------------------------

  # Compute ancestral reconstruction
  recon_method <- "ML" # Method: ML.         Maximum Likelihood.
  ML_significance_threshold <- .875 # ML literature suggests that a ratio of 7:1 suggests a high confidence ancestral reconstruction per node .875/.125 = 7

  if (disc_cont == "continuous"){

    # RECONSTRUCTION -----------------------------------------------------------
    set.seed(3)
    reconstruction <- ace(mat[ , num, drop = TRUE], tr, type = disc_cont, method = recon_method, model = "BM", marginal = FALSE)
    ML_anc_rec <- reconstruction$ace # vector containing reconstruction data for all internal nodes (#51 -99) where tips are 1-50.
    tip_and_node_recon <- c(mat[ , num, drop = TRUE], ML_anc_rec)
    names(tip_and_node_recon) <- c(1:sum(Ntip(tr), Nnode(tr)))

    # CONFIDENCE IN RECONSTRUCTION ---------------------------------------------
    tip_and_node_anc_rec_confidence <- rep(1, length(tip_and_node_recon)) # Ancestral reconstruction for continuous values only gives back a 95% CI. We can't use any of this information to decide which nodes are low confidence so treat all reconstructed values as high confidence.

  } else if (disc_cont == "discrete"){

    # RECONSTRUCTION -----------------------------------------------------------
    set.seed(4)
    recon_model <- pick_recon_model(mat, tr, disc_cont, num, recon_method)
    reconstruction <- ace(mat[ , num, drop = TRUE], tr,
                 model = recon_model,
                 type = disc_cont,
                 method = recon_method,
                 marginal = FALSE) # Marginal = FALSE means that the marginal is in fact calculated. Do not set marginal to TRUE, because this only does the condition (only downstream edges considered). ace never calculates the joint.
    ML_anc_rec <- as.numeric(colnames(reconstruction$lik.anc)[apply(reconstruction$lik.anc, 1, which.max)]) # Extract the mostly likely character state using which.max
    names(ML_anc_rec) <- c((Ntip(tr) + 1):(Ntip(tr) + Nnode(tr)))
    tip_and_node_recon <- c(mat[ , num, drop = TRUE], ML_anc_rec)
    names(tip_and_node_recon) <- c(1:sum(Ntip(tr), Nnode(tr)))

    # CONFIDENCE IN RECONSTRUCTION ---------------------------------------------
    anc_rec_confidence <- apply(reconstruction$lik.anc, 1, max) # Get the highest confidence value at each node
    tip_and_node_anc_rec_confidence <- c(rep(1, Ntip(tr)), anc_rec_confidence) # count all tips as high confidence
    tip_and_node_anc_rec_confidence <- discretize_confidence_using_threshold(tip_and_node_anc_rec_confidence, ML_significance_threshold) # Count anything lower than threshold as low confidence

  } else {
    stop("Ancestral reconstruction cannot run.")
  }

  # AS A MATRIX, IN THE SAME FORMAT AS TREE$EDGE. SO NODE VALUES WILL APPEAR MULTIPLE TIMES IN THE tr.
  # THIS FORMAT WILL MAKE IT MUCH EASIER TO CALCULATE PHENOTYPE CHANGE ON EDGES.
  reconstruction_as_edge_mat <- tr$edge
  for (k in 1:nrow(reconstruction_as_edge_mat)){
    reconstruction_as_edge_mat[k, 1] <- tip_and_node_recon[tr$edge[k, 1]]
    reconstruction_as_edge_mat[k, 2] <- tip_and_node_recon[tr$edge[k, 2]]
  }


  if (Nnode(tr) != length(ML_anc_rec)){
    stop("Ancestral reconstruction of nodes is the wrong length.")
  }

  results <- list("node_anc_rec" = ML_anc_rec,
                  "tip_and_node_rec_conf" = tip_and_node_anc_rec_confidence,
                  "recon_edge_mat" = reconstruction_as_edge_mat,
                  "tip_and_node_recon" = tip_and_node_recon
                  )
  return(results)
} # end ancestral_reconstruction_by_ML()

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

is_tip <- function(node_num, tr){
  # Function description ------------------------------------------------------#
  # Test if a node is a tip or an internal node.
  #
  # Inputs:
  # node_num: Integer. Index of node.
  # tr: phylogenetic tree.
  #
  # Output:
  # Logical. TRUE OR FALSE.
  #
  # Check input ---------------------------------------------------------------#
  check_node_is_in_tree(node_num, tr)
  #
  # Function & return output --------------------------------------------------#
  return(node_num <= Ntip(tr))
} # end is_tip()

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
  # TODO add unit tests
  # TODO deal with cases when there isn't enough data (aren't at least 1 of least t_index and non_t_index)
  # t_index = transition edges
  # not_t_index = non transition edges
  # phenotype_by_edges = change in phenotype on each edge

  # Check input ---------------------------------------------------------------#
  if (length(t_index) < 1 | length(non_t_index) < 1){
    stop("Not enough high confidence transition edges to use for KS test.")
  }


  # subset index to only high confidence edges:
  p_trans_delta     <- calculate_phenotype_change_on_edge(t_index,     phenotype_by_edges)
  p_non_trans_delta <- calculate_phenotype_change_on_edge(non_t_index, phenotype_by_edges)

  #print("ks test")
  #print(length(t_index))
  #print(length(non_t_index))
  #print(p_trans_delta)
  #print(p_non_trans_delta)

  ks_results        <- ks.test(p_trans_delta, p_non_trans_delta)
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
  check_for_root_and_boostrap(tr)
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

  # Check and return output ----------------------------------------------------
  results <- list("pvals" = pvals, "ks_statistics" = empirical_ks_stat_list,
                  "observed_pheno_trans_delta" = observed_pheno_trans_delta,
                  "observed_pheno_non_trans_delta" = observed_pheno_non_trans_delta,
                  "trans_median" = trans_median, "all_edges_median" = all_edges_median,
                  "num_genotypes" = ncol(mat), "observed_ks_stat" = observed_ks_stat) # 2018-11-28
  return(results)
} # end calculate_genotype_significance()


collapse_identical_genotypes <- function(genotype_matrix, output_directory, output_name){
  # TO DO:
  # Add a check of the output format.
  #
  # Function description -------------------------------------------------------
  # Save any genotype with the same absence/presence pattern into an Rdata object
  # and remove those duplicates from the genotype matrix to reduce the penalization
  # induced by multiple test correction.
  #
  # Input:
  # genotype_matrix.  Matrix. Genotypes in the columns, samples in the rows. All entries are either 0 or 1. No NAs.
  # output_directory. Character. Path to output directory.
  # output_name.      Character. Identifier for output files.
  #
  # Output:
  # None.
  #
  # Check input ----------------------------------------------------------------
  check_if_binary_matrix(genotype_matrix)
  check_is_string(output_name)
  check_if_dir_exists(output_directory)

  # Function -------------------------------------------------------------------
  genotype_matrix <- t(genotype_matrix)
  genotype_matrix  <- as.data.frame(genotype_matrix)
  duplicate_factor <- interaction(genotype_matrix, drop = TRUE)  #interaction function requires data.frame
  genotype_groups  <- rep(list(NULL), length(unique(duplicate_factor)))
  non_duplicates   <- NULL
  for (i in 1:length(unique(duplicate_factor))){
    genotype_groups[i]  <- list(row.names(genotype_matrix)[duplicate_factor == unique(duplicate_factor)[i]])
    non_duplicates <- append(non_duplicates, genotype_groups[[i]][1])
  }
  genotype_matrix <- genotype_matrix[row.names(genotype_matrix) %in% non_duplicates, , drop = FALSE]
  genotype_matrix <- as.matrix(genotype_matrix)
  genotype_matrix <- t(genotype_matrix)
  file_name <- create_file_name(output_directory, output_name, "genotypes_with_identical_patterns.rda")

  # Check and return output ----------------------------------------------------
  save(genotype_groups, file = file_name)
  return(genotype_matrix)
} # end collapse_identical_genotypes()

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
  #file_name <-  paste(output_dir, "/", format(Sys.time(), "%Y-%m-%d_"), "phyc_", output_name, "_", other_info, sep = "")
  file_name <-  paste(output_dir, "/", "phyc_", output_name, "_", other_info, sep = "")
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
  set.seed(50)
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
  check_for_root_and_boostrap(tree)
  check_dimensions(phenotype_matrix, Ntip(tree), 2, 1, 1)
  check_dimensions(genotype_matrix, Ntip(tree), 2, NULL, 1)
  check_if_binary_matrix(genotype_matrix)
  check_rownames(phenotype_matrix, tree)
  check_rownames(genotype_matrix, tree)

  results <- list("tree" = tree, "phenotype" = phenotype_matrix, "genotype" = genotype_matrix)
  return(results)
} # end create_test_data()

discretize_confidence_using_threshold <- function(confidence_vector, threshold){
  confidence_vector[confidence_vector  < threshold] <- 0
  confidence_vector[confidence_vector >= threshold] <- 1
  return(confidence_vector)
} # end discretize_confidence_using_threshold()

identify_short_edges <- function(tr){
  # Removes any edges that make up a 10% or more of the total tr length.

  short_edges <- rep(1, Nedge(tr))
  while(max(tr$edge.length[as.logical(short_edges)]) >= (0.1 * sum(tr$edge.length[as.logical(short_edges)]))){
    short_edges[tr$edge.length == max(tr$edge.length[as.logical(short_edges)])] <- 0
  }
  return(short_edges)
} # end identify_short_edges()

reorder_tips_and_nodes_to_edges <- function(tips_and_node_vector, tr){
  ordered_by_edges <- rep(NA, Nedge(tr))
  for (i in 1:Nedge(tr)){
    ordered_by_edges[i] <- tips_and_node_vector[tr$edge[i, 2]]
  }
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
  results                        <- list("hit_pvals" = fdr_corrected_pvals, "sig_pvals" = sig_pvals)

  # Check and return output ----------------------------------------------------
  return(results)
} # end get_sig_hits_while_correcting_for_multiple_testing

keep_at_least_two_high_conf_trans_edges <- function(genotype_transition, genotype_confidence){
  # TO DO:
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
  # TO DO:
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

plot_continuous_phenotype <- function(tr, pheno_vector, pheno_anc_rec){
  plot_p_recon <- contMap(tr, pheno_vector, method = "user", anc.states = pheno_anc_rec, plot = FALSE)
  plot(plot_p_recon,
       add = TRUE,
       ylim = c(-1 / 25 * Ntip(tr), Ntip(tr)),
       colors = plot_p_recon$cols,
       lwd = 4,
       ftype = "off",
       offset = 1.7)
} # end plot_continuous_phenotype()

plot_significant_hits <- function(disc_cont, tr, a, dir, name, pval_all_transition, pheno_vector, annot, perm, results_all_trans, pheno_anc_rec, geno_reconstruction, geno_confidence, geno_transition, geno, pheno_recon_ordered_by_edges, tr_and_pheno_hi_conf, all_trans_sig_hits){
  # Function description -------------------------------------------------------
  #
  # Inputs:
  # tr.                  Phylo.
  # a.                   Number. Alpha.
  # dir.                 Character. Output path.
  # name.                Character. Output name.
  # pheno_vector.        Vector.
  # annot.               Matrix.
  # perm.                Number.

  # Output:
  # None.
  #
  # Check inputs ---------------------------------------------------------------
  check_for_root_and_boostrap(tr)
  check_if_alpha_valid(a)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_if_vector(pheno_vector)
  check_if_permutation_num_valid(perm)

  # Function -------------------------------------------------------------------
  trans_edge_mat <- NULL
  for (i in 1:length(geno_transition)){
    trans_edge_mat <- cbind(trans_edge_mat, geno_transition[[i]]$transition)
  }
  colnames(trans_edge_mat) <- colnames(geno)

  # TODO Update trans_edge_mat to exclude low confidence  edges, currently it includes transition edges (all high and some low confidence transitions)
  for (c in 1:ncol(trans_edge_mat)){
    trans_edge_mat[(1:Nedge(tr))[geno_confidence[[c]] == 0], c] <- NA
  }
  # end update trans_edge_mat
  ph_trans <- abs(pheno_recon_ordered_by_edges[ , 1] - pheno_recon_ordered_by_edges[ , 2])

  print("0")
  p_trans_mat <- matrix(ph_trans, nrow = length(ph_trans), ncol = 1)
  colnames(p_trans_mat) <- "delta_pheno"
  p_trans_mat <- as.data.frame(round(p_trans_mat, 2))

  significant_loci <- data.frame("locus" = rep("not_sig", ncol(trans_edge_mat)), stringsAsFactors = FALSE)
  row.names(significant_loci) <- colnames(trans_edge_mat)
  significant_loci[row.names(significant_loci) %in% row.names(all_trans_sig_hits), ] <- "sig"

  log_p_value <- data.frame(-log(pval_all_transition$hit_pvals))
  column_annot <- cbind(significant_loci, log_p_value)
  print("1")

  row.names(p_trans_mat) <- row.names(trans_edge_mat) <- c(1:Nedge(tr))
  # end heatmap prep
  ann_colors = list(
    locus = c(not_sig = "white", sig = "blue")
  )

  sorted_trans_edge_mat <-           trans_edge_mat[match(row.names(p_trans_mat)[order(p_trans_mat[ , 1])], row.names(trans_edge_mat)), ]
  ordered_by_p_val      <- sorted_trans_edge_mat[ , match(row.names(log_p_value)[order(log_p_value[ , 1])], colnames(sorted_trans_edge_mat))]
  column_annot_ordered_by_p_val <-     column_annot[match(row.names(log_p_value)[order(log_p_value[ , 1])], row.names(column_annot)), ]
  colnames(column_annot) <- colnames(column_annot_ordered_by_p_val) <- c("locus", "-ln(p-val)")

  fname <- create_file_name(dir, name, paste("summary_and_sig_hit_results.pdf", sep = "")) # 2019-02-25 removed counter from pdf file name because combining.
  pdf(fname, width = 16, height = 20)

  make_manhattan_plot(dir, name, pval_all_transition$hit_pvals, a, "transition")
  print("2")
  cell_width_value <- 1.5
  if (ncol(ordered_by_p_val) < 50){
    cell_width_value <- 10
  }


  pheatmap( # Plot the heatmap
    ordered_by_p_val,
    main          = paste0("Edges:\n hi conf trans vs delta pheno"),
    cluster_cols  = FALSE,
    na_col = "grey",
    cluster_rows  = FALSE,
    show_rownames = FALSE,
    color = c("white", "black"),
    annotation_col = column_annot_ordered_by_p_val,
    annotation_row = p_trans_mat,
    annotation_colors = ann_colors,
    show_colnames = TRUE,
    cellwidth = cell_width_value)

  pheatmap( # Plot the heatmap
    sorted_trans_edge_mat,
    main          = paste0("Edges:\n hi conf trans vs delta pheno"),
    cluster_cols  = TRUE,
    cluster_rows  = FALSE,
    show_rownames = FALSE,
    color = c("white", "black"),
    annotation_col = column_annot,
    annotation_row = p_trans_mat,
    annotation_colors = ann_colors,
    show_colnames = TRUE,
    cellwidth = cell_width_value,
    na_col = "grey")
  print("3")
  # ONLY MAKE THE FOLLOWING PLOTS FOR SIGNIFICANT LOCI
  counter <- 0 # TODO can I get rid of counter now? 2019-02-25
  for (j in 1:nrow(pval_all_transition$hit_pvals)){
    if (pval_all_transition$hit_pvals[j, 1] < a){
      print("making significant plots")
      counter <- counter + 1
      par(mfrow = c(3, 3))
      par(mgp   = c(3, 1, 0))
      par(oma   = c(0, 0, 4, 0))
      par(mar = c(4, 4, 4, 4))

      plot_continuous_phenotype(tr, pheno_vector, pheno_anc_rec)
      plot_tree_with_colored_edges(tr, geno_reconstruction, geno_confidence, "grey", "red", paste0(row.names(pval_all_transition$hit_pvals)[j], "\n Genotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", j)
      plot_tree_with_colored_edges(tr, geno_transition,     geno_confidence, "grey", "red", "Genotype transition edge:\n Red = transition; Black = No transition", annot, "trans", j)
      histogram_all_delta_pheno_overlaid_with_high_conf_delta_pheno(p_trans_mat, geno_confidence, tr, j)
      histogram_abs_high_confidence_delta_pheno_highlight_transition_edges(results_all_trans, tr, j, "grey", "red") # 6. Delta phenotype histogram. Note that results_all_trans$observed_pheno_non_trans_delta[[j]] is already subset down to only high confidence edges
      histogram_raw_high_confidence_delta_pheno_highlight_transition_edges(geno_transition, geno_confidence, pheno_recon_ordered_by_edges, tr, j, "grey", "red")

      hist(log(results_all_trans$ks_statistics[[j]]),
           breaks = perm/10, col = "grey", border = FALSE,
           main = paste("Null distribution of KS statistic for all transitions.\n Red = Observed KS statistic.\n p-value = ", round(pval_all_transition$hit_pvals[j, 1], 10), "\np-value rank = ", rank(pval_all_transition$hit_pvals)[j], sep = ""),
           ylab = "Count",
           xlab = "ln(KS statistic)",
           xlim = c(min(log(as.numeric(results_all_trans$observed_ks_stat[j])), log(results_all_trans$ks_statistics[[j]])), 0))
      abline(v = log(as.numeric(results_all_trans$observed_ks_stat[j])), col = "red")

      pheatmap(
        sorted_trans_edge_mat[ , j, drop = FALSE],
        main          = paste0(row.names(pval_all_transition$hit_pvals)[j], "\n Tree edges: hi conf trans vs delta pheno"),
        cluster_cols  = FALSE,
        cluster_rows  = FALSE,
        na_col = "grey",
        show_rownames = FALSE,
        color = c("white", "black"),
        annotation_row = p_trans_mat,
        show_colnames = TRUE,
        cellwidth = 20)
    }
  }

  dev.off()

  print("4")
  trans_dir_edge_mat <- NULL
  for (i in 1:length(geno_transition)){
    trans_dir_edge_mat <- cbind(trans_dir_edge_mat, geno_transition[[i]]$trans_dir)
  }
  print("dim trans_dir_edge_mat")
  print(dim(trans_dir_edge_mat))
  print("geno_transition length")
  print(length(geno_transition))
  print("[[i]]$trans_dir")
  print(geno_transition[[i]]$trans_dir)


  print("5")
  colnames(trans_dir_edge_mat) <- colnames(geno)
  print("6")
  # TODO Update trans_edge_mat to exclude low confidence  edges, currently it includes transition edges (all high and some low confidence transitions)
  for (c in 1:ncol(trans_dir_edge_mat)){
    trans_dir_edge_mat[(1:Nedge(tr))[geno_confidence[[c]] == 0], c] <- NA
  }
  print("7")

  all_tables <- all_lists <- rep(list(NULL), ncol(trans_dir_edge_mat))
  print("8")
  delta_pheno_table <- matrix(0, nrow = 3, ncol = 1)
  row.names(delta_pheno_table) <- c("geno_parent_0_child_1", "geno_parent_1_child_0", "geno_no_change")
  colnames(delta_pheno_table) <- c("sum(|delta_phenotype|)")
  print("9")
  for (i in 1:ncol(trans_dir_edge_mat)){
    temp_table <- delta_pheno_table
    temp_table[1, 1] <- sum(p_trans_mat[which(trans_dir_edge_mat[ , i] == 1),  1], na.rm = TRUE)
    temp_table[2, 1] <- sum(p_trans_mat[which(trans_dir_edge_mat[ , i] == -1), 1], na.rm = TRUE)
    temp_table[3, 1] <- sum(p_trans_mat[which(trans_dir_edge_mat[ , i] == 0),  1], na.rm = TRUE)
    all_tables[[i]] <- temp_table
  }
  print("10")
  names(all_tables) <- colnames(geno)
  delta_pheno_table <- all_tables
  print("11")
  delta_pheno_list <- rep(list(0), 3)
  names(delta_pheno_list) <- c("geno_parent_0_child_1", "geno_parent_1_child_0", "geno_no_change")
  for (i in 1:ncol(trans_dir_edge_mat)){
    temp_list <- delta_pheno_list
    temp_list[[1]] <- p_trans_mat[which(trans_dir_edge_mat[ , i] == 1),  1]
    temp_list[[2]] <- p_trans_mat[which(trans_dir_edge_mat[ , i] == -1), 1]
    temp_list[[3]] <- p_trans_mat[which(trans_dir_edge_mat[ , i] == 0),  1]
    all_lists[[i]] <- temp_list
    names(all_lists)[i] <- colnames(trans_dir_edge_mat)[i]
  }
  print("12")
  delta_pheno_list <- all_lists
  print("done with this function")
  results <- list("trans_dir_edge_mat" = trans_dir_edge_mat, "p_trans_mat" = p_trans_mat, "delta_pheno_table" = delta_pheno_table, "delta_pheno_list" = delta_pheno_list )
  return(results)
} # end plot_significant_hits()

calc_raw_diff <- function(edge_list, ph_edges){
  delta <- rep(NA, length(unique(edge_list)))
  for (j in 1:length(edge_list)){
    delta[j] <- ph_edges[edge_list[j], 1] - ph_edges[edge_list[j], 2]
  }
  return(delta)
} # end calc_raw_diff()

histogram_raw_high_confidence_delta_pheno_highlight_transition_edges <- function(geno_transition, geno_confidence, pheno_recon_ordered_by_edges, tr, index, non_trans_color, trans_color){
  trans_index     <- c(1:Nedge(tr))[as.logical(geno_transition[[index]]$transition)]
  non_trans_index <- c(1:Nedge(tr))[!geno_transition[[index]]$transition]
  hi_conf_trans_index     <- trans_index[    as.logical(geno_confidence[[index]][trans_index    ])]
  hi_conf_non_trans_index <- non_trans_index[as.logical(geno_confidence[[index]][non_trans_index])]
  #hi_conf_pheno_trans_delta     <- calculate_phenotype_change_on_edge(hi_conf_trans_index,     pheno_recon_ordered_by_edges)
  #hi_conf_pheno_non_trans_delta <- calculate_phenotype_change_on_edge(hi_conf_non_trans_index, pheno_recon_ordered_by_edges)
  raw_trans_delta     <- calc_raw_diff(hi_conf_trans_index,     pheno_recon_ordered_by_edges)
  raw_non_trans_delta <- calc_raw_diff(hi_conf_non_trans_index, pheno_recon_ordered_by_edges)

  hist(raw_non_trans_delta,
       main = paste("Raw delta phenotype on only high confidence edges\n # trans edge= ",
                    length(raw_trans_delta), "\n# non trans edge =",
                    length(raw_non_trans_delta), sep = ""),
       breaks = Nedge(tr)/4,
       col = non_trans_color,
       border = FALSE,
       xlim = c(min(raw_trans_delta, raw_non_trans_delta), max(raw_trans_delta, raw_non_trans_delta)),
       ylab = "Count",
       cex = .8,
       xlab = "Raw delta phenotype")

  hist(raw_trans_delta,
       breaks = Nedge(tr)/4,
       col = trans_color,
       border = trans_color,
       add = TRUE,
       xlim = c(min(raw_trans_delta, raw_non_trans_delta), max(raw_trans_delta, raw_non_trans_delta)))
}

histogram_abs_high_confidence_delta_pheno_highlight_transition_edges <- function(results_all_trans, tr, index, non_trans_color, trans_color){
  par(mar = c(4, 4, 4, 4))
  hist(results_all_trans$observed_pheno_non_trans_delta[[index]],
       breaks = Nedge(tr)/4,
       col = non_trans_color,
       border = FALSE,
       ylab = "Count",
       xlab = "|Delta phenotype|",
       xlim = c(0, max(results_all_trans$observed_pheno_trans_delta[[index]], results_all_trans$observed_pheno_non_trans_delta[[index]])),
       cex = .8,
       main = paste("|delta phenotype| on only high confidence edges \n # non-trans edges = ",
                    length(results_all_trans$observed_pheno_non_trans_delta[[index]]),
                    "\n # trans edges = ", length(results_all_trans$observed_pheno_trans_delta[[index]]), sep = ""))


  hist(results_all_trans$observed_pheno_trans_delta[[index]],
       breaks = Nedge(tr)/4,
       col = trans_color,
       border = trans_color,
       add = TRUE,
       xlim = c(0, max(results_all_trans$observed_pheno_trans_delta[[index]], results_all_trans$observed_pheno_non_trans_delta[[index]])))
} # end histogram_abs_high_confidence_delta_pheno_highlight_transition_edges()

histogram_all_delta_pheno_overlaid_with_high_conf_delta_pheno <- function(p_trans_mat, geno_confidence, tr, index){
  edge_num <- length(unlist(p_trans_mat))
  hi_edge_num <- length(unlist(p_trans_mat)[as.logical(geno_confidence[[index]])])
  title <- paste("|delta phenotype| on all edges \n Light Green: all edges = ", edge_num, "\n Grey: high confidence edges = ", hi_edge_num, sep = "")
  delta_phenotype_on_all_edges <- as.numeric(unlist(p_trans_mat))
  hist(delta_phenotype_on_all_edges,
       breaks = Nedge(tr)/4,
       col = rgb(0, 0.5, 0, 0.25),
       border = FALSE,
       ylab = "Count",
       xlab = "Delta phenotype",
       main = title)

  delta_phenotype_on_high_confidence_edges <- as.numeric(unlist(p_trans_mat))[as.logical(geno_confidence[[index]])]
  hist(delta_phenotype_on_high_confidence_edges, # plot phenotype transition only high confidence edges for this genotype
       breaks = Nedge(tr)/4,
       col = rgb(0, 0, 0, 0.25),
       border = FALSE,
       ylab = "Count",
       add = TRUE)
} # end histogram_all_delta_pheno_overlaid_with_high_conf_delta_pheno()


save_data_table <- function(matrix, output_dir, pheno_geno_name, extension){
  write.table(matrix,
              file = paste0(output_dir, "/phyc_", pheno_geno_name, extension),
              sep = "\t",
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE)
} # end save_data_table()

save_results_as_r_object <- function(dir, name, object){
  save(object, file = paste0(dir, "/phyc_", name, ".rda"))
} # end save_results_as_r_object()

read_in_arguments <- function(args){
  # TODO: Update description & reduce if statements down to functions
  # Function description ------------------------------------------------------#
  # Read in the commandline arguments.
  #
  # Inputs:
  # args. Command line inputs.
  #
  # Output:
  # test
  # tree
  # phenotype
  # genotype
  # output_name
  # output_dir
  # alpha
  # discrete_or_continuous
  # annotation
  #
  # Check inputs, function, & check and return outputs -------------------------
  if (length(args) == 2){
    if (args[1] == "test"){
      test <- TRUE
    } else {
      stop("The first argument for a test run should be \"test.\"")
    }
    test_data              <- create_test_data()
    phenotype              <- test_data$phenotype
    tree                   <- test_data$tree
    genotype               <- test_data$genotype
    output_name            <- "test_data"
    output_dir             <- args[2]
    perm                   <- 1000
    alpha                  <- 0.05
    discrete_or_continuous <- check_input_format(phenotype, tree, genotype, output_name, output_dir, perm, alpha)
    annotation             <- NULL
    results <- list("test" = test, "tree" = tree, "phenotype" = phenotype,
                    "genotype" = genotype, "output_name" = output_name,
                    "output_dir" = output_dir, "perm" = perm, "alpha" = alpha,
                    "discrete_or_continuous" = discrete_or_continuous,
                    "annotation" = annotation)
    return(results)
  } else if (length(args) == 9){
    test                   <- FALSE
    phenotype              <- read_in_tsv_matrix(args[1])
    tree                   <- read.tree(args[2])
    if (!is.rooted(tree)){
      tree <- midpoint(tree)
    }
    genotype               <- read_in_tsv_matrix(args[3])
    output_name            <- args[4] # Ex: log_toxin_snp_stop
    output_dir             <- args[5] # Directory in which all output files will be saved
    perm                   <- as.numeric(args[6]) #typically 10,000
    alpha                  <- as.numeric(args[7])
    bootstrap_cutoff       <- as.numeric(args[8]) # typically 0.70
    annotation             <- read_in_tsv_matrix(args[9])
    discrete_or_continuous <- check_input_format(phenotype, tree, genotype, output_name, output_dir, perm, alpha, annotation)
    results <- list("test" = test, "tree" = tree, "phenotype" = phenotype,
                    "genotype" = genotype, "output_name" = output_name,
                    "output_dir" = output_dir, "perm" = perm, "alpha" = alpha,
                    "discrete_or_continuous" = discrete_or_continuous,
                    "annotation" = annotation, "bootstrap_cutoff" = bootstrap_cutoff)
    return(results)
  } else if (length(args) == 8){
    test                   <- FALSE
    phenotype              <- read_in_tsv_matrix(args[1])
    tree                   <- read.tree(args[2])
    # added is.rooted if statement on 2018-09-25 to deal with odd midpoint rooting issue for PSM dataset.
    if (!is.rooted(tree)){
      tree <- midpoint(tree)
    }
    # End of section added on 2018-09-25
    genotype               <- read_in_tsv_matrix(args[3])
    output_name            <- args[4] # Ex: log_toxin_snp_stop
    output_dir             <- args[5] # Directory in which all output files will be saved
    perm                   <- as.numeric(args[6]) #typically 10,000
    alpha                  <- as.numeric(args[7])
    bootstrap_cutoff       <- as.numeric(args[8])
    discrete_or_continuous <- check_input_format(phenotype, tree, genotype, output_name, output_dir, perm, alpha, NULL)
    results <- list("test" = test, "tree" = tree, "phenotype" = phenotype,
                    "genotype" = genotype, "output_name" = output_name,
                    "output_dir" = output_dir, "perm" = perm, "alpha" = alpha,
                    "discrete_or_continuous" = discrete_or_continuous,
                    "annotation" = NULL, "bootstrap_cutoff" = bootstrap_cutoff)
    return(results)
  } else {
    stop("2, 8 or 9 inputs required. \n
         Either: 1. test 2. output directory or \n
         1. Phenotype 2. Tree 3. Genotype 4. Output name 5. Output directory \n
         6. Number of permutations 7. Alpha 8. Bootstrap confidence threshold and optional 9. Annotation")
  }
} # end read_in_arguments()

read_in_tsv_matrix <- function(mat){
  # Function description -------------------------------------------------------
  # Read in the standardized GWAS matrix format: rows correspond to samples, columns correspond to genotypes/phenotypes.
  # Data are tab-separated.
  #
  # Inputs:
  # mat. Character. Path to matrix file.
  #
  # Output:
  # temp. Matrix.
  #
  # Check inputs ---------------------------------------------------------------
  check_is_string(mat)
  check_file_exists(mat)

  # Function -------------------------------------------------------------------
  temp <- read.table(mat,
                     sep = "\t",
                     row.names = 1,
                     header = TRUE,
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
  temp <- as.matrix(temp)

  # Check and return output ----------------------------------------------------
  if (class(temp) != "matrix"){stop("Ouput is incorrectly formatted")}
  return(temp)
} # end read_in_tsv_matrix()



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
  check_for_root_and_boostrap(tr)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_if_binary_matrix(mat)

  # Function -------------------------------------------------------------------
  geno_to_drop <- rep(FALSE, ncol(mat))
  geno_to_drop[colSums(mat) <= 1] <- TRUE
  geno_to_drop[colSums(mat) == (Ntip(tr) - 1)] <- TRUE
  geno_to_drop[colSums(mat) == Ntip(tr)] <- TRUE


  dropped_genotype_names <- colnames(mat)[geno_to_drop]
  #filename <- paste(dir, "/phyc_", name, "_genotypes_dropped_because_convergence_not_possible.txt", sep = "")
  #writeLines(dropped_genotype_names, con = filename, sep = "\n")

  mat <- mat[ , !geno_to_drop, drop = FALSE]

  # Check and return output ----------------------------------------------------
  check_if_binary_matrix(mat)
  check_rownames(mat, tr)
  results <- list("mat" = mat, "dropped_genotype_names" = dropped_genotype_names)
  #return(mat)
  return(results)
} # end reduce_redundancy()

save_hits <- function(hits, output_dir, output_name, pval_name){
  # Function description -------------------------------------------------------
  # Create a file name and save results to that file name.
  #
  # Inputs:
  # hits.        Vector. Pvals.
  # output_dir.  Character.
  # output_name. Character.
  # pval_name.   Character.
  #
  # Output:
  # None
  #
  # Check inputs ---------------------------------------------------------------
  check_is_string(pval_name)
  check_is_string(output_name)
  check_if_dir_exists(output_dir)

  # Function -------------------------------------------------------------------
  if (nrow(hits) > 0){
    file_name <- create_file_name(output_dir, output_name, pval_name)
    hit_and_rank <- cbind(hits, rank(hits))
    colnames(hit_and_rank)[2] <- "rank"
    write.table(x = hit_and_rank, file = paste(file_name, ".tsv", sep = ""), sep = "\t", quote = FALSE, eol = "\n", row.names = TRUE, col.names = TRUE)
  } else { # added 2018-11-13
    file_name <- create_file_name(output_dir, output_name, pval_name)
    empty <- matrix(NA, nrow = 1, ncol = 1)
    colnames(empty) <- "no_sig_hits"
    row.names(empty) <- "no_sig_hits"
    write.table(x = empty, file = paste(file_name, ".tsv", sep = ""), sep = "\t", quote = FALSE, eol = "\n", row.names = TRUE, col.names = TRUE)
  }
} #end save_hits()

make_manhattan_plot <- function(outdir, geno_pheno_name, pval_hits, alpha, trans_or_recon){
  # Create negative log p-values with arbitrary locus numbers
  neg_log_p_value <- data.frame(-log(pval_hits))
  neg_log_p_with_num <- cbind(1:nrow(neg_log_p_value), neg_log_p_value)
  colnames(neg_log_p_with_num)[1] <- "locus"
  sig_temp <-subset(neg_log_p_with_num, neg_log_p_with_num[ , 2] > -log(alpha))
  ymax <- max(-log(0.01), neg_log_p_with_num[ , 2, drop = TRUE])
  #pdf(paste0(outdir, "/phyc_", trans_or_recon, "_", geno_pheno_name, "_manhattan_plot.pdf"))
  with(neg_log_p_with_num,
       plot(x = neg_log_p_with_num[ , 1],
            y = neg_log_p_with_num[ , 2, drop = TRUE],
            type = "p",
            main = paste(trans_or_recon, "phyC", geno_pheno_name, sep = " "),
            col = rgb(0, 0, 0, 0.3),
            pch = 19,
            xlab = "Genetic locus",
            ylim = c(0, ymax),
            ylab = "-ln(p-val)" ))

  abline(h = -log(alpha), col = "red")
  if (nrow(sig_temp) > 0){
    text(x = sig_temp[ , 1], y = sig_temp[ , 2], labels = row.names(sig_temp), pos = 1, cex = 0.7)
  }
  #dev.off()
} #end make_manhattan_plot()

plot_tree_with_colored_edges <- function(tr, edges_to_highlight, geno_confidence, edge_color_na, edge_color_bright, title, annot, trans_or_recon, index){
  edge_color <- rep("black", Nedge(tr))
  if (trans_or_recon == "recon"){
    edge_color[edges_to_highlight[[index]] == 1] <- edge_color_bright
  } else if (trans_or_recon == "trans"){
    edge_color[edges_to_highlight[[index]]$transition == 1] <- edge_color_bright
  }
  edge_color[geno_confidence[[index]] == 0] <- edge_color_na # grey out long edges and low ML bootstrap support
  par(mar = c(4, 4, 4, 4))
  plot(tr,
       font = 1,
       edge.color = edge_color,
       main = title,
       use.edge.length = FALSE,
       label.offset = 3,
       adj = 0)
  tiplabels(pch = 21,
            col = annot[ , 2],
            adj = 2,
            bg = annot[ , 2],
            cex = 0.75)
  if (!is.null(annot)){
    legend("bottomleft",
           legend = unique(annot[ , 1]),
           col = unique(annot[ , 2]),
           lty = 1,
           ncol = length(unique(annot[ , 1])),
           lwd = 5,
           cex = 0.6)
  }
} # end plot_tree_with_colored_edges()



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


plot_sig_hits_summary <- function(heat_tr, tr, g_mat, p_mat, annot, sig_hits, heatmap_title, filename_start, high_conf_edges, recon_or_trans, short, pheno_conf, bootstrap, pheno_recon_or_trans_by_edge, all_high_conf_edges, geno_trans_or_recon, output_name){
  hits <- g_mat[ , colnames(g_mat) %in%  row.names(sig_hits), drop = FALSE]

  if (recon_or_trans == "trans"){
    filename <- paste(filename_start, "tree.pdf", sep = "")
  } else {
    filename <- paste(filename_start, "tree.pdf", sep = "")
  }
  pdf(filename, width = 16, height = 20)
  par(mfrow = c(1, 2))
  par(mgp = c(3, 1, 0))
  par(oma = c(0, 0, 4, 0))
  # 1. High confidence phenotype (reconstruction or transition)
  edge_color <- pheno_recon_or_trans_by_edge
  edge_color[edge_color == 1] <- "green"
  edge_color[high_conf_edges == 0] <- "grey"
  edge_color[edge_color == 0] <- "black"
  if(recon_or_trans == "trans"){
    title <- paste("Phenotype transitions\n Black=no trans Green=trans Grey=low conf", sep = "")
  } else if (recon_or_trans == "recon"){
    title <- paste("Phenotype reconstruction\n Black=absent Blue=present Grey=low conf", sep = "")
  }
  plot(tr,
       main = title,
       type = "phylogram",
       use.edge.length = TRUE,
       edge.width = 2.0,
       cex = 0.4,
       cex.main = 1.0,
       no.margin = FALSE,
       label.offset = 0.0001,
       edge.color = edge_color)

  # 2. Low confidence edges shown (Low confidence due to either: ML phenotype, bootstrap values, and/or long tree edges)
  edge_color <- rep("black", Nedge(tr))
  edge_color[pheno_conf == 0] <- "blue"
  edge_color[bootstrap == 0] <- "purple"
  edge_color[short == 0] <- "red"
  title <- paste("Low confidence edges\nRed=long Purple=bootstrap Blue=pheno\nBlack=high conf all metrics", sep = "")
  plot(tr,
       main = title,
       type = "phylogram",
       use.edge.length = TRUE,
       edge.width = 2.0,
       cex = 0.4,
       cex.main = 1.0,
       no.margin = FALSE,
       label.offset = 0.0001,
       edge.color = edge_color)


  #mtext(output_name,
  #      outer = TRUE,
  #      side = 3,
  #      cex = 1.2,
  #      line = 1)
  #dev.off()

  # Now 2 PDF per significant hit: 1 is the heatmap and 1 is the geno & pheno trans/recon

  counter <- 0
  for (i in 1:ncol(g_mat)){
    # Subset to only significant hits:
    if (colnames(g_mat)[i] %in% row.names(sig_hits)){
      counter <- counter + 1
      hit_name <- colnames(g_mat)[i]
      # First heatmap
      single_hit <- g_mat[ , i, drop = FALSE]
      # Second/Third: Geno recon or trans & Pheno recon or trans
      #if (recon_or_trans == "trans"){
      #  filename <- paste(filename_start, "tree_sig_hit_", counter , ".pdf", sep = "")
      #} else {
      #  filename <- paste(filename_start, "tree_sig_hit_", counter , ".pdf", sep = "")
      #}
      #pdf(filename, width = 16, height = 20)
      par(mfrow = c(1, 2))
      par(mgp = c(3, 1, 0))
      par(oma = c(0, 0, 4, 0))

      # Genotype
      edge_color <- geno_trans_or_recon[[i]]
      edge_color[edge_color == 1] <- "orange"
      edge_color[edge_color == 0] <- "black"
      edge_color[all_high_conf_edges[[i]] == 0] <- "grey"

      if (recon_or_trans == "trans"){
        title <- paste("Genotype transition edges\nBlack=no trans Orange=trans\nGrey=low confidence", sep = "")
      } else {
        title <- paste("Genotype reconstruction\nBlack=absent Orange=present\nGrey=low confidence", sep = "")
      }
      plot(tr,
           main = title,
           type = "phylogram",
           use.edge.length = TRUE,
           edge.width = 2.0,
           cex = 0.5,
           cex.main = 1.0,
           no.margin = FALSE,
           label.offset = 0.0001,
           edge.color = edge_color)

    # Phenotype
    edge_color <- pheno_recon_or_trans_by_edge
    edge_color[edge_color == 1] <- "green"
    edge_color[edge_color == 0] <- "black"
    edge_color[all_high_conf_edges[[i]] == 0] <- "grey"
    if(recon_or_trans == "trans"){
      title <- paste("Phenotype transition edges\n Black=no trans Green=trans\nGrey=low confidence", sep = "")
    } else {
      title <- paste("Phenotype reconstruction\n Black=absent Green=present\nGrey=low confidence", sep = "")
    }
    plot(tr,
         main = title,
         type = "phylogram",
         use.edge.length = TRUE,
         edge.width = 2.0,
         cex = 0.5,
         cex.main = 1.0,
         no.margin = FALSE,
         label.offset = 0.0001,
         edge.color = edge_color)

     # mtext(paste("Sig hit: ", hit_name, sep = ""),
    #        outer = TRUE,
    #        side = 3,
    #        cex = 1.2,
    #        line = 1)
    #  dev.off()
    } # end subset on significant
  } # end for loop
} # end plot_sig_hits_summary()

create_heatmap_compatible_tree <- function(tree){
  heatmap_tree <- tree
  heatmap_tree$edge.length[which(heatmap_tree$edge.length == 0)] <- 0.00001
  heatmap_tree <- chronopl(heatmap_tree,
                           lambda = 0.1,
                           tol = 0)
  heatmap_tree <- as.dendrogram(as.hclust.phylo(heatmap_tree))
  return(heatmap_tree)
} # end create_heatmap_compatible_tree()

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
      set.seed(2)
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
  dropped_genotype_names <- colnames(geno)[!keepers]
  return(dropped_genotype_names)
} # end get_dropped_genotypes

report_num_high_confidence_trans_edge <- function(genotype_transition, high_conf_edges, geno_names){
  num_high_confidence_transition_edges <- rep(0, length(high_conf_edges))
  for (p in 1:length(high_conf_edges)){
    num_high_confidence_transition_edges[p] <- sum(genotype_transition[[p]]$transition * high_conf_edges[[p]])
  }

  names(num_high_confidence_transition_edges) <- geno_names
  return(num_high_confidence_transition_edges)
} # end report_num_high_confidence_trans_edge



discrete_plots <- function(tr, dir, name, a,
                           annot, num_perm, recon_hit_vals,
                           trans_hit_vals, p_recon_edges,
                           g_recon_edges, pheno_anc_rec,
                           recon_perm_obs_results, trans_perm_obs_results,
                           tr_and_pheno_hi_conf, geno_confidence,
                           g_trans_edges, p_trans_edges, snp_in_gene){

  pdf(paste0(dir, "/phyc_",  name, ".pdf"))
  # reconstruction first
  par(mfrow = c(1,1))
  make_manhattan_plot(dir, name, recon_hit_vals, a, "reconstruction")

  # TODO 2019-03-18 fix all references to reconstruction a genotype transition-- because that's what it should be
  g_recon_mat <- matrix(0, nrow = Nedge(tr), ncol = length(g_recon_edges))

  for (i in 1:length(g_recon_edges)){
    g_recon_mat[ , i] <- g_recon_edges[[i]]
    g_recon_mat[geno_confidence[[i]] == 0, i] <- NA
  }

  p_recon_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
  p_mat <- matrix(p_recon_edges, nrow = length(p_recon_edges), ncol = 1)
  colnames(p_mat) <- "pheno_presence"
  phenotype_annotation <- as.data.frame(p_mat)
  row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)

  temp_g_recon_mat <- cbind(phenotype_annotation, g_recon_mat)
  g_recon_mat <- temp_g_recon_mat[order(temp_g_recon_mat[,1], na.last = FALSE, decreasing = FALSE ), 2:ncol(temp_g_recon_mat), drop = FALSE]


  cell_width_value <- 1.5
  if (ncol(g_recon_mat) < 50){
    cell_width_value <- 10
  }

  colnames(g_recon_mat) <- row.names(recon_hit_vals)

  significant_loci <- data.frame("locus" = rep("not_sig", ncol(g_recon_mat)), stringsAsFactors = FALSE)
  row.names(significant_loci) <- row.names(recon_hit_vals)
  log_p_value <- data.frame(-log(recon_hit_vals))
  significant_loci[log_p_value > -log(a)] <- "sig"

  if (!is.null(snp_in_gene)){
    snp_in_gene <- as.data.frame(snp_in_gene, row.names = 1)
    colnames(snp_in_gene) <- "SNPs in gene"
    column_annot <- cbind(significant_loci, log_p_value, snp_in_gene) # TODO add a test to make sure order doesn't matter here
  } else {
    column_annot <- cbind(significant_loci, log_p_value)
  }


  ann_colors = list(
    locus = c(not_sig = "white", sig = "blue"),
    pheno_presence = c( na = "grey", absence = "white", presence = "red")
  )

  ordered_by_p_val      <-           g_recon_mat[ , match(row.names(log_p_value)[order(log_p_value[ , 1])], colnames(g_recon_mat))]
  column_annot_ordered_by_p_val <-     column_annot[match(row.names(log_p_value)[order(log_p_value[ , 1])], row.names(column_annot)), ]

  # reconstruction loci summary heat maps
  pheatmap( # Plot the heatmap
     ordered_by_p_val,
     main          = paste0("Edges:\n Genotype transition with phenotype presence/absence"),
     cluster_cols  = FALSE,
     na_col = "grey",
     cluster_rows  = FALSE,
     show_rownames = FALSE,
     color = c("white", "black"),
     annotation_col = column_annot_ordered_by_p_val,
     annotation_row = phenotype_annotation,
     annotation_colors = ann_colors,
     show_colnames = TRUE,
     cellwidth = cell_width_value)

  # loop through reconstruction sig hits:
  pheno_as_list <- list(p_recon_edges)
  pheno_conf_as_list <- list(tr_and_pheno_hi_conf)
  # TODO break these plots into more functions b/c lots of redundant code between recon and transition plots
  for (j in 1:nrow(recon_hit_vals)){
    if (recon_hit_vals[j, 1] < a){
      par(mfrow = c(3, 2), mgp = c(3, 1, 0), oma = c(0, 0, 4, 0), mar = c(4, 4, 4, 4))
      # pheno
      plot_tree_with_colored_edges(tr, pheno_as_list, pheno_conf_as_list, "grey", "red", paste0("\n Phenotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", 1)
      # geno
      plot_tree_with_colored_edges(tr, g_recon_edges, geno_confidence, "grey", "red", paste0(row.names(recon_hit_vals)[j], "\n Genotype transition:\n Red = transition; Black = no transition"), annot, "recon", j)
      # Permutation test
      max_x <- max(recon_perm_obs_results$permuted_count[[j]], recon_perm_obs_results$observed_overlap[j]) # change to loop through sig hits
      hist(recon_perm_obs_results$permuted_count[[j]],
           breaks = num_perm/10,
           xlim = c(0, max_x),
           col = "grey",
           border = FALSE,
           ylab = "Count",
           xlab = "# edges where genotype-phenotype co-occur",
           main = paste0("Overlap of genotype transition edge\n& phenotype presence \npval=", round(recon_hit_vals[j, 1], 4), "\nRed=observed,Grey=permutations", sep = ""))
      abline(v = recon_perm_obs_results$observed_overlap[j], col = "red")

      p_recon_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
      p_mat <- matrix(p_recon_edges, nrow = length(p_recon_edges), ncol = 1)
      colnames(p_mat) <- "pheno_presence"
      phenotype_annotation <- as.data.frame(p_mat)
      row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)


      temp_g_recon_edges <- g_recon_edges[[j]]
      temp_g_recon_edges[geno_confidence[[j]] == 0] <- NA
      g_mat <- as.matrix(temp_g_recon_edges)
      row.names(g_mat) <- c(1:nrow(g_mat))
      colnames(g_mat) <- "genotype_transition"
      temp_g_mat <- cbind(g_mat, phenotype_annotation)
      g_mat <- temp_g_mat[order(temp_g_mat[,2], temp_g_mat[,1], na.last = FALSE, decreasing = FALSE ), 1, drop = FALSE]

      plotting_logical <- check_if_g_mat_can_be_plotted(g_mat)

      if (plotting_logical){
        ann_colors <- make_ann_colors(g_mat)

        if (!is.null(snp_in_gene)){
          num_snps <- snp_in_gene[row.names(recon_hit_vals)[j], , drop = FALSE]
          row.names(num_snps) <- "genotype_transition"
          colnames(num_snps) <- "SNPs_in_gene"
          ann_colors <- c(ann_colors, list(SNPs_in_gene = c(num_snps_in_gene = "blue")))
        } else {
          num_snps <- NULL
        }

        pheatmap(
          mat               = g_mat,
          main              = paste0(row.names(recon_hit_vals)[j], "\n Tree edges clustered by edge type\n Genotype transition edge\n & phenotype present edge"),
          cluster_cols      = FALSE,
          cluster_rows      = FALSE,
          na_col            = "grey",
          show_rownames     = FALSE,
          color             = c("white", "red"),
          annotation_row    = phenotype_annotation,
          annotation_col    = num_snps,
          annotation_legend = TRUE,
          annotation_colors = ann_colors,
          show_colnames     = TRUE,
          legend            = TRUE,
          cellwidth         = 20)
      }
    }
  }

  par(mfrow = c(1,1))
  make_manhattan_plot(dir, name, trans_hit_vals, a, "transition")

  # start transition heatmaps
  g_trans_mat <- matrix(0, nrow = Nedge(tr), ncol = length(g_recon_edges))

  for (i in 1:length(g_recon_edges)){
    g_trans_mat[ , i] <- g_trans_edges[[i]]
    g_trans_mat[geno_confidence[[i]] == 0, i] <- NA
  }

  p_trans_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
  p_mat <- matrix(p_trans_edges, nrow = length(p_trans_edges), ncol = 1)
  colnames(p_mat) <- "pheno_transitions"
  phenotype_annotation <- as.data.frame(p_mat)
  row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)

  temp_g_trans_mat <- cbind(phenotype_annotation, g_trans_mat)
  g_trans_mat <- temp_g_trans_mat[order(temp_g_trans_mat[,1], na.last = FALSE, decreasing = FALSE ), 2:ncol(temp_g_trans_mat), drop = FALSE]
  colnames(g_trans_mat) <- row.names(trans_hit_vals)

  significant_loci <- data.frame("locus" = rep("not_sig", ncol(g_trans_mat)), stringsAsFactors = FALSE)
  row.names(significant_loci) <- row.names(trans_hit_vals)
  log_p_value <- data.frame(-log(trans_hit_vals))
  significant_loci[log_p_value > -log(a)] <- "sig"


  if (!is.null(snp_in_gene)){
    column_annot <- cbind(significant_loci, log_p_value, snp_in_gene) # TODO add a test to make sure order doesn't matter here
  } else {
    column_annot <- cbind(significant_loci, log_p_value)
  }

  ann_colors = list(
    locus = c(not_sig = "white", sig = "blue"),
    pheno_transitions = c( na = "grey", no_transition = "white", transition = "red")
  )

  ordered_by_p_val      <-           g_trans_mat[ , match(row.names(log_p_value)[order(log_p_value[ , 1])], colnames(g_trans_mat))]
  column_annot_ordered_by_p_val <-     column_annot[match(row.names(log_p_value)[order(log_p_value[ , 1])], row.names(column_annot)), ]

  # Transition loci summary heat maps
  pheatmap( # Plot the heatmap
    ordered_by_p_val,
    main          = paste0("Edges:\n Genotype transitions with phenotype transitions"),
    cluster_cols  = FALSE,
    na_col = "grey",
    cluster_rows  = FALSE,
    show_rownames = FALSE,
    color = c("white", "black"),
    annotation_col = column_annot_ordered_by_p_val,
    annotation_row = phenotype_annotation,
    annotation_colors = ann_colors,
    show_colnames = TRUE,
    cellwidth = cell_width_value)

  # loop through transition sig hits:
  for (j in 1:nrow(trans_hit_vals)){
    if (trans_hit_vals[j, 1] < a){
      par(mfrow = c(3, 2), mgp = c(3, 1, 0), oma = c(0, 0, 4, 0), mar = c(4, 4, 4, 4))
      # Plot pheno
      p_trans_edges_as_list <- list(p_trans_edges)
      plot_tree_with_colored_edges(tr, pheno_as_list,         pheno_conf_as_list, "grey", "red", paste0("\n Phenotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", 1)
      plot_tree_with_colored_edges(tr, p_trans_edges_as_list, pheno_conf_as_list, "grey", "red", paste0("\n Phenotype transitions:\n Red = transition; Black = no change"), annot, "recon", 1)
      # Plot geno
      plot_tree_with_colored_edges(tr, g_recon_edges, geno_confidence, "grey", "red", paste0(row.names(recon_hit_vals)[j], "\n Genotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", j)
      plot_tree_with_colored_edges(tr, g_trans_edges, geno_confidence, "grey", "red", paste0(row.names(trans_hit_vals)[j], "\n Genotype transitions:\n Red = transition; Black = no change"), annot, "recon", j)
      # Permutation test
      max_x <- max(trans_perm_obs_results$permuted_count[[j]], trans_perm_obs_results$observed_overlap[j]) # change to loop through sig hits
      hist(trans_perm_obs_results$permuted_count[[j]],
           breaks = num_perm/10,
           xlim = c(0, max_x),
           col = "grey",
           border = FALSE,
           ylab = "Count",
           xlab = "# edges where genotype-phenotype transitions co-occur",
           main = paste0("Geno & pheno transition overlap\npval=", round(trans_hit_vals[j, 1], 4), "\nRed=observed,Grey=permutations", sep = "")) # TODO add rank pvalue
      abline(v = trans_perm_obs_results$observed_overlap[j], col = "red")

      # edge heatmap - heatmap is tree edges, annotation is phenotype edges
      par(mfrow = c(1,1))
      p_trans_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
      p_mat <- matrix(p_trans_edges, nrow = length(p_trans_edges), ncol = 1)
      colnames(p_mat) <- "pheno_transition"
      phenotype_annotation <- as.data.frame(p_mat)
      row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)

      temp_g_trans_edges <- g_trans_edges[[j]]
      temp_g_trans_edges[geno_confidence[[j]] == 0] <- NA
      g_mat <- as.matrix(temp_g_trans_edges)
      row.names(g_mat) <- c(1:nrow(g_mat))
      colnames(g_mat) <- "genotype_transition"
      temp_g_mat <- cbind(g_mat, phenotype_annotation)
      g_mat<- temp_g_mat[order(temp_g_mat[,2], temp_g_mat[,1], na.last = FALSE, decreasing = FALSE ), 1, drop = FALSE]

      ann_colors = list(pheno_transition = c( na = "grey", no_transition = "white", transition = "red"))

      plotting_logical <- check_if_g_mat_can_be_plotted(g_mat)


      if (plotting_logical){
        ann_colors <- make_ann_colors(g_mat)
        if (!is.null(snp_in_gene)){
          num_snps <- snp_in_gene[row.names(trans_hit_vals)[j], , drop = FALSE]
          row.names(num_snps) <- "genotype_transition"
          colnames(num_snps) <- "SNPs_in_gene"
          ann_colors <- c(ann_colors, list(SNPs_in_gene = c(num_snps_in_gene = "blue")))
        }

        pheatmap(
          g_mat,
          main              = paste0(row.names(trans_hit_vals)[j], "\n Tree edges: genotype & phenotype transitions"),
          cluster_cols      = FALSE,
          cluster_rows      = FALSE,
          na_col            = "grey",
          show_rownames     = FALSE,
          color             = c("white", "red"),
          annotation_row    = phenotype_annotation,
          annotation_col    = num_snps,
          annotation_colors = ann_colors,
          annotation_legend = TRUE,
          show_colnames     = TRUE,
          legend = FALSE,
          cellwidth         = 20)
      }
    }
  }
  dev.off()
} # end discrete_plots()


gene_test_from_snps <- function(){
# goal of this is to rebuild the gene test from snps and keep snp annotation
} # end gene_test_from_snps

build_gene_anc_recon_and_conf_from_snp <- function(tr, geno, g_reconstruction_and_confidence, gene_to_snp_lookup_table){

  # VALIDATE INPUTS -----------------------------------------------------------#
  check_for_root_and_boostrap(tr)
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
  # TODO test this function works.
  gene_list_built_from_snps_just_node_anc_rec <- rep(list(0), length(gene_list$tip_and_node_recon))
  for (m in 1:length(gene_list$tip_and_node_recon)){
    gene_list_built_from_snps_just_node_anc_rec[[m]] <- gene_list[[m]]$tip_and_node_recon[(Ntip(tr) + 1):(Ntip(tr) + Nedge(tr))]
  }
  return(gene_list_built_from_snps_just_node_anc_rec)
} # end build_node_anc_recon_from_gene_list()


build_gene_trans_from_snp_trans <- function(tr, geno, geno_transition, gene_to_snp_lookup_table){

  # VALIDATE INPUTS -----------------------------------------------------------#
  check_for_root_and_boostrap(tr)
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

check_if_g_mat_can_be_plotted <- function(geno_matrix){
  ones <- sum(geno_matrix == 1, na.rm = TRUE) > 1
  zeros <- sum(geno_matrix == 0, na.rm = TRUE) > 1
  nas <- sum(is.na(geno_matrix)) > 1

  plot_logical <- FALSE #
  if (ones == 1 && zeros == 1 && nas == 0) {
    plot_logical <- TRUE
  }
  if (ones + zeros + nas == 3) {
    plot_logical <- TRUE
  }
  return(plot_logical)
}

make_ann_colors <- function(geno_matrix){
  ones <- sum(geno_matrix == 1, na.rm = TRUE) > 1
  zeros <- sum(geno_matrix == 0, na.rm = TRUE) > 1
  nas <- sum(is.na(geno_matrix)) > 1

  if (ones + zeros + nas == 3){
    ann_colors = list(pheno_presence = c(na = "grey", absent = "white", present = "red"))
  } else if (ones == 1 && zeros == 1 && nas == 0) {
    ann_colors = list(pheno_presence = c(absent = "white", present = "red"))
  } else if (ones == 1 && zeros == 0 && nas == 1) {
    ann_colors = list(pheno_presence = c(na = "grey", present = "red"))
  } else if (ones == 0 && zeros == 1 && nas == 1) {
    ann_colors = list(pheno_presence = c(na = "grey", absent = "white"))
  } else if (ones == 0 && zeros == 0 && nas == 1) {
    ann_colors = list(pheno_presence = c(na = "grey"))
  } else if (ones == 0 && zeros == 1 && nas == 0) {
    ann_colors = list(pheno_presence = c(absent = "white"))
  } else if (ones == 1 && zeros == 0 && nas == 0) {
    ann_colors = list(pheno_presence = c(present = "red"))
  } else {
    stop("no ones, zeroes, or NAs present in g_mat")
  }
  return(ann_colors)
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


pick_recon_model <- function(mat, tr, disc_cont, num, recon_method){
  # Note, SYMreconstruction removed because SYM == ER for binary inputs.
  # Use this function to choose the best model for reconstruction.

  # Check inputs --------------------------------------------------------------#
  if (disc_cont != "discrete"){
    stop("Only pick recon model for discrete characters. Continuous characters must be BM.")
  }

  check_if_binary_matrix(mat)

  # Function ------------------------------------------------------------------#
  alpha <- 0.05
  significant_difference_in_AIC <- 2

  # Reference for model testing: https://www.r-phylo.org/wiki/HowTo/Ancestral_State_Reconstruction & http://blog.phytools.org/2017/07/comparing-fitted-discrete-character.html
  # Test ER vs ARD
  ERreconstruction  <- ace(mat[ , num, drop = TRUE], tr, type = disc_cont, method = recon_method, marginal = FALSE, model = "ER")

  # Some ARD models don't work well with the data and given a warning message like:  "In sqrt(diag(solve(h))) : NaNs produced"
  # To ensure the ER model is prefered in this case use the following warning catching:
  error_msg <- "ARD_bad_fit"
  ARDreconstruction <- tryCatch(ace(mat[ , num, drop = TRUE], tr, type = disc_cont, method = recon_method, marginal = FALSE, model = "ARD"), warning = function(x) {error_msg})

  # If ARD gave a warning, pick ER
  best_model <- "ER"
  if (length(ARDreconstruction) == 1){
    best_model <- "ER"
  } else { # Pick best of ER and ARD
    p_ER_ARD  <- 1 - pchisq(2*abs(ERreconstruction$loglik - ARDreconstruction$loglik), 1)
    if (p_ER_ARD < alpha & AIC(ERreconstruction) > (significant_difference_in_AIC + AIC(ARDreconstruction))){
      best_model <- "ARD"
    }
    kER  <- length(ERreconstruction$rates)
    kARD <- length(ARDreconstruction$rates)
    new_p_ER_ARD  <- pchisq(-2*(logLik(ERreconstruction)-logLik(ARDreconstruction)),  df = kARD - kER,  lower.tail = FALSE)

    num_digits <- 4
    if (round(p_ER_ARD, num_digits) != round(new_p_ER_ARD, num_digits)){
      stop("ER_ARD loglikelihood test should be the same from both calculations")
    }
  }

  return(best_model)
} # end pick_recon_model()

find_parent_edge <- function(tr, edge_num){
  # given an edge number, get edge number of parent edge
  parent_node <- tr$edge[edge_num, 1]
  parent_edge <- which(tr$edge[ , 2] == parent_node)

  return(parent_edge)
  # TODO it breaks on 1st edge, need to deal with that, but not sure what's best yet
} # end find_parent_edge

check_if_phenotype_normal <- function(pheno, continuous_or_discrete){
  if (continuous_or_discrete == "continuous"){
    result <- shapiro.test(unlist(pheno))
    alpha <- 0.05
    if (result$p < alpha){
      print("Consider normalizing your phenotype")
    }
  }
} # end check_if_phenotype_normal

check_if_convergence_occurs <- function(pheno, tree, continuous_or_discrete){
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
  # VALIDATION
  if (length(genotype_transition_by_edges[[1]]$transition) != Nedge(tr)){
    stop("Dimension mismatch")
  }
  check_for_root_and_boostrap(tr)
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
