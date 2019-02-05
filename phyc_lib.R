# Katie Saund

# TO DO 
# fix delta phenotype distribution 


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
library(ape)       # ape::ace function (ancestral reconstruction)
library(phytools)  # phylogenetic tree function library
library(ComplexHeatmap) # to make final plots for discrete phenotypes
library(phangorn)
library(pheatmap) # plots for continuous phenotypes
library(grid) # plots for continuous phenotypes
library(gridExtra) # plots for continuous phenotypes
library(ggplot2) # plots for continuous phenotypes

# library(truncnorm) # truncnorm::rtruncnorm function (creates truncated normal distribution) # 2018-08-29 started getting a bug about truncnorm loading....will haveto figure that out later. 

# FUNCTIONS FOR PHYC ----------------------------------------------------------#
ancestral_reconstruction_by_ML <- function(tr, mat, num, disc_cont, confidence_threshold){ 
  # TO DO: 
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
  reconstruction <- integer(Nedge(tr)) # initialize vector of zeroes of length Nedge(tr)
  
  if (disc_cont == "continuous"){
    
    # RECONSTRUCTION -----------------------------------------------------------
    set.seed(3)
    fitER <- ace(mat[ , num, drop = TRUE], tr, type = disc_cont, method = recon_method) 
    ML_anc_rec <- fitER$ace # vector containing reconstruction data for all internal nodes (#51 -99) where tips are 1-50.
    tip_and_node_recon <- c(mat[ , num, drop = TRUE], ML_anc_rec)
    names(tip_and_node_recon) <- c(1:sum(Ntip(tr), Nnode(tr)))
      

    
    # CONFIDENCE IN RECONSTRUCTION ---------------------------------------------
    tip_and_node_anc_rec_confidence <- rep(1, length(tip_and_node_recon))

  } else if (disc_cont == "discrete"){
    
    # RECONSTRUCTION -----------------------------------------------------------
    recon_model <- "ER" # Model:  ER.         Indicates an equal rates model of transition.
    set.seed(4)
    fitER <- ace(mat[ , num, drop = TRUE], tr, model = recon_model, type = disc_cont, method = recon_method) 
    ML_anc_rec <- as.numeric(colnames(fitER$lik.anc)[apply(fitER$lik.anc, 1, which.max)]) # Extract the mostly likely character state using which.max
    names(ML_anc_rec) <- c((Ntip(tr) + 1):(Ntip(tr) + Nnode(tr)))
    tip_and_node_recon <- c(mat[ , num, drop = TRUE], ML_anc_rec)
    names(tip_and_node_recon) <- c(1:sum(Ntip(tr), Nnode(tr)))
    
    # CONFIDENCE IN RECONSTRUCTION ---------------------------------------------
    anc_rec_confidence <- apply(fitER$lik.anc, 1, max) # get confidence values at the nodes
    tip_and_node_anc_rec_confidence <- c(rep(1, Ntip(tr)), anc_rec_confidence) # count all tips as high confidence
    tip_and_node_anc_rec_confidence <- discretize_confidence_using_threshold(tip_and_node_anc_rec_confidence, confidence_threshold) # count anything lower than threshold as low confidence

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
  transition <- transition_direction <- parent_node <- child_node <- integer(Nedge(tr)) # initialize all as zeroes
  
  for (i in 1:Nedge(tr)){
    if (tr$edge[i, 1] > Ntip(tr) & tr$edge[i, 2] <= Ntip(tr)){ # Parent is internal node, child is a tip
      parent_node[i] <- node_recon[tr$edge[i, 1] - Ntip(tr)]   # Assign node value
      child_node[i]  <- mat[ , num][tr$edge[i, 2]]             # Assign tip value 
    }
    if (tr$edge[i, 1] > Ntip(tr) & tr$edge[i, 2] > Ntip(tr)){  # Both parent and child are internal nodes
      parent_node[i] <- node_recon[tr$edge[i, 1] - Ntip(tr)]   # Assign node value
      child_node[i]  <- node_recon[tr$edge[i, 2] - Ntip(tr)]   # Assign node value
    }
    
    transition[i] <- sum(parent_node[i] + child_node[i])
    # transition[i] is either 0, 1, or 2 for discrete traits
    # transition[i] is not to be used when the trait is continuous because all, or very nearly all edges are transition edges.
    
    if (disc_cont == "discrete"){
      if (parent_node[i] == 1 & child_node[i] == 0){
        transition_direction[i] <- -1
      } else if (parent_node[i] == 0 & child_node[i] == 1){
        transition_direction[i] <- 1
      } 
      transition[transition == 2] <- 0 # If parent_node and child_node had same value (1) then no transition occured. 
      
    } else if (disc_cont == "continuous"){
      if (parent_node[i] > child_node[i]){
        transition_direction[i] <- -1
      } else if (parent_node[i] < child_node[i]){
        transition_direction[i] <- 1
      } 
      
      transition <- NA # transition[i] is not to be used when the trait is continuous because all, or very nearly all edges are transition edges.
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
  # t_index = transition edges
  # not_t_index = non transition edges
  # phenotype_by_edges = change in phenotype on each edge
  
  
  # subset index to only high confidence edges: 
  p_trans_delta     <- calculate_phenotype_change_on_edge(t_index,     phenotype_by_edges)
  p_non_trans_delta <- calculate_phenotype_change_on_edge(non_t_index, phenotype_by_edges)
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

check_dimensions <- function(mat, exact_rows, min_rows, exact_cols, min_cols){
  # Function description ------------------------------------------------------- 
  # Check that the matrix is of the specified dimensions. 
  #
  # Input: 
  # mat:        Matrix. 
  # exact_rows. Numeric. Describes expected number of rows in matrix. Can be NULL. 
  # min_rows.   Numeric. Describes minimum number of rows in matrix. Must be specified.  
  # exact_cols. Numeric. Describes expected number of columns in matrix. Can be NULL. 
  # min_cols.   Numeric. Describes minimum number of columns in matrix. Must be specified.  
  #
  # Output: 
  # None. 
  #
  # Check input ----------------------------------------------------------------
  if (!is.null(exact_rows)){check_is_number(exact_rows)}
  check_is_number(min_rows)
  if (!is.null(exact_cols)){check_is_number(exact_cols)}
  check_is_number(min_cols)
  if (class(mat) != "matrix"){stop("input must be a matrix")}
  
  # Function -------------------------------------------------------------------
  if (nrow(mat) < min_rows){
    stop("matrix has too few rows")
  }
  if (ncol(mat) < min_cols){
    stop("matrix has too few columns")
  }
  if (!(is.null(exact_rows))){
    if (nrow(mat) != exact_rows){
      stop("matrix has wrong number of rows")
    }
  }
  if (!(is.null(exact_cols))){
    if (ncol(mat) != exact_cols){
      stop("matrix has wrong number of columns")
    }
  }
} # end check_dimensions()

check_if_alpha_valid <- function(a){
  # Function description ------------------------------------------------------- 
  # Check that the alpha (threshold for significance) is within a valid range (0 < alpha < 0.1). 
  #
  # Input: 
  # a. Numeric. Value of alpha. 
  #
  # Output: 
  # None. 
  #
  # Check input ----------------------------------------------------------------
  check_is_number(a)
  
  # Function -------------------------------------------------------------------
  if (a > .1 | a <= 0){
    stop("Provide a valid alpha.")
  }
} # end check_if_alpha_valid()

check_if_dir_exists <- function(dir){
  # Function description ------------------------------------------------------- 
  # Check that output directory exists so data can be saved in it. 
  #
  # Input: 
  # dir. Character. Path to output directory. 
  #
  # Output: 
  # None. 
  #
  # Check input ----------------------------------------------------------------
  check_is_string(dir)
  
  # Function -------------------------------------------------------------------
  if (!dir.exists(dir)){
    stop("The output directory indicated does not exist.")
  }
} # end check_if_dir_exists()

check_if_p_val_valid <- function(p_values){
  # Function description ------------------------------------------------------- 
  # Check that p values are valid (0 <= p values <= 1). 
  #
  # Input: 
  # p_values. Character. Path to output directory. 
  #
  # Output: 
  # None. 
  #
  # Check input ----------------------------------------------------------------
  if (!is.vector(p_values)){
    if (!is.matrix(p_values)){
      if (!is.data.frame(p_values)){
        stop("p_values is not formatted correctly")
      }
    }
  }
  
  # Function -------------------------------------------------------------------
  if (max(as.numeric(p_values)) > 1 | min(as.numeric(p_values)) < 0){
    stop("hit_values is incorrectly formatted.")
  }
} # end check_if_p_val_valid()

check_if_permutation_num_valid <- function(perm){
  # Function description ------------------------------------------------------- 
  # Check that the permutation number inidicated is valid (1 <= perm). 
  #
  # Input: 
  # perm. Number. Times to shuffle the data on the tree to create a null distribution for the permutation test.  
  #
  # Output: 
  # None. 
  #
  # Check input ----------------------------------------------------------------
  check_is_number(perm)
  
  # Function -------------------------------------------------------------------
  if (perm < 1 ||perm %% 1 != 0){
    stop("The permutation number should be a positive integer indicating the number of null distributions to create.")
  }
} # end check_if_permutation_num_valid()

check_is_string <- function(char){
  # Function description ------------------------------------------------------- 
  # Check that the input is a character string. 
  #
  # Input: 
  # char. Character.
  #
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (!is.character(char)){
    stop("Object must be a character string.")
  }
} # end check_is_string()

check_input_format <- function(pheno, tr, geno, name, dir, perm, a, annot){
  # Function description ------------------------------------------------------- 
  # Check that all of the inputs into phyC are in the valid format.  
  #
  # Input: 
  # pheno. Matrix. Phenotype. 
  # tr.    Phylo. Tree. 
  # geno.  Matrix. Genotype. 
  # name.  Character. Output name. 
  # dir.   Character. Output path. 
  # perm.  Number. Times to shuffle the data on the tree to create a null distribution for the permutation test.  
  # a.     Number. Alpha. 
  # annot. Matrix. Annotation for heatmaps. 
  
  # Output: 
  # discrete_or_continuous. Character. Either "discrete" or "continuous". Describes the input phenotype.  
  #
  # Check input ----------------------------------------------------------------
  check_dimensions(geno, Ntip(tr), 2, NULL, 2) # Genotype matrix should have same rows as tr$tip.label, at least 2 genotypes in the columns
  check_dimensions(pheno, Ntip(tr), 2, 1, 1) # Phnoetype matrix should have same rows as tr$tip.label and exactly 1 column
  check_rownames(geno, tr) # Genotype rownames should be in the same order as the tr$tip.label
  check_rownames(pheno, tr) # Phenotype rownames should be in the same order as the tr$tip.label
  check_for_NA_and_inf(geno) 
  check_for_NA_and_inf(pheno)
  check_for_root_and_boostrap(tr)
  check_if_binary_matrix(geno)
  check_is_string(name) 
  check_if_dir_exists(dir)
  check_if_permutation_num_valid(perm)
  check_if_alpha_valid(a)
  if (!is.null(annot)){
    check_dimensions(annot, Ntip(tr), 2, 2, 2)
  }
  
  # Function -------------------------------------------------------------------
  discrete_or_continuous <- assign_pheno_type(pheno)
  
  # Check and return output ----------------------------------------------------
  check_is_string(discrete_or_continuous)
  return(discrete_or_continuous)
} # end check_input_format()

check_if_vector <- function(vector){
  # Function description ------------------------------------------------------- 
  # Check that input is a list.   
  #
  # Input: 
  # vector. Vector.   
  
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (!is.vector(vector)){
    stop("Input must be a vector")
  }
} # end check_if_vector()

check_for_NA_and_inf <- function(mat){
  # Function description ------------------------------------------------------- 
  # Check that matrix contains no NAs and no +/- infinities.   
  #
  # Input: 
  # mat. Matrix.  
  #
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (class(mat) != "matrix"){
    stop("Input should be a matrix.")
  }
  if (sum(is.na(mat)) > 0){
    stop("Input matrices should not have any NA values.")
  }
  if (sum(mat == -Inf) > 0){
    stop("Inpute matrices should not have any -Inf values.")
  }
  if (sum(mat == Inf) > 0){
    stop("Inpute matrices should not have any -Inf values.")
  }
} # end check_for_NA_and_inf()

check_for_root_and_boostrap <- function(tr){
  # Function description ------------------------------------------------------- 
  # Check that phylogenetic tree is rooted and contains bootstrap values in the node labels. 
  #
  # Input: 
  # tr. Phylo.   
  #
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (class(tr) != "phylo"){
    stop("Tree must be phylo object")
  }
  if (!is.rooted(tr)){
    stop("Tree must be rooted")
  }
  if (is.null(tr$node.label)){
    stop("Tree must have boostrap values in the nodes")
  }
} # end check_for_root_and_boostrap()


check_if_binary_matrix <- function(mat){
  # Function description ------------------------------------------------------- 
  # Check that the matrix only contains values 1 or 0. 
  #
  # Input: 
  # mat. Matrix.   
  #
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (sum(!(mat %in% c(0, 1))) > 0 | class(mat) != "matrix"){
    stop("Genotype matrix should be only 1s and 0s")
  }
} # end check_if_binary_matrix()

check_file_exists <- function(file_name){
  # Function description ------------------------------------------------------- 
  # Check that the file exists. 
  #
  # Input: 
  # file_name. Character.   
  #
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (!file.exists(file_name)){
    stop("File does not exist")
  }
} # end check_file_exists()

check_rownames <- function(mat, tr){
  # Function description ------------------------------------------------------- 
  # Check that phylogenetic tree tip labels are identical to the matrix row.names.  
  #
  # Input: 
  # mat. Matrix. 
  # tr. Phylo.   
  # 
  # Output: 
  # None. 
  #
  # Check input ----------------------------------------------------------------
  if (class(mat) != "matrix" | class(tr) != "phylo"){
    stop("Inputs are incorrectly formatted.")
  }
  
  # Function -------------------------------------------------------------------
  if (sum(row.names(mat) != tr$tip.label) != 0){
    stop("Matrix must be formatted with samples in matrix in the same order as tree$tip.label.")
  } 
} # end check_rownames()

check_is_number <- function(num){
  # Function description ------------------------------------------------------- 
  # Check that input is some type of number. 
  #
  # Input: 
  # num. Number. Could be numeric, double, or integer.  
  # 
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (!is.numeric(num)){
    if (!is.integer(num)){
      if (!is.double(num)){
        stop("Must be a number")
      }
    }
  }
} # end check_is_number()

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
  #tip_and_node_anc_rec_confidence[as.numeric(tip_and_node_anc_rec_confidence)  < 0.70] <- 0 # Use a 70% confidence threshold in the ancestral reconstruction 
  #tip_and_node_anc_rec_confidence[as.numeric(tip_and_node_anc_rec_confidence) >= 0.70] <- 1
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
  # check_if_p_val_valid(hit_values)
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
  # if (nrow(fdr_corrected_pvals) > 0){check_if_p_val_valid(unlist(fdr_corrected_pvals))}
  # if (nrow(sig_pvals) > 0){          check_if_p_val_valid(unlist(sig_pvals))}
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
  # check_if_p_val_valid(pvals$hit_pvals)
  
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

plot_significant_hits <- function(tr, a, dir, name, pval_all_transition, pheno_vector, annot, perm, results_all_trans, pheno_anc_rec, geno_reconstruction, geno_confidence, geno_transition, geno, pheno_recon_ordered_by_edges){
  # Function description ------------------------------------------------------- 
  #
  # Inputs: 
  # tr.                  Phylo. 
  # a.                   Number. Alpha. 
  # dir.                 Character. Output path.
  # name.                Character. Output name. 
  # pval_all_transition. ?
  # pheno_vector.        Vector. 
  # annot.               Matrix. 
  # perm.                Number. 
  # results_all_trans.   ?. 
  # pheno_anc_rec.       ?. 
  # geno_reconstruction. ?
  # geno_confidence.     ?
  # geno_transition.     ?
  
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

  # prep for heatmap 6
  trans_edge_mat <- NULL
  for (i in 1:length(geno_transition)){
    trans_edge_mat <- cbind(trans_edge_mat, geno_transition[[i]]$transition)
  }
  colnames(trans_edge_mat) <- colnames(geno)
  
  # Update trans_edge_mat to exclude low confidence  edges, currently it includes transition edges (all high and some low confidence transitions)
  for (c in 1:ncol(trans_edge_mat)){
    trans_edge_mat[(1:Nedge(tr))[geno_confidence[[c]] == 0], c] <- NA
  }
  # end update trans_edge_mat 
  ph_trans <- abs(pheno_recon_ordered_by_edges[ , 1] - pheno_recon_ordered_by_edges[ , 2])
  
  p_trans_mat <- matrix(ph_trans, nrow = length(ph_trans), ncol = 1)
  colnames(p_trans_mat) <- "delta_pheno"
  p_trans_mat <- as.data.frame(round(p_trans_mat, 2))
  
  significant_loci <- data.frame("locus" = rep("not_sig", ncol(trans_edge_mat)), stringsAsFactors = FALSE)
  row.names(significant_loci) <- colnames(trans_edge_mat)
  significant_loci[row.names(significant_loci) %in% row.names(all_transitions_sig_hits), ] <- "sig"
  
  log_p_value <- data.frame(log(pval_all_transition$hit_pvals))
  column_annot <- cbind(significant_loci, log_p_value)
  
  save_manhattan_plot(dir, name, pval_all_transition$hit_pvals, a, "transition")
  
  row.names(p_trans_mat) <- row.names(trans_edge_mat) <- c(1:Nedge(tr))
  # end heatmap prep
  ann_colors = list(
    locus = c(not_sig = "white", sig = "blue")
  )
  
  sorted_trans_edge_mat <-           trans_edge_mat[match(row.names(p_trans_mat)[order(p_trans_mat[ , 1])], row.names(trans_edge_mat)), ]
  ordered_by_p_val      <- sorted_trans_edge_mat[ , match(row.names(log_p_value)[order(log_p_value[ , 1])], colnames(sorted_trans_edge_mat))]
  column_annot_ordered_by_p_val <-     column_annot[match(row.names(log_p_value)[order(log_p_value[ , 1])], row.names(column_annot)), ]
  
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
    cellwidth = 20, 
    file =  paste0(dir, "/phyc_", name, "_trans_edge_heatmap_with_log_p_ordered.pdf"))
  
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
    cellwidth = 20, 
    na_col = "grey", 
    file =  paste0(dir, "/phyc_", name, "_trans_edge_heatmap_with_log_p.pdf"))

  save_data_table(trans_edge_mat, dir,name, "_trans_edge_matrix.tsv") 
  save_data_table(p_trans_mat, dir, name, "_pheno_trans_edge_matrix.tsv")
  
  # ONLY MAKE THE FOLLOWING PLOTS FOR SIGNIFICANT LOCI
  counter <- 0
  for (j in 1:nrow(pval_all_transition$hit_pvals)){
    if (pval_all_transition$hit_pvals[j, 1] < a){ 

      counter <- counter + 1
      fname <- create_file_name(dir, name, paste("sig_hit_results_", counter, ".pdf", sep = ""))
      
      pdf(fname, width = 16, height = 20)
      par(mfrow = c(3, 3))
      par(mgp   = c(3, 1, 0))
      par(oma   = c(0, 0, 4, 0))
      par(mar = c(4, 4, 4, 4))
      
      plot_continuous_phenotype(tr, pheno_vector, pheno_anc_rec)
      
      plot_tree_with_colored_edges(tr, geno_reconstruction, geno_confidence, "grey", "red", "Genotype reconstruction:\n Red = Variant; Black = WT",                annot, "recon", j)
      
      plot_tree_with_colored_edges(tr, geno_transition,     geno_confidence, "grey", "red", "Genotype transition edge:\n Red = transition; Black = No transition", annot, "trans", j)
      
      histogram_all_delta_pheno_overlaid_with_high_conf_delta_pheno(p_trans_mat, geno_confidence, tr, j)
      
      histogram_abs_high_confidence_delta_pheno_highlight_transition_edges(results_all_trans, tr, j, "grey", "red") # 6. Delta phenotype histogram. Note that results_all_trans$observed_pheno_non_trans_delta[[j]] is already subset down to only high confidence edges

      histogram_raw_high_confidence_delta_pheno_highlight_transition_edges(geno_transition, geno_confidence, pheno_recon_ordered_by_edges, tr, j, "grey", "red")
     
      # 10. KS null distribution
      hist(log(results_all_trans$ks_statistics[[j]]), 
           breaks = perm/10, col = "grey", border = FALSE,
           main = paste("Null distribution of KS statistic for all transitions.\n Red = Observed KS statistic.\n p-value = ", round(pval_all_transition$hit_pvals[j, 1], 10), "\np-value rank = ", rank(pval_all_transition$hit_pvals)[j], sep = ""), 
           ylab = "Count", 
           xlab = "ln(KS statistic)", 
           xlim = c(min(log(as.numeric(results_all_trans$observed_ks_stat[j])), log(results_all_trans$ks_statistics[[j]])), 0))
      abline(v = log(as.numeric(results_all_trans$observed_ks_stat[j])), col = "red")
      
      # 7. heatmap
      pheatmap(
        sorted_trans_edge_mat[ , j, drop = FALSE],
        main          = paste0("Edges:\n hi conf trans vs delta pheno"),
        cluster_cols  = FALSE,
        cluster_rows  = FALSE,
        na_col = "grey", 
        show_rownames = FALSE,
        color = c("white", "black"),
        annotation_row = p_trans_mat,
        show_colnames = TRUE,
        cellwidth = 20)
      
      mtext(paste("Signicant hit:", row.names(pval_all_transition$hit_pvals)[j], sep = ""), 
            outer = TRUE,
            side = 3, 
            cex = 1.2, 
            line = 1)
      
      dev.off()
    }
  }
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

read_in_arguments <- function(args){
  # TO DO:
  # Update output description. 
  # Function description ------------------------------------------------------- 
  # Read in the commandline arguments to phyC_run.R
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
                     stringsAsFactors = FALSE)
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
  filename <- paste(dir, "/phyc_", name, "_genotypes_dropped_because_convergence_not_possible.txt", sep = "")
  writeLines(dropped_genotype_names, con = filename, sep = "\n")
  
  mat <- mat[ , !geno_to_drop, drop = FALSE]
  
  # Check and return output ----------------------------------------------------
  check_if_binary_matrix(mat)
  check_rownames(mat, tr)
  return(mat)
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
  #if (nrow(hits) > 0){
  #check_if_p_val_valid(hits)
  #}
  
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

save_manhattan_plot <- function(outdir, geno_pheno_name, pval_hits, alpha, trans_or_recon){
  # Create negative log p-values with arbitrary locus numbers
  neg_log_p_value <- data.frame(-log(pval_hits))
  neg_log_p_with_num <- cbind(1:nrow(neg_log_p_value), neg_log_p_value)
  colnames(neg_log_p_with_num)[1] <- "locus"
  sig_temp <-subset(neg_log_p_with_num, neg_log_p_with_num[ , 2] > -log(alpha))
  pdf(paste0(outdir, "/phyc_", trans_or_recon, "_", geno_pheno_name, "_manhattan_plot.pdf"))
  with(neg_log_p_with_num,   
       plot(x = neg_log_p_with_num[ , 1], 
            y = neg_log_p_with_num[ , 2, drop = TRUE], 
            type = "p", 
            main = paste(trans_or_recon, "phyC", geno_pheno_name, sep = " "), 
            col = rgb(0, 0, 0, 0.3), 
            pch = 19, 
            xlab = "Genetic locus",
            ylab = "-ln(p-value) after FDR" ))
  
  abline(h = -log(alpha), col = "red")
  text(x = sig_temp[ , 1], y = sig_temp[ , 2], labels = row.names(sig_temp), pos = 1, cex = 0.7)
  dev.off()
} #end save_manhattan_plot()

plot_tree_with_colored_edges <- function(tr, edges_to_highlight, geno_confidence, edge_color_na, edge_color_bright, title, annot, trans_or_recon, index){
  edge_color <- rep("black", Nedge(tr))
  if (trans_or_recon == "recon"){
    edge_color[edges_to_highlight[[index]] == 1] <- edge_color_bright
  } else if (trans_or_recon == "trans"){
    edge_color[edges_to_highlight[[index]]$transition == 1] <- "red"
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
  # TO DO: update description
  # 1) We're now workingwith ordered by edges rather than ordered by tips then noes
  # 2) rename combined/genotype confidence. it's confusing.
  # 3) stop returning combined_confidence- don't need it. 
  
  # genotype_reconstruction is a list of vectors. Each vector corresponds with 1 genotype from the geno_mat. 
  # Each entry in the vector corresponds to a node or tip. 
  # phenotype_reconstruction is a vector where each entry corresponds to a node or tip. 
  # genotype_confidence is a list of vectors. Each vector corresponds with 1 genotype from the geno_mat.
  # Each entry in the vector corresponds to a node or tip. 1 means high confidence in the node, 0 means low confidence. 
  # phenotype_confidence is a vector where each entry corresponds to a node of tip. Same encoding as above. 

  
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
  

  mtext(output_name, 
        outer = TRUE,
        side = 3, 
        cex = 1.2, 
        line = 1)
  dev.off()
  
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
      if (recon_or_trans == "trans"){
        filename <- paste(filename_start, "tree_sig_hit_", counter , ".pdf", sep = "")
      } else {
        filename <- paste(filename_start, "tree_sig_hit_", counter , ".pdf", sep = "")
      } 
      pdf(filename, width = 16, height = 20)
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
      
      mtext(paste("Sig hit: ", hit_name, sep = ""), 
            outer = TRUE,
            side = 3, 
            cex = 1.2, 
            line = 1)
      dev.off()
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

      # count only times when empirical, permuted overlap of genotype and phenotype is more common than obsered (both_present[i])
      # and when the empirical, permuted genotype does not overlap with phenotype is less common than observed (only_geno_present[i])
      new_counter <- sum((empirical_both_present >= both_present[i]) * (empirical_only_geno_present <= only_geno_present[i]))
      temp_pval <- ((new_counter + 1)/(permutations + 1))
      
      # made changes 2018-01-17
      if (sort(redistributed_both_present, decreasing = FALSE)[(alpha * permutations)] == 0 & both_present[i] == 0){
        pval <- 1
      } else if (temp_pval == 0 | temp_pval == 1){
        pval <- 2/(permutations + 1)
      } else if (temp_pval > 0.5){
        pval <- ((1 - temp_pval)  * 2)
      } else if (temp_pval <= 0.5){
        pval <- (temp_pval * 2)
      }
      
      hit_pvals[i] <- format(round(pval, 20), nsmall = 20)
    }
  }
  names(hit_pvals) <- colnames(mat)
  return(hit_pvals)
} # end calculate_hit_pvals_corrected

save_dropped_genotypes <- function(geno, keepers, save_dir, save_name){
  dropped_genotype_names <- colnames(geno)[!keepers]
  filename <- paste(save_dir, "/phyc_", save_name, "_genotypes_dropped_because_too_few_high_conf_trans_edges.txt", sep = "")
  writeLines(dropped_genotype_names, con = filename, sep = "\n")
} # end save_dropped_genotypes

report_num_high_confidence_trans_edge <- function(genotype_transition, high_conf_edges, geno_names, outdir, outname){
  num_high_confidence_transition_edges <- rep(0, length(genotype_transition))
  for (p in 1:length(genotype_transition)){
    num_high_confidence_transition_edges[p] <- sum(genotype_transition[[p]]$transition * high_conf_edges[[p]])
  }
  names(num_high_confidence_transition_edges) <- geno_names
  write.table(num_high_confidence_transition_edges, file = paste(outdir, "/phyc_", outname, "_num_high_conf_trans_edges_per_geno.tsv", sep = ""), sep = "\t", quote= FALSE, col.names = FALSE, row.names = TRUE)
} # end report_num_high_confidence_trans_edge

# END OF SCRIPT ---------------------------------------------------------------#
