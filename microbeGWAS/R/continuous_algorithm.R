run_ks_test <- function(t_index, non_t_index, phenotype_by_edges){
  # TODO deal with cases when there isn't enough data (aren't at least 1 of least t_index and non_t_index)
  # TODO -- If i were really fancy -- I would remove these from the genotype and then report that they're bad and continue without breaking out of the algorithm
  # right now it will cause a stop error --- which won't be good for lots of things.
  # Can't seem to recapitulate the error when there are not enough samples to run the test,but all indices are list of at least one. Hmm.

  # Function description -------------------------------------------------------
  # Run a Kolmogorov-Smirnov test on the continuous phenotype. The phenotype is
  # divided into two groups: transition edges and non-transition edges.
  #
  # Input:
  # t_index.     Vector. Each number is the index phenotype transition edges on the tree.
  # non_t_index. Vector.  Each number is the index phenotype transition edges on the tree.
  # phenotype_by_edges. Numeric matrix. Phenotype values are each node stored in a matrix. Each row corresponds to an edge. 1st column is parent node. 2nd column is daughter node.
  #
  # Outputs:
  # A list of 4 elements:
  #   $pval. Numeric. KS Test p-value.
  #   $statistic. Numeric. KS Test statistic.
  #   $pheno_trans_delta. Numeric vector. The value of the delta of the phenotype on transition edges.
  #   $pheno_non_trans_delta. Numeric vector. The value of the delta of the phenotype on all non-transition edges.
  #
  # Check input ----------------------------------------------------------------
  if (length(t_index) < 1 | length(non_t_index) < 1){
    stop("Not enough high confidence transition edges to use for KS test.")
  }

  if (!is.vector(t_index)){stop("Transition index must be a vector")}
  check_is_number(t_index[1])
  if (max(t_index) > nrow(phenotype_by_edges)){stop("Transition index must represent a tree edge")}

  if (!is.vector(non_t_index)){stop("Non-transition index must be a vector")}
  check_is_number(non_t_index[1])
  if (max(non_t_index) > nrow(phenotype_by_edges)){stop("Non-transition index must represent a tree edge")}


  # Function -------------------------------------------------------------------
  p_trans_delta     <- calculate_phenotype_change_on_edge(t_index,     phenotype_by_edges)
  p_non_trans_delta <- calculate_phenotype_change_on_edge(non_t_index, phenotype_by_edges)
  set.seed(1)
  ks_results        <- ks.test(p_trans_delta, p_non_trans_delta)

  # Return output --------------------------------------------------------------
  results <- list("pval"      = round(ks_results$p.value, digits = 20),
                  "statistic" = round(ks_results$statistic, digits = 20),
                  "pheno_trans_delta"     = p_trans_delta,
                  "pheno_non_trans_delta" = p_non_trans_delta)
  return(results)
} # end run_ks_test()

get_hi_conf_tran_indices <- function(geno_tran, geno_conf, index, tr){
  # Function description -------------------------------------------------------
  #
  # Input:
  # geno_tran. List of lists. Length(genotype_transition_list) == number of genotypes.
  #            For each genotype there are two lists: $transition and $trans_dir
  #            Length($transition) == length($trans_dir) == Nedge(tr).
  #            Each of these are either 0, 1, or -1.
  # tr. Phylo.
  # index. Number. Just i from the loop outside of this function.
  # geno_conf.  List. Length(list) = ncol(mat) == number of genotypes.
  #                   Each entry is a vector with length == Nedge(tr). All
  #                   entries are 0 (low confidence) or 1 (high confidence).
  #
  # Outputs:
  # results. List of 2.
  #   "trans_index" = hi_conf_trans_index. Vector of numbers. Each number
  #                                        corresponds to an edge that is both a
  #                                       transition and high confidence.
  #   "non_trans_index" = hi_conf_non_trans_index. Vector of numbers. Each number
  #                                        corresponds to an edge that is both a
  #                                        not a transition and high confidence.
  #
  # Check input ----------------------------------------------------------------
  if (length(geno_tran) != length(geno_conf)){
    stop("One genotype should correspond to each entry of geno_tran and geno_conf")
  }
  check_is_number(index)
  if (index > length(geno_tran) | index < 1){
    stop("Index is a counter from 1:number of genomes.")
  }
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)


    # Function -------------------------------------------------------------------
  # GRAB THE IDS OF THE TRANSITION EDGES
  trans_index     <- c(1:Nedge(tr))[as.logical(geno_tran[[index]]$transition)]
  non_trans_index <- c(1:Nedge(tr))[!geno_tran[[index]]$transition]
  # [1] EX:  8 12 13 16 19 26 27 31 37 44 52 56 64 67 68 76 77 80 89 92 97 98
  # THESE EDGES ARE DEFINED BY THE NODES IN THE CORRESPONDING ROWS OF tr$EDGE

  # SUBSET TRANSITION / NON-TRANSITION EDGES TO ONLY HIGH CONFIDENCE ONES
  hi_conf_trans_index     <- trans_index[    as.logical(geno_conf[[index]][trans_index    ])]
  hi_conf_non_trans_index <- non_trans_index[as.logical(geno_conf[[index]][non_trans_index])]

  # Return output --------------------------------------------------------------
  return(list("trans_index" = hi_conf_trans_index,
              "non_trans_index" = hi_conf_non_trans_index))
}

continuous_permutation <- function(index_obj, tr, geno_conf, perm, num_i){
  # Note on implementation of this permutation. I've tested this as a for()
  # loop, an apply statement, and using the replicate() function.
  # for() loop was more than 10x faster than the apply statement.
  # for() loop was slightly faster than the replicate function.

  # do the permutation part
  num_trans_edges          <- length(index_obj$trans_index)
  list_of_all_edges        <- c(1:Nedge(tr))
  hi_conf_edges            <- list_of_all_edges[as.logical(geno_conf[[num_i]])]
  num_hi_conf_edges        <- sum(geno_conf[[num_i]])
  if (num_hi_conf_edges != length(hi_conf_edges)){
    stop("these should be the same")
  }
  permuted_trans_index_mat <- matrix(   nrow = perm, ncol = num_trans_edges)
  set.seed(1) # for reproducability of the sample() function
  for (j in 1:perm){  # create a random sample of the tr
    permuted_trans_index_mat[j, ] <- sample(1:num_hi_conf_edges,
                                        size = num_trans_edges,
                                        replace = TRUE,
                                        prob = tr$edge.length[hi_conf_edges]/sum(tr$edge.length[hi_conf_edges]))
  } # end for (j)

  # permuted_trans_index_mat is my new, permuted "indices$trans_index" where
  # each row is a new, fake list of transition genotype branches.
  # BUT CAVEAT: these are just fake/null transitions and some of them are
  # probably actually touching! If I wanted to be super legit I would recreate
  # as many hits, calculate new transitions, and then use those in my permutation
  # test, somehow controlling for variable numbers of transitions.

  # Notes on variable names: on May 9, 2019 which_branches became hi_conf_edges
  # and all_sampled_branches became permuted_trans_index_mat to clarify what
  # each variable means.
  return(list("hi_conf_edges" = hi_conf_edges,
              "permuted_trans_index_mat" = permuted_trans_index_mat))
}

# genotype, args$perm, geno_trans, args$tree, pheno_recon_edge_mat, high_conf_ordered_by_edges, geno_recon_ordered_by_edges

calculate_empirical_pheno_delta <- function(perm, permuted_trans_mat, hi_conf_edge_nums, pheno_by_edges){
  # Calculate permuted (aka empirical) pheno deltas
  empirical_ks_stat <- rep(NA, perm)
  for (k in 1:perm){
    permuted_trans_index     <- unique(permuted_trans_mat[k, ])
    permuted_non_trans_index <- c(1:length(hi_conf_edge_nums))[!(c(1:length(hi_conf_edge_nums)) %in% unique(permuted_trans_mat[k, ]))]
    empirical_results        <- run_ks_test(permuted_trans_index, permuted_non_trans_index, pheno_by_edges)
    empirical_ks_stat[k]     <- empirical_results$statistic
  } # end for (k)
  return(empirical_ks_stat)
}

calculate_genotype_significance <- function(mat, permutations, genotype_transition_list, tr, pheno_recon_ordered_by_edges, genotype_confidence, genotype_reconstruction){
  # Function description -------------------------------------------------------
  #
  # Input:
  # mat. Matrix. Nrow(mat) == number of genotypes. Each column is a variant. Matrix is binary.
  # permutations. Integer. Number of times to run the permutation test.
  # genotype_transition_list. List of lists. Length(genotype_transition_list) == number of genotypes.
  #                           For each genotype there are two lists: $transition and $trans_dir
  #                           Length($transition) == length($trans_dir) == Nedge(tr).
  #                           Each of these are either 0, 1, or -1.
  # tr. Phylo.
  # pheno_recon_ordered_by_edges.  Matrix. Nrow = Nedge(tr). Ncol = 2. Values
  #                                   are the phenotype reconstruction at each
  #                                   node, as ordered by edges. It's the same
  #                                   organization as tr$edge.
  # genotype_confidence. List. Length(list) = ncol(mat) == number of genotypes.
  #                      Each entry is a vector with length == Nedge(tr). All
  #                      entries are 0 (low confidence) or 1 (high confidence).
  #
  # Outputs:
  # results. List of 8.
  #          $pvals. Named numeric vector. Length == number of genotypes. Values between 1 and 0. Names are genotype names.
  #          $ks_statistics. List of numeric vectors. Length of list == number of genotypes. Each vector has length == number of permutations. Values between 1 and 0.
  #          $observed_pheno_trans_delta. List of numeric vectors. Length of list == number of genotypes. Vectors are of variable length because length is the number of transition edges for that particular genotype. Vectors are numeric.
  #          $observed_pheno_non_trans_delta. List of numeric vectors. Length of list == number of genotypes. Vectors are of variable length because length is the number of non-transition edges for that particular genotype. Vectors are numeric.
  #          $trans_median. Numberic. Vector. Length = number of genotypes. Describes median delta phenotype on all transition edges.
  #          $all_edges_median. Numeric vector. Length = number of genotypes. Describes median delta phenotype on all edges.
  #          $num_genotypes. Integer. The number of genotypes.
  #          $observed_ks_stat. Numeric Vector. Length = number of genotypes. Values between 1 and 0.
  #

  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_permutation_num_valid(permutations)
  check_if_binary_matrix(mat)
  check_if_binary_vector(genotype_confidence[[1]])
  check_dimensions(mat = mat, exact_rows = Ntip(tr), min_rows = Ntip(tr), exact_cols = NULL, min_cols = 1)
  check_dimensions(mat = pheno_recon_ordered_by_edges, exact_rows = Nedge(tr), min_rows = Nedge(tr), exact_cols = 2, min_cols = 2)
  if(length(genotype_transition_list[[1]]$transition) != Nedge(tr)){
    stop("genotype$transition incorrectly formatted")
  }
  if(length(genotype_transition_list[[1]]$trans_dir) != Nedge(tr)){
    stop("genotype$trans_dir incorrectly formatted")
  }
  if(length(genotype_transition_list) != ncol(mat)){
    stop("genotype transition doesn't a list for each genotype")
  }
  if (length(genotype_confidence[[1]]) != Nedge(tr)){
    stop("genotype_confidence is incorrectly formatted")
  }

  # Function -------------------------------------------------------------------
  num_genotypes <- ncol(mat)
  pvals <- trans_median <- all_edges_median <- observed_ks_stat <- rep(NA, num_genotypes)
  names(pvals) <- colnames(mat)
  empirical_ks_stat_list <- observed_pheno_trans_delta <- observed_pheno_non_trans_delta <- rep(list(0), num_genotypes)

  for (i in 1:num_genotypes){
    indices <- get_hi_conf_tran_indices(genotype_transition_list, genotype_confidence, i, tr)

    # Run KS test to find out if the phenotype change on transition edges is significantly different from phenotype change on non-transition edges
    observed_results <- run_ks_test(indices$trans_index, indices$non_trans_index, pheno_recon_ordered_by_edges)

    # Save these for reporting / plots
    observed_pheno_non_trans_delta[[i]] <- observed_results$pheno_non_trans_delta
    observed_pheno_trans_delta[[i]]     <- observed_results$pheno_trans_delta
    observed_ks_stat[i]                 <- observed_results$statistic
    trans_median[i]                     <- median(observed_results$pheno_trans_delta)
    all_edges_median[i]                 <- median(c(observed_results$pheno_trans_delta, observed_results$pheno_non_trans_delta))

    # Create permuted transition index matrix and get edge numbersof the high confidence edges
    perm_results <- continuous_permutation(indices, tr, genotype_confidence, permutations, i)

    # Calculate permuted (aka empirical) pheno deltas
    empirical_ks_stat <- calculate_empirical_pheno_delta(permutations, perm_results$permuted_trans_index_mat, perm_results$hi_conf_edges, pheno_recon_ordered_by_edges)

    # empirical p value calculation here: (1 + more extreme observations) / (1 + permutations)
    pvals[i] <- (sum(1 + sum(empirical_ks_stat > observed_ks_stat[i]))/(permutations + 1))

    # Save these for reporting / plots
    empirical_ks_stat_list[[i]] <- empirical_ks_stat
  } # end for (i)

  # Return output --------------------------------------------------------------
  results <- list("pvals" = pvals, "ks_statistics" = empirical_ks_stat_list,
                  "observed_pheno_trans_delta" = observed_pheno_trans_delta,
                  "observed_pheno_non_trans_delta" = observed_pheno_non_trans_delta,
                  "trans_median" = trans_median, "all_edges_median" = all_edges_median,
                  "num_genotypes" = ncol(mat), "observed_ks_stat" = observed_ks_stat) # 2018-11-28
  return(results)
} # end calculate_genotype_significance()

get_sig_hits_while_correcting_for_multiple_testing <- function(hit_values, fdr){
  # Function description -------------------------------------------------------
  # Given p-values that have not yet been corrected for multiple testing, apply
  # a false discovery rate. Using FDR instead of FWER/bonferroni because these
  # approaches tend to suffer from low power. FDR control increases power while
  # bounding error.
  #
  # Inputs:
  # hit_values: Numeric. Vector of the empirical p-value for each genotype. Each entry corresponds to respective column in genotype_matrix. Should be between 0 and 1.
  # fdr: Numeric. False discovery rate (typically ~10% for discovery work). Between 0 and 1.
  #
  # Output:
  # pvals.  List of 2.
  #         $hit_pvals. Dataframe. 1 column. Nrow = number of genotypes. Row.names = genotypes. Column name = "fdr_corrected_pvals". Values between 1 and 0.
  #         $sig_pvals. Dataframe. 1 column. Nrow = number of genotypes that are significant after FDR correction. Column name = "fdr_corrected_pvals[fdr_corrected_pvals < fdr]". Row.names = genotypes. Nrow = is variable-- could be between 0 and max number of genotypes. It will only have rows if the corrected p-value is less than the fdr value.
  #
  # Check inputs ---------------------------------------------------------------
  check_num_between_0_and_1(fdr)

  # Function -------------------------------------------------------------------
  fdr_corrected_pvals            <- p.adjust(hit_values, method = "fdr")
  sig_pvals                      <- as.data.frame(fdr_corrected_pvals[fdr_corrected_pvals < fdr])
  sig_pvals                      <- as.data.frame(sig_pvals)
  row.names(sig_pvals)           <- names(hit_values)[fdr_corrected_pvals < fdr]
  fdr_corrected_pvals            <- as.data.frame(fdr_corrected_pvals)
  row.names(fdr_corrected_pvals) <- names(hit_values)

  # Check and return output ----------------------------------------------------
  results                        <- list("hit_pvals" = fdr_corrected_pvals, "sig_pvals" = sig_pvals)
  return(results)
} # end get_sig_hits_while_correcting_for_multiple_testing
