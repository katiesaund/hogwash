run_ks_test <- function(t_index, non_t_index, phenotype_by_edges){
  # TODO deal with cases when there isn't enough data (aren't at least 1 of least t_index and non_t_index)
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
  if(!is.vector(t_index)){stop("Transition index must be a vector")}
  check_is_number(t_index[1])
  if(max(t_index) > nrow(phenotype_by_edges)){stop("Transition index must represent a tree edge")}

  if(!is.vector(non_t_index)){stop("Non-transition index must be a vector")}
  check_is_number(non_t_index[1])
  if(max(non_t_index) > nrow(phenotype_by_edges)){stop("Non-transition index must represent a tree edge")}

  if (length(t_index) < 1 | length(non_t_index) < 1){
    stop("Not enough high confidence transition edges to use for KS test.")
  }
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

calculate_genotype_significance <- function(mat, permutations, genotype_transition_list, tr, pheno_recon_ordered_by_edges, genotype_confidence, genotype_reconstruction){
  # Function description -------------------------------------------------------
  #
  # Input:
  # mat is a genotype_matrix. Each row is a genotype. Each column is a 0/1 value on an edge.
  # permutations is a numeric integer.
  # genotype_transition_list is a list of length = ncol(mat) and each list has length = Nedge(tr). All entries are 0 or 1.
  # tr object of phylo.
  # phenotype_recon_ordered_by_edges is a matrix. Nrow = Nedge(tr). Ncol =2. where the phenotype reconstruction is ordered by edges. IN THE SAME FORMAT AS tr$EDGE. SO NODE VALUES WILL APPEAR MULTIPLE TIMES IN THE tr.
  #            THIS FORMAT WILL MAKE IT MUCH EASIER TO CALCULATE PHENOTYPE CHANGE ON EDGES.
  # genotype_confidence is a list of length = ncol(mat) and each list has length = Nedge(tr). All entries are 0 (low confidence) or 1 (high confidence).
  #            NOTE: genotype_confidence lists the confidence in each edge. High confidence means the edge is high confidence by genotype reconstruction, phenotype reconstruction, bootstrap value, and edge length.
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



get_sig_hits_while_correcting_for_multiple_testing <- function(hit_values, fdr){
  # Function description -------------------------------------------------------
  # Given p-values that have not yet been corrected for multiple testing, apply
  # a false discovery rate. Using FDR instead of FWER/bonferroni because these
  # approaches tend to suffer from low power. FDR control increases power while
  # bounding error.
  # TODO pretty sure this is too stringent as currently implemented-- not really doing FDR. TAlk to Krishna. I should have to set the false discovery rate somewhere.
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
