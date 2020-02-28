#' Run KS Test
#'
#' @description Run a Kolmogorov-Smirnov test on the continuous phenotype. The
#'   phenotype is divided into two groups: transition edges and non-transition
#'   edges.
#' @param t_index Vector. Each number is the index phenotype transition edges on
#'   the tree.
#' @param non_t_index Vector. Each number is the index phenotype transition
#'   edges on the tree.
#' @param phenotype_by_edges Numeric matrix. Phenotype values are each node
#'   stored in a matrix. Each row corresponds to an edge. 1st column is parent
#'   node. 2nd column is daughter node.
#'
#' @return A list of 4 elements:
#'   \describe{
#'     \item{pval.}{Numeric. KS Test p-value.}
#'     \item{statistic.}{Numeric. KS Test statistic.}
#'     \item{pheno_trans_delta.}{Numeric vector. The value of the delta of the
#'     phenotype on transition edges.}
#'     \item{pheno_non_trans_delta.}{Numeric vector. The value of the delta of
#'     the phenotype on all non-transition edges.}
#'   }
#'
#' @noRd
run_ks_test <- function(t_index, non_t_index, phenotype_by_edges){
  # Check inputs ---------------------------------------------------------------
  if (length(t_index) < 1 | length(non_t_index) < 1) {
    stop("Not enough high confidence transition edges to use for KS test.")
  }
  if (!is.vector(t_index)) {
    stop("Transition index must be a vector")
  }
  check_is_number(t_index[1])
  if (max(t_index) > nrow(phenotype_by_edges)) {
    stop("Transition index must represent a tree edge")
  }
  if (!is.vector(non_t_index)) {
    stop("Non-transition index must be a vector")
  }
  check_is_number(non_t_index[1])
  if (max(non_t_index) > nrow(phenotype_by_edges)) {
    stop("Non-transition index must represent a tree edge")
  }

  # Function -------------------------------------------------------------------
  p_trans_delta <-
    calculate_phenotype_change_on_edge(t_index, phenotype_by_edges)
  p_non_trans_delta <-
    calculate_phenotype_change_on_edge(non_t_index, phenotype_by_edges)

  set.seed(1)
  ks_results <-  withCallingHandlers(stats::ks.test(p_trans_delta,
                                                    p_non_trans_delta),
                      warning = function(w) {
                        if (grepl("cannot compute exact p-value with ties",
                                  w$message))
                          invokeRestart("muffleWarning")
                      })
  # Return output --------------------------------------------------------------
  results <- list("pval"      = round(ks_results$p.value, digits = 20),
                  "statistic" = round(ks_results$statistic, digits = 20),
                  "pheno_trans_delta"     = p_trans_delta,
                  "pheno_non_trans_delta" = p_non_trans_delta)
  return(results)
}

#' Get the indices of high confidence transition edges and non-transition edges
#'
#' @description Returns a list of the indices of each high confidence transition
#'   tree edge and a list of the indices that are high confidence, but
#'   non-transition tree edges.
#'
#' @param geno_tran List of lists. Length(genotype_transition_list) == number of
#'   genotypes. For each genotype there are two lists: $transition and
#'   $trans_dir. Length($transition) == length($trans_dir) == Nedge(tr). Each of
#'   these are either 0, 1, or -1.
#' @param geno_conf List. Length(list) = ncol(mat) == number of genotypes. Each
#'   entry is a vector with length == Nedge(tr). All entries are 0 (low
#'   confidence) or 1 (high confidence).
#' @param index Number. Just i from the loop outside of this function.
#' @param tr Phylo.
#'
#' @return results. List of 2.
#'   \describe{
#'     \item{trans_index}{Vector of numbers. Each number corresponds to an edge
#'     that is both a transition and high confidence.}
#'     \item{non_trans_index}{hi_conf_non_trans_index. Vector of numbers. Each
#'     number corresponds to an edge that is both a not a transition and high
#'     confidence.}
#'   }
#' @noRd
#'
get_hi_conf_tran_indices <- function(geno_tran, geno_conf, index, tr){
  # Check input ----------------------------------------------------------------
  check_equal(length(geno_tran), length(geno_conf))
  check_is_number(index)
  if (index > length(geno_tran) | index < 1) {
    stop("Index is a counter from 1:number of genomes.")
  }
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)

  # Function -------------------------------------------------------------------
  # GRAB THE IDS OF THE TRANSITION EDGES
  trans_index <- c(1:ape::Nedge(tr))[as.logical(geno_tran[[index]]$transition)]
  non_trans_index <- c(1:ape::Nedge(tr))[!geno_tran[[index]]$transition]
  # [1] EX:  8 12 13 16 19 26 27 31 37 44 52 56 64 67 68 76 77 80 89 92 97 98
  # THESE EDGES ARE DEFINED BY THE NODES IN THE CORRESPONDING ROWS OF tr$edge

  # SUBSET TRANSITION / NON-TRANSITION EDGES TO ONLY HIGH CONFIDENCE ONES
  hi_conf_trans_index <-
    trans_index[as.logical(geno_conf[[index]][trans_index])]
  hi_conf_non_trans_index <-
    non_trans_index[as.logical(geno_conf[[index]][non_trans_index])]

  # Return output --------------------------------------------------------------
  return(list("trans_index" = hi_conf_trans_index,
              "non_trans_index" = hi_conf_non_trans_index))
}

#' Perform permutation for continuous data
#'
#' @description Perform permutation of which genotype edges are considered high
#'   confidence transitions.
#'
#' @param geno_no_tran_index_list List. Length(list) == number of genotypes.
#'   Each entry has a vector of integers which correspond to which tree edges
#'   are genotype non_transition edges.
#' @param tr Phylo.
#' @param hi_conf_index_list. List. Length(list) = number of genotypes. Each
#'   entry has a vector of integers which correspond to which tree edges are
#'   high confidence (genotype anceestral reconstruction, phenotype
#'   reconstruction, and tree (bootstrap, edge length)).
#' @param perm Integer. Number of times to run the permutation test.
#' @param num_i Integer. The i from the loop this function is set within.
#' @param pheno_delta_list List. Length(list) == number of genotypes. Each
#'   entry is the phenotype delta for each tree edge.
#'
#' @return List of 2:
#'   \describe{
#'     \item{hi_conf_edges}{List of edge numbers. Corresponds to edges in tree
#'     with high confidence.}
#'     \item{permuted_trans_index_mat}{Matrix. Nrow = number of permutations.
#'     Ncol = Number of transition edges. Each entry is an edge number.}
#'   }
#'
#' @noRd
continuous_permutation <- function(geno_no_tran_index_list,
                                   tr,
                                   hi_conf_index_list,
                                   perm,
                                   num_i,
                                   pheno_delta_list){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_is_number(perm)
  check_if_permutation_num_valid(perm)
  check_is_number(num_i)
  check_equal(length(geno_no_tran_index_list), length(hi_conf_index_list))
  check_equal(length(pheno_delta_list), length(hi_conf_index_list))
  check_equal(length(pheno_delta_list[[1]]), ape::Nedge(tr))

  # Function -------------------------------------------------------------------
  # Note on implementation of this permutation. I've tested this as a for()
  # loop, an apply statement, and using the replicate() function.
  # for() loop was more than 10x faster than the apply statement.
  # for() loop was slightly faster than the replicate function.


  # Recall we're working within one genotype (i) at a time.
  curr_delta_pheno <- pheno_delta_list[[num_i]]
  curr_hi_conf_geno_pheno_tr <- hi_conf_index_list[[num_i]]
  curr_geno_no_tran_index <- geno_no_tran_index_list[[num_i]]
  num_hi_conf_edges        <- length(curr_hi_conf_geno_pheno_tr)


  # Subset delta phenotypes, non-trans indices, and trans indices to just those
  #   with hi conf for geno, pheno, and tr
  hi_conf_delta_pheno <- curr_delta_pheno[curr_hi_conf_geno_pheno_tr]
  hi_conf_geno_no_tran_index <- curr_geno_no_tran_index[curr_geno_no_tran_index %in% curr_hi_conf_geno_pheno_tr]
  hi_conf_geno_yes_tran_index <- curr_hi_conf_geno_pheno_tr[!curr_hi_conf_geno_pheno_tr %in% hi_conf_geno_no_tran_index]
  check_equal(length(hi_conf_delta_pheno), num_hi_conf_edges)
  check_equal(sum(length(hi_conf_geno_yes_tran_index), length(hi_conf_geno_no_tran_index)), num_hi_conf_edges)

  num_trans_edges <- length(hi_conf_geno_yes_tran_index)
  permuted_trans_edge_index <- rep(NA, num_trans_edges)
  set.seed(1) # for reproducability of the sample() function
  # do the permutation part
  # create a random sample of the tr
  null_gamma_values <- rep(NA, perm)
  for (j in 1:perm) {
    permuted_trans_edge_index <-
      sample(curr_hi_conf_geno_pheno_tr,
             size = num_trans_edges,
             replace = FALSE,
             prob = tr$edge.length[curr_hi_conf_geno_pheno_tr] /
               sum(tr$edge.length[curr_hi_conf_geno_pheno_tr]))

    permuted_non_trans_edge_index <- curr_hi_conf_geno_pheno_tr[!curr_hi_conf_geno_pheno_tr %in% permuted_trans_edge_index]
    permuted_trans_delta_pheno <- hi_conf_delta_pheno[permuted_trans_edge_index]
    permuted_non_trans_delta_pheno <- hi_conf_delta_pheno[permuted_non_trans_edge_index]
    check_equal(sum(length(permuted_trans_delta_pheno), length(permuted_non_trans_delta_pheno)), length(hi_conf_delta_pheno))

    median_non_trans_delta_pheno <- median(permuted_non_trans_delta_pheno)
    null_gamma_values[j] <- sum(permuted_trans_delta_pheno > median_non_trans_delta_pheno)
  }
  # Return output --------------------------------------------------------------
  return(null_gamma_values)
}

#' Calculate empirical delta phenotype
#'
#' @description Return Kolmogorov-Smirnov test statistic on the permuted
#'   (empirical) phenotype deltas per tree edge.
#'
#' @param perm Number. Number of permutations.
#' @param permuted_trans_mat Matrix. Nrow = number of permutations. Ncol =
#'   Number of transition edges. Each entry is an edge number.
#' @param hi_conf_edge_nums Vector of edge numbers. Corresponds to edges in tree
#'   with high confidence.
#' @param pheno_by_edges Matrix. Nrow = Nedge(tr). Ncol = 2. Values are the
#'   phenotype reconstruction at each node, as ordered by edges. It's the same
#'   organization as tr$edge.
#'
#' @return empirical_ks_stat
#'
#' @noRd
calculate_empirical_pheno_delta <- function(perm,
                                            permuted_trans_mat,
                                            hi_conf_edge_nums,
                                            pheno_by_edges){
  # Check input ----------------------------------------------------------------
  check_if_permutation_num_valid(perm)
  check_dimensions(permuted_trans_mat,
                   exact_rows = perm,
                   min_rows = perm,
                   exact_cols = NULL,
                   min_cols = 1)
  check_is_number(hi_conf_edge_nums[1])
  check_class(pheno_by_edges, "matrix")
  check_is_number(pheno_by_edges[1, 1])

  # Function -------------------------------------------------------------------
  empirical_ks_stat <- rep(NA, perm)
  for (k in 1:perm) {
    permuted_trans_index <- unique(permuted_trans_mat[k, ])
    permuted_non_trans_index <-
      c(1:length(hi_conf_edge_nums))[!(c(1:length(hi_conf_edge_nums))
                                       %in% unique(permuted_trans_mat[k, ]))]
    empirical_results <- run_ks_test(permuted_trans_index,
                                     permuted_non_trans_index,
                                     pheno_by_edges)
    empirical_ks_stat[k] <- empirical_results$statistic
  }
  # Check and return output ----------------------------------------------------
  return(empirical_ks_stat)
}

#' Calculate significance for continuous data
#'
#' @description Calculate beta(phenotype) and beta(genotype) for each phenotype-
#'   genotype pair. Find their intersections (these are the observed values).
#'   Then permute the genotypes and find the null intersections. Calculate
#'   significance based on the extremity of the observed values as compared to
#'   the null distributions.
#'
#' @param hi_conf_list
#'   \describe {
#'     \item{genotype}{Matrix. Nrow(geno_mat) == number of samples. Each column
#'     is a variant. Matrix is binary.}
#'     \item{genotype_transition}{List of lists.
#'     Length(genotype_transition_list) == num genotypes. For each genotype
#'     there are two lists: $transition & $trans_dir. Length($transition) ==
#'     length($trans_dir) == Nedge(tr). Each of these are either 0, 1, or -1.}
#'     \item{high_conf_ordered_by_edges}{List. Length(list) = ncol(mat) ==
#'     number of genotypes. Each entry is a vector with length == Nedge(tr). All
#'     entries are 0 (low confidence) or 1 (high confidence).}
#' @param permutations Integer. Number of times to run the permutation test.
#' @param tr Phylo.
#' @param pheno_recon_edge_mat Matrix. Nrow = Nedge(tr). Ncol = 2.
#'   Values are the phenotype reconstruction at each node, as ordered by edges.
#'   It's the same organization as tr$edge.
#'
#' @return List of 8:
#'   \describe{
#'     \item{pvals}{Named numeric vector. Length == number of genotypes. Values
#'     between 1 and 0. Names are genotype names.}
#'     \item{null_gammas}{List of numeric vectors. Length of list == number
#'     of genotypes. Each vector has length == number of permutations. Values
#'     are integers}
#'     \item{observed_pheno_trans_delta}{List of numeric vectors. Length of
#'     list == number of genotypes. Vectors are of variable length because
#'     length is the number of transition edges for that particular genotype.
#'     Vectors are numeric.}
#'     \item{observed_pheno_non_trans_delta}{List of numeric vectors. Length of
#'     list == number of genotypes. Vectors are of variable length because
#'     length is the number of non-transition edges for that particular
#'     genotype. Vectors are numeric.}
#'     \item{non_trans_median}{Numberic. Vector. Length = number of genotypes.
#'     Describes median delta phenotype on all genotype non-transition edges.}
#'     \item{num_genotypes}{Integer. The number of genotypes.}
#'     \item{observed_gamma}{Numeric Vector. Length = number of genotypes.
#'     Integers.}
#'   }
#'
#' @noRd
calc_sig <- function(hi_conf_list,
                     permutations,
                     tr,
                     pheno_recon_edge_mat){
  # Check input ----------------------------------------------------------------
  geno_mat <- high_conf_list$genotype
  genotype_transition_list <- high_conf_list$genotype_transition
  genotype_confidence <- high_conf_list$high_conf_ordered_by_edges
  tr_pheno_confidence <- high_conf_list$tr_and_pheno_hi_conf

  check_for_root_and_bootstrap(tr)
  check_if_permutation_num_valid(permutations)
  check_if_binary_matrix(geno_mat)
  check_if_binary_vector(genotype_confidence[[1]])
  check_dimensions(mat = geno_mat,
                   exact_rows = ape::Ntip(tr),
                   min_rows = ape::Ntip(tr),
                   exact_cols = NULL, min_cols = 1)
  check_dimensions(mat = pheno_recon_edge_mat,
                   exact_rows = ape::Nedge(tr),
                   min_rows = ape::Nedge(tr),
                   exact_cols = 2,
                   min_cols = 2)
  check_equal(length(genotype_transition_list[[1]]$transition), ape::Nedge(tr))
  check_equal(length(genotype_transition_list), ncol(geno_mat))
  check_equal(length(genotype_confidence[[1]]), ape::Nedge(tr))

  # Function -------------------------------------------------------------------
  num_genotypes <- ncol(geno_mat)
  num_tree_edge <- nrow(pheno_recon_edge_mat)
  pvals <- non_trans_median <- observed_gamma_value <- pheno_beta <-
    geno_beta <- rep(NA, num_genotypes)
  names(pvals) <- colnames(geno_mat)
  null_gamma_list <- observed_pheno_trans_delta <- geno_non_trans_index_list <-
    high_conf_index_list <- observed_pheno_non_trans_delta <-
    observed_pheno_delta_list <- rep(list(0), num_genotypes)

  for (i in 1:num_genotypes) {
    geno_non_trans_index_list[[i]] <-
      which(genotype_transition_list[[i]]$transition == 0)
    high_conf_index_list[[i]] <-
      which(genotype_confidence[[i]] == 1 & tr_pheno_confidence == 1)
  }

  for (i in 1:num_genotypes) {
    observed_pheno_delta_list[[i]] <-
      calculate_phenotype_change_on_edge(1:num_tree_edge,
                                         pheno_recon_edge_mat)
    observed_pheno_non_trans_delta[[i]] <-
      calculate_phenotype_change_on_edge(geno_non_trans_index_list[[i]],
                                         pheno_recon_edge_mat)
  }

  for (i in 1:num_genotypes) {
    non_trans_median[i] <-
      stats::median(observed_pheno_non_trans_delta[[i]])
  }

  for (i in 1:num_genotypes) {
    pheno_beta[i] <- sum(observed_pheno_delta_list[[i]] >
                           non_trans_median[i] &
                           tr_pheno_confidence[[i]] == 1)
    geno_beta[i] <- sum(genotype_transition_list[[i]]$transition == 1 &
                          genotype_confidence[[i]] == 1)
    observed_gamma_value[i] <-
      sum(observed_pheno_delta_list[[i]] > non_trans_median[i] &
            genotype_transition_list[[i]]$transition == 1 &
            high_conf_index_list == 1)
  }

   # LEFT OFF HERE
  for (i in 1:num_genotypes) {
    # Create permuted transition index matrix and get edge numbers of the high
    #  confidence edges
    null_gamma_distribution <- continuous_permutation(geno_non_trans_index_list,
                                                      tr,
                                                      high_conf_index_list,
                                                      permutations,
                                                      i,
                                                      observed_pheno_delta_list)

    # Nothing should have to change after this point
    # empirical p value calculation here:
    # (1 + more extreme observations) / (1 + permutations)
    pvals[i] <- calculate_permutation_based_p_value(null_gamma_distribution,
                                                    observed_gamma_value[i],
                                                    permutations)

    # Save these for reporting / plots
    null_gamma_list[[i]] <- null_gamma_distribution
  } # end for (i)

  # Return output --------------------------------------------------------------
  results <-
    list("pvals" = pvals,
         "observed_gamma" = observed_gamma_value,
         "null_gamma" = null_gamma_list,
         "observed_pheno_trans_delta" = observed_pheno_trans_delta,
         "observed_pheno_non_trans_delta" = observed_pheno_non_trans_delta,
         "non_trans_median" = non_trans_median,
         "num_genotypes" = num_genotypes)
  return(results)
}

#' Identify the significant loci and apply mulitiple test correction
#'
#' @description Given p-values that have not yet been corrected for multiple
#'   testing, apply a false discovery rate. Using FDR instead of FWER/bonferroni
#'   because these approaches tend to suffer from low power. FDR control
#'   increases power while bounding error.
#'
#' @details Note, implementation of FDR is different from the choice of
#'   bonferroni in the phyC paper.
#'
#' @param hit_values Numeric. Vector of the empirical p-value for each genotype.
#'   Each entry corresponds to respective column in genotype_matrix. Should be
#'   between 0 and 1.
#' @param fdr Numeric. False discovery rate (typically ~10% for discovery work).
#'   Between 0 and 1.
#'
#' @return pvals. List of two data.frames:
#'   \describe{
#'     \item{hit_pvals}{Dataframe. 1 column. Nrow = number of genotypes.
#'     Row.names = genotypes. Column name = "fdr_corrected_pvals". Values
#'     between 1 and 0.}
#'     \item{sig_pvals}{Dataframe. 1 column. Nrow = number of genotypes that
#'     are significant after FDR correction. Column name =
#'     "fdr_corrected_pvals[fdr_corrected_pvals < fdr]". Row.names = genotypes.
#'     Nrow = is variable-- could be between 0 and max number of genotypes. It
#'     will only have rows if the corrected p-value is less than the fdr value.}
#'   }
#'
#' @noRd
get_sig_hit_and_mult_test_corr <- function(hit_values, fdr){
  # Check inputs ---------------------------------------------------------------
  check_num_between_0_and_1(fdr)

  # Function -------------------------------------------------------------------
  fdr_corrected_pvals <- stats::p.adjust(hit_values, method = "fdr")
  sig_pvals <- as.data.frame(fdr_corrected_pvals[fdr_corrected_pvals < fdr])
  sig_pvals <- as.data.frame(sig_pvals)
  row.names(sig_pvals) <- names(hit_values)[fdr_corrected_pvals < fdr]
  fdr_corrected_pvals <- as.data.frame(fdr_corrected_pvals)
  row.names(fdr_corrected_pvals) <- names(hit_values)

  colnames(fdr_corrected_pvals) <- "fdr_corrected_pvals"
  colnames(sig_pvals) <- "fdr_corrected_pvals"

  # Check and return output ----------------------------------------------------
  results <- list("hit_pvals" = fdr_corrected_pvals,
                  "sig_pvals" = sig_pvals)
  return(results)
}
