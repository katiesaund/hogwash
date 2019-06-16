calculate_permutation_based_p_value <- function(empirical_statistic, observed_statistic, num_perm){
  # Function description -------------------------------------------------------
  # Given all of the empirical statistics derived from permutations, count how
  # many of the empirical/permuted test statistics are greater than or equal to
  # the observed/real test statistic.
  # Adding one in the numerator and denominator accounts for the observed value.
  #
  # Inputs:
  # empirical_statistic. Numeric vector. Length = num_perm.
  # observed_statistic. Number. Length = 1.
  # num_perm. Number.
  #
  # Outputs:
  # pval = Number. Length = 1. Value between 0 and 1.
  #
  # Check input ----------------------------------------------------------------
  check_equal(length(empirical_statistic), num_perm)
  check_is_number(observed_statistic)
  check_is_number(num_perm)
  check_if_permutation_num_valid(num_perm)
  check_is_number(empirical_statistic[1])

  # Function -------------------------------------------------------------------
  pval <- (sum(1 + sum(empirical_statistic >=  observed_statistic)) / (num_perm + 1))

  # Check and return output ----------------------------------------------------
  check_num_between_0_and_1(pval)
  return(pval)
}

count_hits_on_edges <- function(genotype_transition_edges, phenotype_reconstruction, high_confidence_edges, tr){
  # Function description -------------------------------------------------------
  # For each genotype, return a count of the number of edges for which the
  # genotype is a transition and the phenotype is present (when using
  # reconstruction)/phenotype is a transition (when using transition). Also,
  # return a count of the number of edges for which the just the genotype is a
  # transition (phenotype is either a 0 reconstruction or not a transition).
  # TODO: This approach assumes that the genotype_transition_edges being fed into it are all 0 > 1 for original phyC. Not sure if there is an assumption in place for the overlap test.
  # Inputs:
  # genotype_transition_edges. List of numeric vectors. Each vector corresponds with one genotype from the geno_mat. The numbers in each vector correspond to an edge on the tree.
  # phenotype_reconstruction. Numeric vector. Each entry corresponds to an edge.
  # high_confidence_edges. List of vectors. Each vector corresponds to the confidence of 1 genotype from the geno_mat. Each entry in the vector corresponds to an edge on the tree.
  #                        1 means high confidence in the node, 0 means low confidence.
  #
  # Outputs:
  # "both_present" = both_present. Length = number of genotypes.
  # "only_geno_present" = only_geno_present. Length = number of genotypes.
  #
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_is_number(genotype_transition_edges[[1]][1])
  check_is_number(high_confidence_edges[[1]][1])
  check_is_number(phenotype_reconstruction[1])
  check_equal(length(genotype_transition_edges[[1]]), ape::Nedge(tr))
  check_equal(length(phenotype_reconstruction), ape::Nedge(tr))
  check_equal(length(high_confidence_edges[[1]]), ape::Nedge(tr))

  # Function -------------------------------------------------------------------
  both_present <- sapply(1:length(high_confidence_edges), function(x) {
    sum(phenotype_reconstruction[as.logical(high_confidence_edges[[x]])] + genotype_transition_edges[[x]][as.logical(high_confidence_edges[[x]])] == 2)
  })

  only_geno_present <- sapply(1:length(high_confidence_edges), function(x) {
    sum(genotype_transition_edges[[x]][as.logical(high_confidence_edges[[x]])]) - both_present[x]
  })

  # Check and return output ----------------------------------------------------
  hit_counts <- list("both_present" = both_present, "only_geno_present" = only_geno_present)
  return(hit_counts)
} #end count_hits_on_edges()

discrete_calculate_pvals <- function(genotype_transition_edges, phenotype_reconstruction, tr, mat, permutations, fdr, high_confidence_edges){
  # Function description -------------------------------------------------------
  # discrete_calculate_pvals is the "meat" of the discrete phyC algorithm.
  # discrete_calculate_pvals returns the empirical p-value for each observed genotype.
  # Algorithm overview:
  # 1. Subset tree edges to those with high confidence (as determined by phenotype & genotype reconstructions as well as tree bootstrap values).
  # 2. For each genotype from the genotype_matrix:
  #    i.   Exclude any edges on which the genotype is not present
  #    ii.  Create a matrix where each row is a random sampling of the high confidence edges of the tree where probability of choosing edges is proportional to length of edge. Number of
  #         edges selected for the permuted data set is the number of times the empirical genotype appears.
  #    iii. Calculate when the randomly permuted genotype overlap with the high confidence phenotype (these values create the null distribution for the permuted genotypes)
  #    iv.  Calculate empirical p-values.
  #
  # Inputs:
  # genotype_transition_edges. List of vectors. Length(genotype_transition_list) == number of genotypes.
  #                           Length(each vector) == Nedge(tr).
  #                           Each of these are either 0, 1.
  # phenotype_reconstruction. Vector. Length() == Nedge(tr).
  #   Note: Phenotype reconstruction could instead by phenotype transition vector.
  # tr. Phylo.   Rooted phylogenetic tree.
  # mat. Numeric matrix.  The genotype matrix. Nrow(mat) == number of isolates == Ntip(tr). Each column is a variant. Matrix is binary (0 or 1).
  # permutations. Integer.  The number of times to run the permutation test.
  # fdr. Number.  Significance threshold.
  # high_confidence_edges. List. Length(list) = ncol(mat) == number of genotypes.
  #                      Each entry is a vector with length == Nedge(tr). All
  #                      entries are 0 (low confidence) or 1 (high confidence).
  #
  # Outputs:
  # "hit_pvals" = hit_pvals.
  # "permuted_count" = record_of_redistributed_both_present.
  # "observed_overlap" = both_present.
  #
  # Check input ----------------------------------------------------------------
  check_equal(ncol(mat), length(genotype_transition_edges))
  check_equal(length(genotype_transition_edges[[1]]), ape::Nedge(tr))
  check_equal(length(phenotype_reconstruction), ape::Nedge(tr))
  check_if_permutation_num_valid(permutations)
  check_for_root_and_bootstrap(tr)
  check_num_between_0_and_1(fdr)
  check_equal(ncol(mat), length(high_confidence_edges))
  check_equal(length(high_confidence_edges[[1]]), ape::Nedge(tr))

  # Function -------------------------------------------------------------------
  # Calculate observed values
  # (Convergence of 0 -> 1 on phenotype present edges (variable: both present)
  # versus phenotype absent edges (only_geno_present)).
  # Phenotype present means 1 in reconstruction test or phenotype is a transition (1) in the overlap/transition test. Equivalent to Farhat et al. "Resistant" branches.
  # Phenotype absent (only_geno_present) is equivalen to Farhat et al. "Sensitive" branches.
  observed_result   <- count_hits_on_edges(genotype_transition_edges, phenotype_reconstruction, high_confidence_edges, tr)
  both_present      <- observed_result$both_present
  only_geno_present <- observed_result$only_geno_present

  # initialize some values
  num_genotypes                         <- ncol(mat)
  num_edges_with_geno_trans             <- both_present + only_geno_present
  hit_pvals                             <- rep(NA, num_genotypes)
  num_hi_conf_edges                     <- sapply(high_confidence_edges, function(x) sum(x))
  list_of_all_edges                     <- c(1:ape::Nedge(tr))
  hi_conf_edges                        <- list(rep(0, num_genotypes))

  # 1. Subset tr edges to those with high confidence (as determined by phenotype & genotype reconstructions as well as tr bootstrap values).
  for (i in 1:num_genotypes) {
    hi_conf_edges[[i]] <- list_of_all_edges[as.logical(high_confidence_edges[[i]])]
  }

  # 2. For each genotype from the genotype_matrix:
  record_of_redistributed_both_present <- rep(list(0), num_genotypes)
  for (i in 1:num_genotypes) { # looping over each genotype in the genotype_matrix
    if (num_edges_with_geno_trans[i] > num_hi_conf_edges[i]) {
      stop("Too many hits on the branches")
    }
    if ((both_present[[i]] + only_geno_present[[i]]) < 2) {
      # If there are 1 or 0 high confidence edges with the genotype present then the p-value should be reported as 1.0;
      # both_present and only_geno_present are made up of only high confidence branches as defined in count_hits_on_edges
      # which isn't quite true but will indicate that we cannot calculate a p-value because we cannot detect any
      # convergence with one or fewer affected branches. And that we should skip to the next genotype to run the permutation test.
      hit_pvals[i] <- 1.0
      # This should never get triggered because we should be filtering the genotype before this step, but it's being kept in as a fail safe.
    } else {

      permuted_geno_trans_edges <- discrete_permutation(tr, permutations, num_edges_with_geno_trans, num_hi_conf_edges, ape::Nedge(tr), hi_conf_edges, i)
      # Now calculate both_present and only_geno_present with the permuted data in the same fashion as with the observed data
      empirical_both_present <- count_empirical_both_present(permuted_geno_trans_edges, phenotype_reconstruction, high_confidence_edges, i)
      empirical_only_geno_present <- count_empirical_only_geno_present(permuted_geno_trans_edges, empirical_both_present)

      # Note on nomeclature from the original PhyC paper supplement page 7: https://media.nature.com/original/nature-assets/ng/journal/v45/n10/extref/ng.2747-S1.pdf
      # X -> G on R is the same as sum(empirical_both_present >= both_present[i])
      # X -> G on S is the same as sum(empirical_only_geno_present <= only_geno_present[i])
      # Note: sum(empirical_both_present >= both_present[i]) always equals sum(empirical_only_geno_present <= only_geno_present[i])
      #       so we only need to include one in the p-value calculation.
      #       In an older version I used the following block of code, but it's equivalent to the single line below using the calculate_permutation_based_p_value() function.

              # count only times when permuted (empirical) overlap of genotype and phenotype is more common than obsered (both_present[i])
              # and when the permuted (empirical) genotype does not overlap with phenotype is less common than observed (only_geno_present[i])
              # permutation_as_or_more_extreme <- sum((empirical_both_present >= both_present[i]) * (empirical_only_geno_present <= only_geno_present[i]))
              # temp_pval <- ((permutation_as_or_more_extreme + 1)/(permutations + 1))

      pval <- calculate_permutation_based_p_value(empirical_both_present, both_present[i], permutations)
      hit_pvals[i] <- format(round(pval, 20), nsmall = 20)
      record_of_redistributed_both_present[[i]] <- empirical_both_present
    }
  }
  names(hit_pvals) <- colnames(mat)

  # Return output --------------------------------------------------------------
  results <- list("hit_pvals" = hit_pvals, "permuted_count" = record_of_redistributed_both_present, "observed_overlap" = both_present)
  return(results)
} # end discrete_calculate_pvals

discrete_permutation <- function(tr, num_perm, number_edges_with_geno_trans, number_hi_conf_edges, number_edges, high_conf_edges, index){
  # Function description -------------------------------------------------------
  # Perform a permutation of the edges selected to be genotype transition edges.
  # permuted_geno_trans_mat is a matrix where each row corresponds to a permuted set of edges and each column corresponds to one high confidence edge.
  #   To get the edges selected in a permutation, take that row from permuted_geno_trans_mat.
  #   We could use this format only for later analyses, but a more convenient format is to conver to redistributed_hits.
  # redistributed_hits is a matrix with nrow = number of perms (same as permuted_geno_trans_mat, but each row now corresponds to a genotype)
  #
  # Inputs:
  # tr. Phylo.
  # num_perm. Number of permutations.
  # number_edges_with_geno_trans. Numeric vector. Max value should be num edge in tr.
  # number_hi_conf_edges. Numeric vector. Max value should be num edge in tr.
  # number_edges. Nedge(tr).
  # high_conf_edges. Numeric vector of edges that high confidence. Max value should be num edge in tr.
  # index. Integer. "i" from loop. Max value is number of genotypes. Min value == 1.
  #
  # Outputs:
  #
  # Check input ----------------------------------------------------------------
  check_if_permutation_num_valid(num_perm)
  check_for_root_and_bootstrap(tr)
  check_is_number(number_edges_with_geno_trans[index])
  check_is_number(number_hi_conf_edges[index])
  check_is_number(number_edges)
  check_is_number(index)
  check_is_number(high_conf_edges[[index]][1])
  check_equal(number_edges, ape::Nedge(tr))
  if (index < 1) {
    stop("loop index must be positive")
  }
  if (max(high_conf_edges[[index]]) > ape::Nedge(tr) | number_hi_conf_edges[index] > ape::Nedge(tr) | number_edges_with_geno_trans[index] > ape::Nedge(tr)) {
    stop("max value should be Nedge(tr")
  }


  # Function -------------------------------------------------------------------
  permuted_geno_trans_mat <- matrix(nrow = num_perm, ncol = number_edges_with_geno_trans[index])
  redistributed_hits   <- matrix(0, nrow = num_perm, ncol = number_edges)
  set.seed(1)
  for (j in 1:num_perm) {
    permuted_geno_trans_mat[j, ] <- sample(1:number_hi_conf_edges[index],
                                           size = number_edges_with_geno_trans[index],
                                           replace = TRUE,
                                           prob = tr$edge.length[high_conf_edges[[index]]]/sum(tr$edge.length[high_conf_edges[[index]]]))
  }

  #this if statement deals with situation when number_edges_with_geno_trans = 1; which should never happen now that we prefilter genotypes.
  if (nrow(permuted_geno_trans_mat) != num_perm) {
    permuted_geno_trans_mat <- t(permuted_geno_trans_mat)
  }

  for (m in 1:nrow(permuted_geno_trans_mat)) {# or 1:nrow(permuted_geno_trans_mat) is the same as 1:permutations
    redistributed_hits[m, ][permuted_geno_trans_mat[m, ]] <- 1
  }
  # Check and return output ----------------------------------------------------
  check_dimensions(redistributed_hits, exact_rows = num_perm, min_rows = num_perm, exact_cols = ape::Nedge(tr), min_cols = ape::Nedge(tr))
  return(redistributed_hits)
} # end discrete_permutation()

count_empirical_both_present <- function(permuted_mat, pheno_vec, hi_conf_edge, index){
  # Function description -------------------------------------------------------
  # Find number of edges on the tree where the permuted genotype is a transition
  # AND the phenotype is also a transition (or a phenotype is present in the
  # reconstruction; it depends on the type of test being run).
  #
  # Inputs:
  # permuted_mat. Matrix. Nrow = number of permutations. Ncol = number of tr edges. Either 0 or 1.
  # pheno_vec. Numeric vector.
  # hi_conf_edge.
  # index. Number. The "i" of the loop this function is run within.
  #
  # Outputs:
  # result. Numeric vector. Length = num perm.
  #
  # Check input ----------------------------------------------------------------
  check_is_number(index)
  check_if_binary_matrix(permuted_mat)
  check_equal(length(pheno_vec),ncol(permuted_mat))
  check_equal(length(pheno_vec), length(as.logical(hi_conf_edge[[index]])))

  # Function -------------------------------------------------------------------
  result <- sapply(1:nrow(permuted_mat), function(x) {
    sum((pheno_vec + permuted_mat[x, ]) == 2 & hi_conf_edge[[index]] == 1)
  })

  # Check and return output ----------------------------------------------------
  if (nrow(permuted_mat) != length(result)) {
    stop("result dimension mismatch")
  }
  return(result)
} # end count_empirical_both_present()

count_empirical_only_geno_present <- function(permuted_mat, emp_both_present){
  # Function description -------------------------------------------------------
  # Find number of edges on the tree where the permuted genotype is a transition
  # but the phenotype is not (either a transition phenotype edge or phenotype
  # reconstruction present -- it depends on the type of test being run).
  #
  # Inputs:
  # permuted_mat. Matrix. Nrow = number of permutations. Ncol = number of edges. Either 0 or 1.
  # emp_both_present. Numeric vector.
  #
  # Outputs:
  # result. Numeric vector. Length = num perm.
  #
  # Check input ----------------------------------------------------------------
  check_if_binary_matrix(permuted_mat)
  check_equal(length(emp_both_present), nrow(permuted_mat))

  # Function -------------------------------------------------------------------
  result <- sapply(1:nrow(permuted_mat), function(x) {
    sum(permuted_mat[x, ]) - emp_both_present[x]
  })
  # Check & return output --------------------------------------------------------------
  if (nrow(permuted_mat) != length(result)) {
    stop("result dimension mismatch")
  }

  return(result)
} # end count_empirical_only_geno_present()

