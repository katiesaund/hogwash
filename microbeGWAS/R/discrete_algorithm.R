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
  if (length(empirical_statistic) != num_perm){
    stop("Number of empirical test statistics should be equal to the number of permutations")
  }
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

  if(length(genotype_transition_edges[[1]]) != Nedge(tr)){
    stop("Genotype transition edge vector must have an entry for each tree edge.")
  }
  if(length(phenotype_reconstruction) != Nedge(tr)){
    stop("Phenotype transition or reconstruction edge vector must have an entry for each tree edge.")
  }
  if(length(high_confidence_edges[[1]]) != Nedge(tr)){
    stop("Reconstruction confidence and must have an entry for each tree edge.")
  }
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
  if (ncol(mat) != length(genotype_transition_edges)){
    stop("Genotype transition edges should have a vector for each genotype")
  }
  if (length(genotype_transition_edges[[1]]) != Nedge(tr)){
    stop("Genotype transition edges should be made of vectors of length == Nedge(tree)")
  }
  if (length(phenotype_reconstruction) != Nedge(tr)){
    stop("Phenotype reconstruction or transition vector should have length == Nedge(tree)")
  }
  check_if_permutation_num_valid(permutations)
  check_for_root_and_bootstrap(tr)
  check_num_between_0_and_1(fdr)
  if (ncol(mat) != length(high_confidence_edges)){
    stop("Confidence list should have a vector for each genotype.")
  }
  if (length(high_confidence_edges[[1]]) != Nedge(tr)){
    stop("Confidence list should be made of vectors of length == Nedge(tree)")
  }


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
  distribution_of_hits_from_permutation <- rep(0, num_genotypes)
  hit_pvals                             <- rep(NA, num_genotypes)
  num_hi_conf_edges                     <- sapply(high_confidence_edges, function(x) sum(x))
  list_of_all_edges                     <- c(1:Nedge(tr))
  hi_conf_edges                        <- list(rep(0, num_genotypes))

  # 1. Subset tr edges to those with high confidence (as determined by phenotype & genotype reconstructions as well as tr bootstrap values).
  for (i in 1:num_genotypes){
    hi_conf_edges[[i]] <- list_of_all_edges[as.logical(high_confidence_edges[[i]])]
  }

  # 2. For each genotype from the genotype_matrix:
  record_of_redistributed_both_present <- rep(list(0), num_genotypes)
  for (i in 1:num_genotypes){ # looping over each genotype in the genotype_matrix
    if (num_edges_with_geno_trans[i] > num_hi_conf_edges[i]){
      stop("Too many hits on the branches")
    }
    if((both_present[[i]] + only_geno_present[[i]]) < 2){
      # If there are 1 or 0 high confidence edges with the genotype present then the p-value should be reported as 1.0;
      # both_present and only_geno_present are made up of only high confidence branches as defined in count_hits_on_edges
      # which isn't quite true but will indicate that we cannot calculate a p-value because we cannot detect any
      # convergence with one or fewer affected branches. And that we should skip to the next genotype to run the permutation test.
      hit_pvals[i] <- 1.0
      # This should never get triggered because we should be filtering the genotype before this step, but it's being kept in as a fail safe.
    } else {
      # initialize some variables for this loop
      permuted_geno_trans_mat <- matrix(nrow = permutations, ncol = num_edges_with_geno_trans[i])
      redistributed_hits   <- matrix(0, nrow = permutations, ncol = num_genotypes)

      # create a random sample of the tr
      print("1:num_hi_conf_edges[i]")
      print(1:num_hi_conf_edges[i])

      print("num_edges_with_geno_trans[i]")
      print(num_edges_with_geno_trans[i])

      print("hi_conf_edges[[i]]]")
      print(hi_conf_edges[[i]])

      set.seed(1)
      for(j in 1:permutations){
        permuted_geno_trans_mat[j, ] <- sample(1:num_hi_conf_edges[i],
                                            size = num_edges_with_geno_trans[i],
                                            replace = TRUE,
                                            prob = tr$edge.length[hi_conf_edges[[i]]]/sum(tr$edge.length[hi_conf_edges[[i]]]))
      }

      #this if statement deals with situation when num_edges_with_geno_trans = 1; which should never happen now that we prefilter genotypes.
      if(nrow(permuted_geno_trans_mat) != permutations){
        permuted_geno_trans_mat <- t(permuted_geno_trans_mat)
      }

      for (m in 1:nrow(permuted_geno_trans_mat)){ # or 1:nrow(permuted_geno_trans_mat is the same as 1:permutations
        redistributed_hits[m, ][permuted_geno_trans_mat[m, ]] <- 1
      }

      print("red")
      print(redistributed_hits)
      print("perm")
      print(permuted_geno_trans_mat)

      # right now redistributed hits is simply a matrix marking which edges we hit during
      # the permutation step above- because the point of the permutation is to pretend as if
      # we're assigning the presence of the genotype to random edges in the tree.


      # Now calculate both_present and only_geno_present with the permuted data in the same fashion as with the observed data
      empirical_both_present <- sapply(1:nrow(redistributed_hits), function(x) {
        sum(phenotype_reconstruction[as.logical(high_confidence_edges[[i]])] + redistributed_hits[x, ] == 2)
      })

      empirical_only_geno_present <- sapply(1:nrow(redistributed_hits), function(x) {
        sum(redistributed_hits[x, ]) - empirical_both_present[x]
      })

      print("empirical_both_present")
      print(empirical_both_present)
      print("both_present[i]")
      print(both_present[i])
      print("empirical_only_geno_present")
      print(empirical_only_geno_present)
      print("only_geno_present[i]")
      print(only_geno_present[i])


      # I'm pretty sure that redistirbuted_both_present is the same thing as the empirical_both_present so I'm going to remove redistributed_both_present
      # x_on_r <- sum(empirical_both_present >= both_present[i])
      # y_on_s <- sum(empirical_only_geno_present <= only_geno_present[i])

      # count only times when permuted (empirical) overlap of genotype and phenotype is more common than obsered (both_present[i])
      # and when the permuted (empirical) genotype does not overlap with phenotype is less common than observed (only_geno_present[i])
      new_counter <- sum((empirical_both_present >= both_present[i]) * (empirical_only_geno_present <= only_geno_present[i]))
      temp_pval <- ((new_counter + 1)/(permutations + 1))


      print("sort(empirical_both_present, decreasing = FALSE)[(fdr * permutations)]")
      print(sort(empirical_both_present, decreasing = FALSE)[(fdr * permutations)])
      print(sort(empirical_both_present, decreasing = FALSE))
      print("fdr * permutations")
      print(fdr * permutations)

      print("both_present[i")
      print(both_present[i])

      if (sort(empirical_both_present, decreasing = FALSE)[(fdr * permutations)] == 0 & both_present[i] == 0){ # i have no idea why this line exists
        pval <- 1
      } else if (temp_pval == 0 | temp_pval == 1){
        pval <- 2/(permutations + 1)
      } else if (temp_pval > 0.5){
        pval <- ((1 - temp_pval)  * 2)
      } else if (temp_pval <= 0.5){
        pval <- (temp_pval * 2)
      }

      record_of_redistributed_both_present[[i]] <- empirical_both_present
      hit_pvals[i] <- format(round(pval, 20), nsmall = 20)
    }
  }
  names(hit_pvals) <- colnames(mat)

  # Return output --------------------------------------------------------------
  results <- list("hit_pvals" = hit_pvals, "permuted_count" = record_of_redistributed_both_present, "observed_overlap" = both_present)

  return(results)
} # end discrete_calculate_pvals


# 2019-05-15 why the heck are the pvals so weird for discrete?

# weird_pval <- function(temp_pval, permutations){
#
#   pval <- NA
#   if (temp_pval == 0 | temp_pval == 1){
#     pval <- 2/(permutations + 1)
#   } else if (temp_pval > 0.5){
#     pval <- ((1 - temp_pval)  * 2)
#   } else if (temp_pval <= 0.5){
#     pval <- (temp_pval * 2)
#   }
#   return(pval)
# }
#
#
# pval_seq <- seq(from = 0, to  = 1, by = 0.01)
#
# weird_value <- rep(NA, length(pval_seq))
# for (i in 1:length(pval_seq)){
#   weird_value[i] <- weird_pval(pval_seq[i], 10000)
# }
#
# plot(weird_value, pval_seq, xlab = "output p-value", ylab = "input p-value")
#
# pval <- rep(NA, 1000)
# for (j in 1:1000){
#   pval[j] <- (j + 1)/(1000 + 1)
# }
# plot(pval)
#
# plot(weird_value)

