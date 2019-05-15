library(microbeGWAS)
context("Discrete algorithm") -------------------------------------------------#

# test calculate_permutation_based_p_value
test_that("calculate_permutation_based_p_value returns a significant p-value when statistics are much lower than observed", {
  nperm <- 1000
  permuted_tests <- rep(0.1, nperm)
  real_test <- 0.9
  alpha <- 0.01
  expect_true(calculate_permutation_based_p_value(permuted_tests, real_test, nperm) < alpha)
})

test_that("calculate_permutation_based_p_value returns a non-significant p-value when statistics are cenetered around observed", {
  nperm <- 1000
  permuted_tests <- rnorm(n = nperm, mean = 0)
  real_test <- 0
  alpha <- 0.01
  expect_true(calculate_permutation_based_p_value(permuted_tests, real_test, nperm) > alpha)
})

# test count_hits_on_edges
test_that("count_hits_on_edges returns 3 edges shared and 7 edges only with genotype given this test data", {
  num_samples <- 6
  temp_tree <- rtree(num_samples)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  num_edge <- Nedge(temp_tree)
  temp_geno_recon <- temp_hi_conf_edges <- rep(list(0), num_samples)
  for (k in 1:num_samples){
    temp_geno_recon[[k]] <- temp_hi_conf_edges[[k]] <- rep(1, num_edge)
  }
  num_pheno_and_geno_present <- 3
  num_just_geno_present <- num_edge - num_pheno_and_geno_present
  temp_pheno_recon <- c(rep(1, num_pheno_and_geno_present), rep(0, num_just_geno_present))

  results <- count_hits_on_edges(temp_geno_recon, temp_pheno_recon, temp_hi_conf_edges, temp_tree)
  expect_equal(results$both_present[1], num_pheno_and_geno_present)
  expect_equal(results$only_geno_present[1], num_just_geno_present)
})

test_that("count_hits_on_edges returns 0 edges shared and 0 edges only with genotype given this all absent genotype", {
  num_samples <- 6
  temp_tree <- rtree(num_samples)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  num_edge <- Nedge(temp_tree)
  temp_geno_recon <- temp_hi_conf_edges <- rep(list(0), num_samples)
  for (k in 1:num_samples){
    temp_geno_recon[[k]] <- temp_hi_conf_edges[[k]] <- rep(0, num_edge)
  }
  temp_pheno_recon <- rep(1, num_edge)
  results <- count_hits_on_edges(temp_geno_recon, temp_pheno_recon, temp_hi_conf_edges, temp_tree)
  expect_equal(results$both_present[1], 0)
  expect_equal(results$only_geno_present[1], 0)
})


# test discrete_calculate_pvals
test_that("discrete_calculate_pvals returns X given Y", {
  num_samples <- 6
  set.seed(1)
  temp_tree <- rtree(num_samples)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  num_edge <- Nedge(temp_tree)
  temp_geno_trans <- temp_hi_conf_edges <- rep(list(0), num_samples)
  for (k in 1:num_samples){
    temp_geno_trans[[k]] <- c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0)
    temp_hi_conf_edges[[k]] <- rep(1, num_edge)
  }
  temp_pheno_trans <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  temp_geno <- matrix(1, ncol = num_samples, nrow = Ntip(temp_tree)) # doesn't match recon or transition, just made up for now.
  temp_perm <- 1000
  temp_fdr <- 0.25
  disc_trans_results <- discrete_calculate_pvals(temp_geno_trans, temp_pheno_trans, temp_tree, temp_geno, temp_perm, temp_fdr, temp_hi_conf_edges)

  expect_equal(round(as.numeric(disc_trans_results$hit_pvals[1]), 3), 0.763)
  expect_equal(disc_trans_results$observed_overlap[1], 1)
  expect_equal(length(disc_trans_results$permuted_count[[1]]), temp_perm)
})


# calculate_hit_pvals_corrected <- function(hit_counts, phenotype_reconstruction, tr, mat, permutations, alpha, high_confidence_edges){
#   # calculate_genotype_pvals is the "meat" of the phyC algorithm.
#   # calculate_genotype_pvals returns the empirical p-value for each observed genotype.
#   # Algorithm overview:
#   # 1. Subset tree edges to those with high confidence (as determined by phenotype & genotype reconstructions as well as tree bootstrap values).
#   # 2. For each genotype from the genotype_matrix:
#   #    i.   Exclude any edges on which the genotype is not present
#   #    ii.  Create a matrix where each row is a random sampling of the high confidence edges of the tree where probability of choosing edges is proportional to length of edge. Number of
#   #         edges selected for the permuted data set is the number of times the empirical genotype appears.
#   #    iii. Calculate when the randomly permuted genotype overlap with the high confidence phenotype (these values create the null distribution for the permuted genotypes)
#   #    iv.  Calculate empirical p-values.
#   #    v.   Plot observed values on the null distribution.
#
#   # Inputs:
#   # hit_counts.               List.    2 part list containing results of count_hits_on_edges
#   # phenotype_reconstruction. Numeric. 0/1. Length of Nedge(tree)
#   # tree:                     Phylo.   Rooted phylogenetic tree.
#   # mat:                      Matrix.  The genotype matrix. All entries are either 0 or 1. No NAs.
#   # permutations:             Number.  The number of permutations.
#   # alpha.                    Number.  Significance threshold.
#
#   # Function description -------------------------------------------------------
#   # TODO
#   # Compute ancestral reconstruction from a continuous or discrete input.
#   #
#   # Inputs:
#   # Varname. Var class. Description.
#   #
#   # Outputs:
#   # "anc_rec" = ML_anc_rec. Vector. Description.
#   #
#   # Check input ----------------------------------------------------------------
#
#   # Function -------------------------------------------------------------------
#
#   # Return output --------------------------------------------------------------
#
#   # read in variables
#   both_present      <- hit_counts$both_present
#   only_geno_present <- hit_counts$only_geno_present
#
#   # initialize some values
#   num_sample                            <- both_present + only_geno_present
#   distribution_of_hits_from_permutation <- rep(0, length(both_present))
#   hit_pvals                             <- rep(NA, ncol(mat))
#   num_branch                            <- sapply(high_confidence_edges, function(x) sum(x))
#   all_edges                             <- c(1:Nedge(tr))
#   which_branches                        <- list(rep(0, length(high_confidence_edges)))
#
#   # 1. Subset tr edges to those with high confidence (as determined by phenotype & genotype reconstructions as well as tr bootstrap values).
#   for (i in 1:length(high_confidence_edges)){
#     if (length(all_edges) == length(high_confidence_edges[[i]])){
#       which_branches[[i]] <- all_edges[as.logical(high_confidence_edges[[i]])]
#     } else {
#       stop("which_branches incorrect")
#     }
#   }
#
#   # 2. For each genotype from the genotype_matrix:
#   record_of_redistributed_both_present <- rep(list(0), length(both_present))
#   for (i in 1:length(both_present)){ # looping over each genotype in the genotype_matrix
#     if (num_sample[i] > num_branch[i]){
#       stop("Too many hits on the branches")
#     }
#     if((both_present[[i]] + only_geno_present[[i]]) < 2){
#       # If there are 1 or 0 high confidence edges with the genotype present then the p-value should be reported as 1.0;
#       # both present and just hit present are made up of only high confidence branches as defined in count_hits_on_edges
#       # which isn't quite true but will indicate that we cannot calculate a p-value because we cannot detect any
#       # convergence with one or fewer affected branches. And that we should skip to the next genotype to run the permutation test.
#       hit_pvals[i] <- 1.0
#     } else {
#       # initialize some counters/variables for this loop
#       counter              <- 0
#       resampling_record    <- numeric(0)
#       all_sampled_branches <- matrix(nrow = permutations, ncol = num_sample[i])
#       redistributed_hits   <- matrix(0, nrow = permutations, ncol = length(which_branches[[i]]))
#       # create a random sample of the tr
#       set.seed(1)
#       for(j in 1:permutations){
#         curr_num_branch <- num_branch[i]
#         all_sampled_branches[j, ] <- sample(1:curr_num_branch,
#                                             size = num_sample[i],
#                                             replace = TRUE,
#                                             prob = tr$edge.length[which_branches[[i]]]/sum(tr$edge.length[which_branches[[i]]]))
#       }
#
#       #this if statement deals with situation when num_sampled = 1
#       if(nrow(all_sampled_branches) != permutations){
#         all_sampled_branches <- t(all_sampled_branches)
#       }
#
#       for (m in 1:nrow(all_sampled_branches)){
#         redistributed_hits[m, ][all_sampled_branches[m, ]] <- 1
#       }
#       # right now redistributed hits is simply a matrix marking which edges we hit during
#       # the permutation step above- because the point of the permutation is to pretend as if
#       # we're assigning the presence of the genotype to random edges in the tree.
#
#       empirical_both_present <- sapply(1:nrow(redistributed_hits), function(x) {
#         sum(phenotype_reconstruction[as.logical(high_confidence_edges[[i]])] + redistributed_hits[x, ] == 2)
#       })
#
#       empirical_only_geno_present <- sapply(1:nrow(redistributed_hits), function(x) {
#         sum(redistributed_hits[x, ]) - empirical_both_present[x]
#       })
#
#       # now beginning calculating when the randomly permuted "genotypes" overlap with the hi confidence phenotype
#       temp <- matrix(0, nrow = permutations, ncol = length(which_branches[[i]])) #used to be length(Nedge(tree))
#       redistributed_both_present <- rep(0, permutations)
#       two <- rep(2, length(which_branches[[i]]))
#
#       for (p in 1:permutations){
#         if (length(redistributed_hits[p, ]) == length(phenotype_reconstruction[as.logical(high_confidence_edges[[i]])])){
#           temp[p, ] <- redistributed_hits[p, ] + phenotype_reconstruction[as.logical(high_confidence_edges[[i]])]
#         } else {
#           stop("Dimension mismatch redistributed_hits and phenotype_reconstruction.")
#         }
#       }
#
#       # temp ranges from 0 to 2.
#       # A value of 0 means that the edge (as indicated by column) was not selected in that permutation (as indicated by row) & it's a low confidence pheno edge
#       # A value of 1 means that either the edge was not selected in that permutation but is a low confidence edge or vice versa
#       # A value of 2 means that the edge was selected in that permutation and it's a high confidence edge.
#       for (q in 1:permutations){
#         redistributed_both_present[q] <- sum(temp[q, ] == two)
#       }
#
#       # x_on_r <- sum(empirical_both_present >= both_present[i])
#       # y_on_s <- sum(empirical_only_geno_present <= only_geno_present[i])
#
#       # count only times when permuted (empirical) overlap of genotype and phenotype is more common than obsered (both_present[i])
#       # and when the permuted (empirical) genotype does not overlap with phenotype is less common than observed (only_geno_present[i])
#       new_counter <- sum((empirical_both_present >= both_present[i]) * (empirical_only_geno_present <= only_geno_present[i]))
#       temp_pval <- ((new_counter + 1)/(permutations + 1))
#       # Katie TODO 2019-02-26: save the new counter for each permutation, then plot that with the real as a red bar, add p val on top.
#
#       if (sort(redistributed_both_present, decreasing = FALSE)[(alpha * permutations)] == 0 & both_present[i] == 0){
#         pval <- 1
#       } else if (temp_pval == 0 | temp_pval == 1){
#         pval <- 2/(permutations + 1)
#       } else if (temp_pval > 0.5){
#         pval <- ((1 - temp_pval)  * 2)
#       } else if (temp_pval <= 0.5){
#         pval <- (temp_pval * 2)
#       }
#
#       record_of_redistributed_both_present[[i]] <- redistributed_both_present
#       hit_pvals[i] <- format(round(pval, 20), nsmall = 20)
#     }
#   }
#   names(hit_pvals) <- colnames(mat)
#   results <- list("hit_pvals" = hit_pvals, "permuted_count" = record_of_redistributed_both_present, "observed_overlap" = both_present)
#   # return(hit_pvals)
#   return(results)
# } # end calculate_hit_pvals_corrected
