library(microbeGWAS)
context("Continuous algorithm") -----------------------------------------------#

# test run_ks_test
test_that("run_ks_test returns a known test statistic for given data.", {
  set.seed(1)
  temp_pheno <- temp_pheno <- matrix(rnorm(20), ncol = 2)
  num_edges <- nrow(temp_pheno)
  trans_index <- c(1:8)
  non_trans_index <- c(9:10)
  answer <- 0.375
  names(answer) <- "D"
  result <- run_ks_test(trans_index, non_trans_index, temp_pheno)
  expect_equivalent(answer, result$statistic)
  expect_error(run_ks_test(trans_index, non_trans_index, temp_pheno), NA)
})

test_that("run_ks_test gives an error when the indices are too few to run ks.test()", {
  temp_pheno <- temp_pheno <- matrix(1:20, ncol = 2)
  num_edges <- nrow(temp_pheno)
  trans_index <- c(1:9)
  non_trans_index <- ""
  expect_error(run_ks_test(trans_index, non_trans_index, temp_pheno))
})

# test get_sig_hits_while_correcting_for_multiple_testing
test_that("get_sig_hits_while_correcting_for_multiple_testing gives known adjusted p-values for given data", {
  set.seed(1)
  fake_p <- rnorm(n = 10, mean = 0.2, sd = 0.2)
  fake_fdr <- 0.1
  results <- get_sig_hits_while_correcting_for_multiple_testing(fake_p, fake_fdr)
  expect_identical(round(results$hit_pvals$fdr_corrected_pvals, 7), c(0.2490308, 0.3862944, 0.1795316, 0.5190562, 0.3862944, 0.1795316, 0.3862944, 0.3862944, 0.3862944, 0.3473058))
  expect_equal(nrow(results$sig_pvals), 0)
})

test_that("calculate_genotype_significance does X given Y", {
  num_isolates <- 5
  num_loci <- 8
  set.seed(1)
  temp_tree <- rtree(num_isolates)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_geno <- matrix(c(0,1), nrow = num_isolates, ncol = num_loci)
  temp_perm <- 100
  temp_geno_trans <- temp_conf <- temp_geno_recon <- rep(list(NULL), num_loci)
  for (i in 1:num_loci){
    temp_geno_trans[[i]]$transition <- c(1,0,1,0,1,0,1,0)
    temp_geno_trans[[i]]$trans_dir <- c(1,0,-1,0,1,0,-1,0)
    temp_conf[[i]] <- c(1, 1, 1, 0, 0, 1, 1, 1)
    temp_geno_recon[[i]] <- rep(1, Nedge(temp_tree))
  }
  set.seed(1)
  temp_pheno <- matrix(rnorm(Nedge(temp_tree) * 2), ncol = 2, nrow = Nedge(temp_tree))
  results <- calculate_genotype_significance(temp_geno, temp_perm, temp_geno_trans, temp_tree, temp_pheno, temp_conf, temp_geno_recon)
  expect_equal(results$num_genotypes, num_loci)
  expect_equal(round(results$observed_ks_stat[1], 3), 0.333)
  expect_equal(round(results$all_edges_median[1], 3), 0.993)
  expect_equal(round(results$trans_median[1], 2), 1.2)
  expect_equal(results$observed_pheno_non_trans_delta[[1]], c(0.4890317, 1.3942315, 0.7832583))
  expect_equal(round(results$observed_pheno_trans_delta[[1]], 3), c(1.202, 2.347, 0.638))
  expect_equal(round(results$ks_statistics[[1]][1:5], 3), c(0.50, 0.25, 0.75, 0.80, 0.50))
  expect_equal(round(results$pvals, 3), rep(0.812, num_loci))
})

test_that("get_hi_conf_tran_indices returns only high confidence transition edges given this test data", {
  num_isolates <- 5
  num_loci <- 8
  set.seed(1)
  temp_tree <- rtree(num_isolates)
  temp_tree$node.label <- rep(100, Nnode(temp_tree))
  temp_geno_trans <- temp_conf <- rep(list(NULL), num_loci)
  for (i in 1:num_loci){
    temp_geno_trans[[i]]$transition <- c(1,0,1,0,1,0,1,0)
    temp_geno_trans[[i]]$trans_dir <- c(1,0,-1,0,1,0,-1,0)
    temp_conf[[i]] <- c(1, 1, 1, 0, 0, 1, 1, 1)
  }
  index <- 1
  indices <- get_hi_conf_tran_indices(temp_geno_trans, temp_conf, index, temp_tree)
  expect_equal(indices$trans_index, c(1, 3, 7))
  expect_equal(indices$non_trans_index, c(2, 6, 8))
})

# calculate_genotype_significance <- function(mat, permutations, genotype_transition_list, tr, pheno_recon_ordered_by_edges, genotype_confidence, genotype_reconstruction){
#   # Function description -------------------------------------------------------
#   #
#   # Input:
#   # mat is a genotype_matrix. Each row is a genotype. Each column is a 0/1 value on an edge.
#   # permutations is a numeric integer.
#   # genotype_transition_list is a list of length = ncol(mat) and each list has length = Nedge(tr). All entries are 0 or 1.
#   # tr object of phylo.
#   # phenotype_recon_ordered_by_edges is a matrix. Nrow = Nedge(tr). Ncol =2. where the phenotype reconstruction is ordered by edges. IN THE SAME FORMAT AS tr$EDGE. SO NODE VALUES WILL APPEAR MULTIPLE TIMES IN THE tr.
#   #            THIS FORMAT WILL MAKE IT MUCH EASIER TO CALCULATE PHENOTYPE CHANGE ON EDGES.
#   # genotype_confidence is a list of length = ncol(mat) and each list has length = Nedge(tr). All entries are 0 (low confidence) or 1 (high confidence).
#   #            NOTE: genotype_confidence lists the confidence in each edge. High confidence means the edge is high confidence by genotype reconstruction, phenotype reconstruction, bootstrap value, and edge length.
#   #
#   # Outputs:
#   # results. List of 8.
#   #          $pvals. Named numeric vector. Length == number of genotypes. Values between 1 and 0. Names are genotype names.
#   #          $ks_statistics. List of numeric vectors. Length of list == number of genotypes. Each vector has length == number of permutations. Values between 1 and 0.
#   #          $observed_pheno_trans_delta. List of numeric vectors. Length of list == number of genotypes. Vectors are of variable length because length is the number of transition edges for that particular genotype. Vectors are numeric.
#   #          $observed_pheno_non_trans_delta. List of numeric vectors. Length of list == number of genotypes. Vectors are of variable length because length is the number of non-transition edges for that particular genotype. Vectors are numeric.
#   #          $trans_median. Numberic. Vector. Length = number of genotypes. Describes median delta phenotype on all transition edges.
#   #          $all_edges_median. Numeric vector. Length = number of genotypes. Describes median delta phenotype on all edges.
#   #          $num_genotypes. Integer. The number of genotypes.
#   #          $observed_ks_stat. Numeric Vector. Length = number of genotypes. Values between 1 and 0.
#   #
#
#   # Check input ----------------------------------------------------------------
#   check_for_root_and_bootstrap(tr)
#   check_if_permutation_num_valid(permutations)
#
#   # Function -------------------------------------------------------------------
#   num_genotypes <- ncol(mat)
#   pvals <- observed_ks_pval <- trans_median <- all_edges_median <- observed_ks_stat <- rep(NA, num_genotypes) # 2018-11-12 added observed_ks_stat
#   names(observed_ks_pval) <- names(pvals) <- colnames(mat)
#   empirical_ks_pval_list <- empirical_ks_stat_list <- observed_pheno_trans_delta <- observed_pheno_non_trans_delta <- rep(list(0), num_genotypes) # 2018-11-12 added empirical_ks_stat_list
#
#   for (i in 1:num_genotypes){
#     # GRAB THE IDS OF THE TRANSITION EDGES:
#     trans_index     <- c(1:Nedge(tr))[as.logical(genotype_transition_list[[i]]$transition)]
#     non_trans_index <- c(1:Nedge(tr))[!genotype_transition_list[[i]]$transition]
#     # [1] EX:  8 12 13 16 19 26 27 31 37 44 52 56 64 67 68 76 77 80 89 92 97 98
#     # THESE EDGES ARE DEFINED BY THE NODES IN THE CORRESPONDING ROWS OF tr$EDGE
#
#     # SUBSET TRANSITION / NON-TRANSITION EDGES TO ONLY HIGH CONFIDENCE ONES
#     hi_conf_trans_index     <- trans_index[    as.logical(genotype_confidence[[i]][trans_index    ])]
#     hi_conf_non_trans_index <- non_trans_index[as.logical(genotype_confidence[[i]][non_trans_index])]
#
#     # Run KS test to find out if the phenotype change on transition edges is significantly different from phenotype change on non-transition edges
#     observed_results <- run_ks_test(hi_conf_trans_index, hi_conf_non_trans_index, pheno_recon_ordered_by_edges)
#     observed_ks_pval[i] <- observed_results$pval
#     observed_ks_stat[i] <- observed_results$statistic
#
#     # save these for reporting / plots
#     observed_pheno_trans_delta[[i]]     <- observed_results$pheno_trans_delta
#     observed_pheno_non_trans_delta[[i]] <- observed_results$pheno_non_trans_delta
#     trans_median[i]     <- median(observed_results$pheno_trans_delta)
#     all_edges_median[i] <- median(c(observed_results$pheno_trans_delta, observed_results$pheno_non_trans_delta))
#     #
#
#     # do the permutation part
#     num_sample           <- length(hi_conf_trans_index)
#     all_edges            <- c(1:Nedge(tr))
#     which_branches       <- all_edges[as.logical(genotype_confidence[[i]])]
#     all_sampled_branches <- matrix(   nrow = permutations, ncol = num_sample)
#     redistributed_hits   <- matrix(0, nrow = permutations, ncol = length(which_branches))
#
#     # ii. Create a matrix where each row is a random sampling of the high
#     #     confidence edges of the tr where probability of choosing edges is
#     #     proportional to length of edge. Number of edges selected for the
#     #     permuted data set is the number of times the empirical genotype
#     #     appears.
#
#     set.seed(1) # for reproducability of the sample() function
#     for (j in 1:permutations){  # create a random sample of the tr
#       curr_num_branch <- length(which_branches)
#       all_sampled_branches[j, ] <- sample(1:curr_num_branch,
#                                           size = num_sample,
#                                           replace = TRUE,
#                                           prob = tr$edge.length[which_branches]/sum(tr$edge.length[which_branches]))
#     } # end for (j)
#
#     # all_sampled_branches is my new, permuted "hi_conf_trans_index" where each row is a new list of transition genotype branches
#     # BUT CAVEAT: these are just fake/null transitions and some of them are probably actually touching! If I wanted to be
#     # Super legit I would recreate as many hits, calculate new transitions, and then use those in my permutation test, somehow
#     # controlling for variable numbers of transitions. But not doing that for now.
#
#     # calculate permuted pheno deltas
#     empirical_ks_pval <- empirical_ks_stat <- rep(NA, permutations)
#     for (k in 1:nrow(all_sampled_branches)){
#       permuted_trans_index     <- unique(all_sampled_branches[k, ])
#       permuted_non_trans_index <- c(1:length(which_branches))[!(c(1:length(which_branches)) %in% unique(all_sampled_branches[k, ]))]
#       empirical_results        <- run_ks_test(permuted_trans_index, permuted_non_trans_index, pheno_recon_ordered_by_edges)
#       empirical_ks_pval[k]     <- empirical_results$pval
#       empirical_ks_stat[k]     <- empirical_results$statistic
#     } # end for (k)
#
#     # the observed ks.test statistic: observed_ks_stat[i] (fixed this 2018-11-12; before I had it wrong with observed_ks_pval)
#     # empirical p value caluclation here: (1 + more extreme observations) / (1 + permutations)
#     empirical_ks_pval_list[[i]] <- empirical_ks_pval
#     empirical_ks_stat_list[[i]] <- empirical_ks_stat
#     pvals[i] <- (sum(1 + sum(empirical_ks_stat > observed_ks_stat[i]))/(permutations + 1))
#   } # end for (i)
#
#   # Return output --------------------------------------------------------------
#   results <- list("pvals" = pvals, "ks_statistics" = empirical_ks_stat_list,
#                   "observed_pheno_trans_delta" = observed_pheno_trans_delta,
#                   "observed_pheno_non_trans_delta" = observed_pheno_non_trans_delta,
#                   "trans_median" = trans_median, "all_edges_median" = all_edges_median,
#                   "num_genotypes" = ncol(mat), "observed_ks_stat" = observed_ks_stat) # 2018-11-28
#   return(results)
# } # end calculate_genotype_significance()

