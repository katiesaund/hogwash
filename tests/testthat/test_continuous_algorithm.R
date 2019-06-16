# library(ape)
#
# context("Continuous algorithm") #----------------------------------------------#
#
# # test run_ks_test
# test_that("run_ks_test returns a known test statistic for given data.", {
#   set.seed(1)
#   temp_pheno <- temp_pheno <- matrix(rnorm(20), ncol = 2)
#   num_edges <- nrow(temp_pheno)
#   trans_index <- c(1:8)
#   non_trans_index <- c(9:10)
#   answer <- 0.375
#   names(answer) <- "D"
#   result <- run_ks_test(trans_index, non_trans_index, temp_pheno)
#   expect_equivalent(answer, result$statistic)
#   expect_error(run_ks_test(trans_index, non_trans_index, temp_pheno), NA)
# })
#
# test_that("run_ks_test gives an error when the indices are too few to run ks.test()", {
#   temp_pheno <- temp_pheno <- matrix(1:20, ncol = 2)
#   num_edges <- nrow(temp_pheno)
#   trans_index <- c(1:9)
#   non_trans_index <- ""
#   expect_error(run_ks_test(trans_index, non_trans_index, temp_pheno))
# })
#
# # test get_sig_hits_while_correcting_for_multiple_testing
# test_that("get_sig_hits_while_correcting_for_multiple_testing gives known adjusted p-values for given data", {
#   set.seed(1)
#   fake_p <- rnorm(n = 10, mean = 0.2, sd = 0.2)
#   fake_fdr <- 0.1
#   results <- get_sig_hits_while_correcting_for_multiple_testing(fake_p, fake_fdr)
#   expect_identical(round(results$hit_pvals$fdr_corrected_pvals, 7), c(0.2490308, 0.3862944, 0.1795316, 0.5190562, 0.3862944, 0.1795316, 0.3862944, 0.3862944, 0.3862944, 0.3473058))
#   expect_equal(nrow(results$sig_pvals), 0)
# })
#
# # test calculate_genotype_significance
# test_that("calculate_genotype_significance gives expected results given valid inputs", {
#   num_isolates <- 5
#   num_loci <- 8
#   set.seed(1)
#   temp_tree <- rtree(num_isolates)
#   temp_tree$node.label <- rep(100, Nnode(temp_tree))
#   temp_geno <- matrix(c(0,1), nrow = num_isolates, ncol = num_loci)
#   temp_perm <- 100
#   temp_geno_trans <- temp_conf <- temp_geno_recon <- rep(list(NULL), num_loci)
#   for (i in 1:num_loci) {
#     temp_geno_trans[[i]]$transition <- c(1,0,1,0,1,0,1,0)
#     temp_geno_trans[[i]]$trans_dir <- c(1,0,-1,0,1,0,-1,0)
#     temp_conf[[i]] <- c(1, 1, 1, 0, 0, 1, 1, 1)
#     temp_geno_recon[[i]] <- rep(1, Nedge(temp_tree))
#   }
#   set.seed(1)
#   temp_pheno <- matrix(rnorm(Nedge(temp_tree) * 2), ncol = 2, nrow = Nedge(temp_tree))
#   results <- calculate_genotype_significance(temp_geno, temp_perm, temp_geno_trans, temp_tree, temp_pheno, temp_conf, temp_geno_recon)
#   expect_equal(results$num_genotypes, num_loci)
#   expect_equal(round(results$observed_ks_stat[1], 3), 0.333)
#   expect_equal(round(results$all_edges_median[1], 3), 0.993)
#   expect_equal(round(results$trans_median[1], 2), 1.2)
#   expect_equal(results$observed_pheno_non_trans_delta[[1]], c(0.4890317, 1.3942315, 0.7832583))
#   expect_equal(round(results$observed_pheno_trans_delta[[1]], 3), c(1.202, 2.347, 0.638))
#   expect_equal(round(results$ks_statistics[[1]][1:5], 3), c(0.50, 0.25, 0.75, 0.80, 0.50))
#   expect_equal(round(results$pvals, 3), rep(0.98, num_loci))
# })
#
#
# test_that("calculate_genotype_significance returns a significant p-value when
#           phenotype change is noticeably higher on transition edges than on
#           non-transition tree edges", {
#   num_isolates <- 5
#   num_loci <- 2
#   set.seed(1)
#   temp_tree <- rtree(num_isolates)
#   temp_tree$node.label <- rep(100, Nnode(temp_tree))
#   temp_geno <- matrix(c(0,1), nrow = num_isolates, ncol = num_loci)
#   temp_perm <- 100
#   temp_geno_trans <- temp_conf <- temp_geno_recon <- rep(list(NULL), num_loci)
#   for (i in 1:num_loci) {
#     temp_geno_trans[[i]]$transition <- c(0, 0, 0, 1, 1, 1, 1, 0)
#     temp_geno_trans[[i]]$trans_dir  <- c(1, 1, 1, 1, 1, 1, 1, 1)
#     temp_conf[[i]]                  <- c(1, 1, 1, 1, 1, 1, 1, 1)
#     temp_geno_recon[[i]]            <- c(1, 1, 1, 1, 1, 1, 1, 1)
#   }
#   set.seed(1)
#   temp_pheno <- matrix(c(-1, -1, 10, -.5, 10, 10, 0, 0, 0, -.5, 11, 30, 100, 40, 150, 0.1), ncol = 2, nrow = Nedge(temp_tree))
#   # Transition edges deltas: 10, 90, 30
#   # Non-transition edge deltas: 3, 1, 1
#   results <- calculate_genotype_significance(temp_geno, temp_perm, temp_geno_trans, temp_tree, temp_pheno, temp_conf, temp_geno_recon)
#   alpha <- 0.06
#   expect_true(results$pvals[1] < alpha)
#   expect_equal(results$all_edges_median[[1]], median(abs(temp_pheno[ , 1] - temp_pheno[ , 2])))
#   expect_equal(results$trans_median[[1]], median(abs(temp_pheno[4:7 , 1] - temp_pheno[4:7 , 2])))
#   expect_equal(round(results$observed_ks_stat[1], 2), 1)
#   expect_equal(results$observed_pheno_non_trans_delta[[1]], abs(temp_pheno[c(1:3,8), 1] - temp_pheno[c(1:3,8), 2]))
#   expect_equal(round(results$observed_pheno_trans_delta[[1]], 3), abs(temp_pheno[4:7, 1] - temp_pheno[4:7, 2]))
#   expect_equal(round(results$ks_statistics[[1]][1:5], 3), c(0.40, 0.50, 1.00, 0.25, 0.25))
# })
#
# test_that("calculate_genotype_significance returns a non-significant p-value when
#           phenotype change is identical on transition and non-transition edges", {
#   num_isolates <- 4
#   num_loci <- 2
#   set.seed(1)
#   temp_tree <- rtree(num_isolates)
#   temp_tree$node.label <- rep(100, Nnode(temp_tree))
#   temp_geno <- matrix(c(0,1), nrow = num_isolates, ncol = num_loci)
#   temp_perm <- 100
#   temp_geno_trans <- temp_conf <- temp_geno_recon <- rep(list(NULL), num_loci)
#   for (i in 1:num_loci) {
#     temp_geno_trans[[i]]$transition <- c(0, 0, 0, 1, 1, 1)
#     temp_geno_trans[[i]]$trans_dir  <- c(1, 1, 1, 1, 1, 1)
#     temp_conf[[i]]                  <- c(1, 1, 1, 1, 1, 1)
#     temp_geno_recon[[i]]            <- c(1, 1, 1, 1, 1, 1)
#   }
#   set.seed(1)
#   temp_pheno <- matrix(c(-15, -15, 10, 0, 10, 10, -10, 0, 7, 5, 25, 7), ncol = 2, nrow = Nedge(temp_tree))
#   # Transition edges deltas: 5 15 3
#   # Non-transition edge deltas: 5 15 3
#   results <- calculate_genotype_significance(temp_geno, temp_perm, temp_geno_trans, temp_tree, temp_pheno, temp_conf, temp_geno_recon)
#   alpha <- 0.01
#   expect_true(results$pvals[1] > alpha)
#   expect_equal(results$all_edges_median[[1]], median(abs(temp_pheno[ , 1] - temp_pheno[ , 2])))
#   expect_equal(results$trans_median[[1]], median(abs(temp_pheno[4:6 , 1] - temp_pheno[4:6 , 2])))
#   expect_equal(round(results$observed_ks_stat[1], 2), 0)
#   expect_equal(results$observed_pheno_non_trans_delta[[1]], abs(temp_pheno[1:3 , 1] - temp_pheno[1:3 , 2]))
#   expect_equal(round(results$observed_pheno_trans_delta[[1]], 3), abs(temp_pheno[4:6 , 1] - temp_pheno[4:6 , 2]))
# })
#
# test_that("calculate_genotype_significance returns an error when the only
#           confident edges are transition edges. (All non-transition edges are
#           low confidence)", {
#   num_isolates <- 4
#   num_loci <- 2
#   set.seed(1)
#   temp_tree <- rtree(num_isolates)
#   temp_tree$node.label <- rep(100, Nnode(temp_tree))
#   temp_geno <- matrix(c(0,1), nrow = num_isolates, ncol = num_loci)
#   temp_perm <- 100
#   temp_geno_trans <- temp_conf <- temp_geno_recon <- rep(list(NULL), num_loci)
#   for (i in 1:num_loci) {
#     temp_geno_trans[[i]]$transition <- c(0, 0, 0, 1, 1, 1)
#     temp_geno_trans[[i]]$trans_dir  <- c(1, 1, 1, 1, 1, 1)
#     temp_conf[[i]]                  <- c(0, 0, 0, 1, 1, 1)
#     temp_geno_recon[[i]]            <- c(1, 1, 1, 1, 1, 1)
#   }
#   set.seed(1)
#   temp_pheno <- matrix(c(-15, -15, 10, 0, 10, 10, -10, 0, 7, 5, 25, 7), ncol = 2, nrow = Nedge(temp_tree))
#   # Transition edges deltas: 5 15 3
#   # Non-transition edge deltas: 5 15 3
#   expect_error(calculate_genotype_significance(temp_geno, temp_perm, temp_geno_trans, temp_tree, temp_pheno, temp_conf, temp_geno_recon))
# })
#
#
#
# # test get_hi_conf_tran_indices
# test_that("get_hi_conf_tran_indices returns only high confidence transition edges given this test data", {
#   num_isolates <- 5
#   num_loci <- 8
#   set.seed(1)
#   temp_tree <- rtree(num_isolates)
#   temp_tree$node.label <- rep(100, Nnode(temp_tree))
#   temp_geno_trans <- temp_conf <- rep(list(NULL), num_loci)
#   for (i in 1:num_loci) {
#     temp_geno_trans[[i]]$transition <- c(1,0,1,0,1,0,1,0)
#     temp_geno_trans[[i]]$trans_dir <- c(1,0,-1,0,1,0,-1,0)
#     temp_conf[[i]] <- c(1, 1, 1, 0, 0, 1, 1, 1)
#   }
#   index <- 1
#   indices <- get_hi_conf_tran_indices(temp_geno_trans, temp_conf, index, temp_tree)
#   expect_equal(indices$trans_index, c(1, 3, 7))
#   expect_equal(indices$non_trans_index, c(2, 6, 8))
# })
#
#
# # test continuous_permutation
# test_that("continuous_permutation is gives consistent results with this test set", {
#   num_isolates <- 40
#   num_loci <- 80
#   set.seed(1)
#   temp_tree <- rtree(num_isolates)
#   temp_tree$node.label <- rep(100, Nnode(temp_tree))
#   temp_conf <- rep(list(NULL), num_loci)
#   for (i in 1:num_loci) {
#     temp_conf[[i]] <- rep(c(1, 0), Nedge(temp_tree)/2)
#   }
#   num <- 1
#   index <- NULL
#   index$trans_index <- c(1:40)
#   perm <- 100000
#   results <- continuous_permutation(index, temp_tree, temp_conf, perm, num)
#   odd_numbers <- seq(from = 1, to = 77, by = 2)
#   expect_equal(results$hi_conf_edges, odd_numbers)
#   expect_equal(nrow(results$permuted_trans_index_mat), perm)
#   expect_equal(ncol(results$permuted_trans_index_mat), length(index$trans_index))
#   expect_equal(results$permuted_trans_index_mat[1, 1:10], c(18, 8, 19, 5, 23, 25, 20, 21, 15, 38))
# })
#
# # test calculate_empirical_pheno_delta
# test_that("calculate_empirical_pheno_delta returns a list of valid ks-test results of length perm", {
#   perm <- 1000
#   num_isolates <- 40
#   nedge <- 78
#   set.seed(1)
#   temp_permuted_mat <- matrix(round(rnorm(perm * num_isolates, mean = 30, sd = 5), 0), nrow = perm, ncol = num_isolates)
#   set.seed(1)
#   temp_pheno <- matrix(rnorm(2 * nedge), nrow = nedge, ncol = 2)
#   temp_conf <- seq(from = 1, to = max(temp_permuted_mat), by = 2)
#   results <- calculate_empirical_pheno_delta(perm, temp_permuted_mat, temp_conf, temp_pheno)
#   length(results)
#   expect_equal(length(results), perm)
#   expect_true(max(results) < 1)
#   expect_true(min(results) > 0)
# })
#