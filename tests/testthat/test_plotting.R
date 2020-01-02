context("Plots") #----------------------------------------------#

# test plot_continuous_phenotype()
test_that("plot_continuous_phenotype works for valid inputs", {
  ntip <- 10
  temp_tree <- ape::rtree(ntip)
  temp_vec <- rnorm(ntip, mean = 0, sd = 10)
  names(temp_vec) <- temp_tree$tip.label
  temp_anc_rec <-
    ape::ace(x = temp_vec, phy = temp_tree, type = "continuous", method = "ML")
  temp_anc_rec_at_node <- temp_anc_rec$ace
  expect_error(plot_continuous_phenotype(temp_tree,
                                         temp_vec,
                                         temp_anc_rec_at_node),
               NA)
})

test_that("plot_continuous_phenotype works gives error for invalid inputs
  (no names on phenotype vector)", {
  ntip <- 10
  temp_tree <- ape::rtree(ntip)
  temp_vec <- rnorm(ntip, mean = 0, sd = 10)
  temp_anc_rec <-
    ape::ace(x = temp_vec, phy = temp_tree, type = "continuous", method = "ML")
  temp_anc_rec_at_node <- temp_anc_rec$ace
  expect_error(plot_continuous_phenotype(temp_tree,
                                         temp_vec,
                                         temp_anc_rec_at_node))
})

test_that("plot_continuous_phenotype gives error for invalid inputs", {
  ntip <- 10
  temp_tree <- ape::rtree(ntip)
  temp_vec <- rnorm(ntip, mean = 0, sd = 10)
  names(temp_vec) <- temp_tree$tip.label
  temp_anc_rec <-
    ape::ace(x = temp_vec, phy = temp_tree, type = "continuous", method = "ML")
  temp_anc_rec_at_node <- temp_anc_rec$ace[1:3]
  expect_error(plot_continuous_phenotype(temp_tree,
                                         temp_vec,
                                         temp_anc_rec_at_node))
})

test_that("plot_continuous_phenotype gives error for invalid inputs", {
  ntip <- 10
  temp_tree <- ape::rtree(ntip)
  temp_vec <- rnorm(ntip, mean = 0, sd = 10)
  names(temp_vec) <- temp_tree$tip.label
  temp_anc_rec <-
    ape::ace(x = temp_vec, phy = temp_tree, type = "continuous", method = "ML")
  temp_anc_rec_at_node <- temp_anc_rec$ace
  temp_vec <- temp_vec[1:5]
  expect_error(plot_continuous_phenotype(temp_tree,
                                         temp_vec,
                                         temp_anc_rec_at_node))
})

# test hist_raw_conf_delta_pheno
test_that("hist_raw_conf_delta_pheno works for valid inputs", {
  temp_trans <- temp_conf <- NULL
  temp_trans[[1]]$transition <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_trans[[1]]$trans_dir <- c(1, 1, 1, 1, 1, 1, 1, 1)
  temp_conf[[1]] <- c(0, 0, 1, 1, 1, 1, 1, 1)
  temp_pheno <- matrix(c(c(1, 2, 3, 4, 1, 2, 3, 4),
                         c(2, 3, 4, 5, 0, 1, 2, 3)),
                       ncol = 2,
                       nrow = 8)
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 1
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"

  expect_error(hist_raw_hi_conf_delta_pheno(temp_trans,
                                            temp_conf,
                                            temp_pheno,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color),
               NA)
})

test_that("hist_raw_conf_delta_pheno works gives error for invalid inputs", {
  temp_trans <- temp_conf <- NULL
  temp_trans[[1]]$transition <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_trans[[1]]$trans_dir <- c(1, 1, 1, 1, 1, 1, 1, 1)
  temp_conf[[1]] <- c(0, 0, 1, 1, 1, 1, 1, 1)
  temp_pheno <- matrix(c(c(1, 2, 3, 4, 1, 2, 3, 4),
                         c(2, 3, 4, 5, 0, 1, 2, 3)),
                       ncol = 2,
                       nrow = 8)
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 1
  temp_non_trans_color <- "black"
  temp_trans_color <- "black"

  expect_error(hist_raw_hi_conf_delta_pheno(temp_trans,
                                            temp_conf,
                                            temp_pheno,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

test_that("hist_raw_conf_delta_pheno works gives error for invalid inputs", {
  temp_trans <- temp_conf <- NULL
  temp_trans[[1]]$transition <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_trans[[1]]$trans_dir <- c(1, 1, 1, 1, 1, 1, 1, 1)
  temp_conf[[1]] <- c(0, 0, 1, 1, 1, 1, 1, 1)
  temp_pheno <- matrix(c(c(1, 2, 3, 4, 1, 2, 3, 4),
                         c(2, 3, 4, 5, 0, 1, 2, 3)),
                       ncol = 2,
                       nrow = 8)
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 10
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"

  expect_error(hist_raw_hi_conf_delta_pheno(temp_trans,
                                            temp_conf,
                                            temp_pheno,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

test_that("hist_raw_conf_delta_pheno works gives error for invalid inputs", {
  temp_trans <- temp_conf <- NULL
  temp_trans[[1]]$transition <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_trans[[1]]$trans_dir <- c(1, 1, 1, 1, 1, 1, 1, 1)
  temp_conf[[1]] <- c(0, 0, 1, 1, 1, 1, 1, 1)
  temp_pheno <- matrix(c(c(1, 2, 3, 4, 1, 2, 3, 4),
                         c(2, 3, 4, 5, 0, 1, 2, 3)),
                       ncol = 2,
                       nrow = 8)
  temp_tree <- "foobar"
  temp_index <- 1
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"

  expect_error(hist_raw_hi_conf_delta_pheno(temp_trans,
                                            temp_conf,
                                            temp_pheno,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

test_that("hist_raw_conf_delta_pheno works gives error for invalid inputs", {
  temp_trans <- temp_conf <- NULL
  temp_trans[[1]]$transition <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_trans[[1]]$trans_dir <- c(1, 1, 1, 1, 1, 1, 1, 1)
  temp_conf[[1]] <- c(0, 0, 1, 1)
  temp_pheno <- matrix(c(c(1, 2, 3, 4, 1, 2, 3, 4),
                         c(2, 3, 4, 5, 0, 1, 2, 3)),
                       ncol = 2,
                       nrow = 8)
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 1
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"

  expect_error(hist_raw_hi_conf_delta_pheno(temp_trans,
                                            temp_conf,
                                            temp_pheno,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

test_that("hist_raw_conf_delta_pheno works gives error for invalid inputs", {
  temp_trans <- temp_conf <- NULL
  temp_trans[[1]]$transition <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_trans[[1]]$trans_dir <- c(1, 1, 1, 1, 1, 1, 1, 1)
  temp_conf[[1]] <- c(0, 0, 1, 1, 1, 1, 1, 1)
  temp_pheno <- matrix(c(c(1, 2, 3, 4, 1, 2),
                         c(2, 3, 4, 5, 0, 1)),
                       ncol = 2,
                       nrow = 6)
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 1
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"

  expect_error(hist_raw_hi_conf_delta_pheno(temp_trans,
                                            temp_conf,
                                            temp_pheno,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

test_that("hist_raw_conf_delta_pheno works gives error for invalid inputs", {
  temp_trans <- temp_conf <- NULL
  temp_trans[[1]]$transition <- c(1, 1)
  temp_trans[[1]]$trans_dir <- c(1, 1, 1, 1, 1, 1, 1, 1)
  temp_conf[[1]] <- c(0, 0, 1, 1, 1, 1, 1, 1)
  temp_pheno <- matrix(c(c(1, 2, 3, 4, 1, 2, 3, 4),
                         c(2, 3, 4, 5, 0, 1, 2, 3)),
                       ncol = 2,
                       nrow = 8)
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 1
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"

  expect_error(hist_raw_hi_conf_delta_pheno(temp_trans,
                                            temp_conf,
                                            temp_pheno,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

# test hist_abs_hi_conf_delta_pheno
test_that("hist_abs_hi_conf_delta_pheno works for valid inputs", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 1
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"
  temp_trans <- NULL
  temp_trans$observed_pheno_non_trans_delta[[1]] <- c(0, 1, 1, 2)
  temp_trans$observed_pheno_trans_delta[[1]] <- c(2, 3, 4, 5)
  expect_error(hist_abs_hi_conf_delta_pheno(temp_trans,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color),
               NA)
})

test_that("hist_abs_hi_conf_delta_pheno gives errors for invalid inputs", {
  set.seed(1)
  temp_tree <- "tree"
  temp_index <- 1
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"
  temp_trans <- NULL
  temp_trans$observed_pheno_non_trans_delta[[1]] <- c(0, 1, 1, 2)
  temp_trans$observed_pheno_trans_delta[[1]] <- c(2, 3, 4, 5)
  expect_error(hist_abs_hi_conf_delta_pheno(temp_trans,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

test_that("hist_abs_hi_conf_delta_pheno gives errors for invalid inputs", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 10
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"
  temp_trans <- NULL
  temp_trans$observed_pheno_non_trans_delta[[1]] <- c(0, 1, 1, 2)
  temp_trans$observed_pheno_trans_delta[[1]] <- c(2, 3, 4, 5)
  expect_error(hist_abs_hi_conf_delta_pheno(temp_trans,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

test_that("hist_abs_hi_conf_delta_pheno gives errors for invalid inputs", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 1
  temp_non_trans_color <- "orange"
  temp_trans_color <- "orange"
  temp_trans <- NULL
  temp_trans$observed_pheno_non_trans_delta[[1]] <- c(0, 1, 1, 2)
  temp_trans$observed_pheno_trans_delta[[1]] <- c(2, 3, 4, 5)
  expect_error(hist_abs_hi_conf_delta_pheno(temp_trans,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

test_that("hist_abs_hi_conf_delta_pheno gives errors for invalid inputs", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_index <- 1
  temp_non_trans_color <- "black"
  temp_trans_color <- "orange"
  temp_trans <- NULL
  temp_trans$observed_pheno_non_trans_delta[[2]] <- c(0, 1, 1, 2)
  temp_trans$observed_pheno_trans_delta[[2]] <- c(2, 3, 4, 5)
  expect_error(hist_abs_hi_conf_delta_pheno(temp_trans,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
})

# test hist_abs_delta_pheno_all_edges
test_that("hist_abs_delta_pheno_all_edges works for valid inputs", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_pheno_mat <- matrix(abs(rnorm(ape::Nedge(temp_tree), 0, 10)), ncol = 1)
  temp_index <- 1
  temp_conf <- NULL
  temp_conf[[1]] <- c(0, 0, 0, 1, 1, 1, 1, 1)
  expect_error(hist_abs_delta_pheno_all_edges(temp_pheno_mat,
                                              temp_conf,
                                              temp_tree,
                                              temp_index),
               NA)
})

test_that("hist_abs_delta_pheno_all_edges gives error for invalid inputs", {
  temp_tree <- "tree"
  temp_pheno_mat <- matrix(abs(rnorm(8, 0, 10)), ncol = 1)
  temp_index <- 1
  temp_conf <- NULL
  temp_conf[[1]] <- c(0, 0, 0, 1, 1, 1, 1, 1)
  expect_error(hist_abs_delta_pheno_all_edges(temp_pheno_mat,
                                              temp_conf,
                                              temp_tree,
                                              temp_index))
})


test_that("hist_abs_delta_pheno_all_edges gives error for invalid inputs", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_pheno_mat <- c(0, 0, 1)
  temp_index <- 1
  temp_conf <- NULL
  temp_conf[[1]] <- c(0, 0, 0, 1, 1, 1, 1, 1)
  expect_error(hist_abs_delta_pheno_all_edges(temp_pheno_mat,
                                              temp_conf,
                                              temp_tree,
                                              temp_index))
})

test_that("hist_abs_delta_pheno_all_edges gives error for invalid inputs", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_pheno_mat <- matrix(abs(rnorm(ape::Nedge(temp_tree), 0, 10)), ncol = 1)
  temp_index <- 10
  temp_conf <- NULL
  temp_conf[[1]] <- c(0, 0, 0, 1, 1, 1, 1, 1)
  expect_error(hist_abs_delta_pheno_all_edges(temp_pheno_mat,
                                              temp_conf,
                                              temp_tree,
                                              temp_index))
})

test_that("hist_abs_delta_pheno_all_edges gives error for invalid inputs", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_pheno_mat <- matrix(abs(rnorm(ape::Nedge(temp_tree), 0, 10)), ncol = 1)
  temp_index <- 1
  temp_conf <- NULL
  temp_conf[[1]] <- 0
  expect_error(hist_abs_delta_pheno_all_edges(temp_pheno_mat,
                                              temp_conf,
                                              temp_tree,
                                              temp_index))
})

# test make_manhattan_plot
test_that("make_manhattan_lot works for valid input", {
  # Set up
  temp_name <- "test"
  temp_pval <- seq(from = 0.05, to = 0.95, by = 0.05)
  temp_fdr <- 0.15
  temp_type <- "synchronous"

  # Test
  expect_error(make_manhattan_plot(temp_name, temp_pval, temp_fdr, temp_type),
               NA)
})

test_that("make_manhattan_lot gives error for invalid input", {
  # Set up
  temp_name <- 5
  temp_pval <- seq(from = 0.05, to = 0.95, by = 0.05)
  temp_fdr <- 0.15
  temp_type <- "synchronous"

  # Test
  expect_error(make_manhattan_plot(temp_name, temp_pval, temp_fdr, temp_type))
})

test_that("make_manhattan_lot gives error for invalid input", {
  # Set up
  temp_name <- "test"
  temp_pval <- letters[1:5]
  temp_fdr <- 0.15
  temp_type <- "synchronous"

  # Test
  expect_error(make_manhattan_plot(temp_name, temp_pval, temp_fdr, temp_type))
})

test_that("make_manhattan_lot gives error for invalid input", {
  # Set up
  temp_name <- "test"
  temp_pval <- seq(from = 0.05, to = 0.95, by = 0.05)
  temp_fdr <- "foobar"
  temp_type <- "synchronous"

  # Test
  expect_error(make_manhattan_plot(temp_name, temp_pval, temp_fdr, temp_type))
})


test_that("make_manhattan_lot gives error for invalid input", {
  # Set up
  temp_name <- "test"
  temp_pval <- seq(from = 0.05, to = 0.95, by = 0.05)
  temp_fdr <- 0.15
  temp_type <- 5

  # Test
  expect_error(make_manhattan_plot(temp_name, temp_pval, temp_fdr, temp_type))
})

# test plot_tr_w_color_edges
test_that("plot_tr_w_color_edges works for valid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "recon"
  temp_index <- 1
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other),
               NA)
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- "foobar"
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "recon"
  temp_index <- 1
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "recon"
  temp_index <- 1
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "recon"
  temp_index <- 1
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "purple"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "recon"
  temp_index <- 1
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- 10
  temp_title <- "test_title"
  temp_type <- "recon"
  temp_index <- 1
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- matrix(0, 1, 1)
  temp_type <- "recon"
  temp_index <- 1
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "foobar"
  temp_index <- 1
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "recon"
  temp_index <- 10
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "recon"
  temp_index <- 1
  temp_legend_base <- c(1, 0, 1)
  temp_legend_other <- "hi conf presence"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges gives error for invalid inputs on recon", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]] <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "recon"
  temp_index <- 1
  temp_legend_base <- "hi conf absence"
  temp_legend_other <- matrix(0, 1, 1)

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other))
})

test_that("plot_tr_w_color_edges works for valid inputs on trans", {
  # Set up
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_edges_to_highlight <- NULL
  temp_edges_to_highlight[[1]]$transition <- c(1, 1, 1, 1, 0, 0, 0, 0)
  temp_conf <- NULL
  temp_conf[[1]] <- c(1, 1, 0, 0, 1, 1, 0, 0)
  temp_low_conf_color <- "grey"
  temp_bright_color <- "purple"
  temp_title <- "test_title"
  temp_type <- "trans"
  temp_index <- 1
  temp_legend_base <- "hi conf non-transition"
  temp_legend_other <- "hi conf transition"

  # Test
  expect_error(plot_tr_w_color_edges(temp_tree,
                                     temp_edges_to_highlight,
                                     temp_conf,
                                     temp_low_conf_color,
                                     temp_bright_color,
                                     temp_title,
                                     temp_type,
                                     temp_index,
                                     temp_legend_base,
                                     temp_legend_other),
               NA)
})

# blank_plot
test_that("blank_plot gives no errors", {
  expect_error(blank_plot(), NA)
})
