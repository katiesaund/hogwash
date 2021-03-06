context("Plots") #----------------------------------------------#

# test plot_continuous_phenotype()
test_that("plot_continuous_phenotype works for valid inputs", {
  ntip <- 10
  temp_tree <- ape::rtree(ntip)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_vec <- rnorm(ntip, mean = 0, sd = 10)
  names(temp_vec) <- temp_tree$tip.label
  temp_tr_type <- "phylogram"
  temp_anc_rec <-
    suppressWarnings(ape::ace(x = temp_vec,
                              phy = temp_tree,
                              type = "continuous",
                              method = "REML"))
  temp_anc_rec_at_node <- temp_anc_rec$ace
  expect_error(plot_continuous_phenotype(temp_tree,
                                         temp_vec,
                                         temp_anc_rec_at_node,
                                         temp_tr_type,
                                         strain_key = NULL),
               NA)


  temp_strain_key <- matrix(c("one", "one", "two", "two", "three",
                              "four", "four", "four", "four", "four"),
                            nrow = ntip,
                            ncol = 1)
  row.names(temp_strain_key) <- temp_tree$tip.label
  expect_error(plot_continuous_phenotype(temp_tree,
                                         temp_vec,
                                         temp_anc_rec_at_node,
                                         temp_tr_type,
                                         strain_key = temp_strain_key),
               NA)
})

test_that("plot_continuous_phenotype works gives error for invalid inputs
  (no names on phenotype vector)", {
  ntip <- 10
  temp_tree <- ape::rtree(ntip)
  temp_vec <- rnorm(ntip, mean = 0, sd = 10)
  temp_tr_type <- "phylogram"
  temp_anc_rec <-
    suppressWarnings(ape::ace(x = temp_vec,
                              phy = temp_tree,
                              type = "continuous",
                              method = "REML"))
  temp_anc_rec_at_node <- temp_anc_rec$ace
  expect_error(plot_continuous_phenotype(temp_tree,
                                         temp_vec,
                                         temp_anc_rec_at_node,
                                         temp_tr_type,
                                         strain_key = NULL))
})

test_that("plot_continuous_phenotype gives error for invalid inputs", {
  ntip <- 10
  temp_tree <- ape::rtree(ntip)
  temp_vec <- rnorm(ntip, mean = 0, sd = 10)
  names(temp_vec) <- temp_tree$tip.label
  temp_anc_rec <-
    suppressWarnings(ape::ace(x = temp_vec,
                              phy = temp_tree,
                              type = "continuous",
                              method = "REML"))
  temp_anc_rec_at_node <- temp_anc_rec$ace[1:3]
  temp_tr_type <- "phylogram"
  expect_error(plot_continuous_phenotype(temp_tree,
                                         temp_vec,
                                         temp_anc_rec_at_node,
                                         temp_tr_type,
                                         strain_key = NULL))
})

test_that("plot_continuous_phenotype gives error for invalid inputs", {
  ntip <- 10
  temp_tree <- ape::rtree(ntip)
  temp_vec <- rnorm(ntip, mean = 0, sd = 10)
  names(temp_vec) <- temp_tree$tip.label
  temp_anc_rec <-
    suppressWarnings(ape::ace(x = temp_vec,
                              phy = temp_tree,
                              type = "continuous",
                              method = "REML"))
  temp_anc_rec_at_node <- temp_anc_rec$ace
  temp_vec <- temp_vec[1:5]
  temp_tr_type <- "phylogram"

  expect_error(plot_continuous_phenotype(temp_tree,
                                         temp_vec,
                                         temp_anc_rec_at_node,
                                         temp_tr_type,
                                         strain_key = NULL))
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
  temp_trans$num_genotypes <- 1
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
  temp_trans$num_genotypes <- 1
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
  temp_trans$num_genotypes <- 1
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
  temp_trans$num_genotypes <- 1
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
  temp_trans$num_genotypes <- 2
  expect_error(hist_abs_hi_conf_delta_pheno(temp_trans,
                                            temp_tree,
                                            temp_index,
                                            temp_non_trans_color,
                                            temp_trans_color))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE),
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     NULL,
                                     FALSE))
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
  temp_tr_type <- "phylogram"

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
                                     temp_legend_other,
                                     temp_tr_type,
                                     strain_key = NULL,
                                     tip_name_log = TRUE),
               NA)


  temp_strain_key <- matrix(c("one", "one", "two", "two", "three"),
                            nrow = 5,
                            ncol = 1)
  row.names(temp_strain_key) <- temp_tree$tip.label
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
                                     temp_legend_other,
                                     temp_tr_type,
                                     strain_key = temp_strain_key,
                                     tip_name_log = TRUE),
               NA)
})

# blank_plot
test_that("blank_plot gives no errors", {
  expect_error(blank_plot(), NA)
})
