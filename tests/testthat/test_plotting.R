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
