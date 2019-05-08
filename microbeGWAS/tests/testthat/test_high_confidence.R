library(microbeGWAS)
context("High confidence") #----------------------------------------------#
# TODO write unit tests for all three functions

# test discretize_confidence_using_threshold
test_that("discretize_confidence_using_threshold should give this expected result", {
  temp_vector <- c(1:100)
  temp_threshold <- 50
  expect_equal(discretize_confidence_using_threshold(temp_vector, temp_threshold), c(rep(0, 49), rep(1, 51)))
})

test_that("discretize_confidence_using_threshold should give this expected result - now with fractions", {
  temp_vector <- seq(from = 0, to = 1, by = 0.1)
  temp_threshold <- 0.5
  expect_equal(discretize_confidence_using_threshold(temp_vector, temp_threshold), c(rep(0, 5), rep(1, 6)))
})


# TODO fix warning by figuring how to create better test data.
# test report_num_high_confidence_trans_edge
test_that("report_num_high_confidence_trans_edge returns expected outcome for this test set", {
  fake_geno_names <- letters[1:3]
  fake_trans <- NULL
  fake_trans <- fake_hi_conf_edges <- rep(list(0), length(fake_geno_names))
  fake_trans[[1]]$transition <- c(1, 1, 1)
  fake_trans[[2]]$transition <- c(0, 0, 0)
  fake_trans[[3]]$transition <- c(1, 1, 0)

  fake_hi_conf_edges[[1]] <- fake_hi_conf_edges[[2]] <- fake_hi_conf_edges[[3]] <- c(1, 0, 1)

  report_num_high_confidence_trans_edge(fake_trans, fake_hi_conf_edges,fake_geno_names)
  expected_result <- c(2, 0, 1)
  names(expected_result) <- fake_geno_names
  expect_equal(report_num_high_confidence_trans_edge(fake_trans, fake_hi_conf_edges,fake_geno_names), expected_result)
})


# TODO test assign_high_confidence_to_transition_edges
# assign_high_confidence_to_transition_edges <- function(tr, all_confidence_by_edge, genotype_transition_by_edges, geno){
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
#   # VALIDATION
#   if (length(genotype_transition_by_edges[[1]]$transition) != Nedge(tr)){
#     stop("Dimension mismatch")
#   }
#   check_for_root_and_bootstrap(tr)
#   if (length(all_confidence_by_edge[[1]]) != Nedge(tr)){
#     stop("Dimension mismatch")
#   }
#
#   # FUNCTION ----------------------------------------------------------------#
#
#   # Identify all edges for which the edge and the parent edge are both high confidence
#   edge_and_parent_both_confident <- edge_and_parent_confident_and_trans_edge <- rep(list(rep(0, Nedge(tr))), ncol(geno))
#   for (ge in 1:ncol(geno)){
#     for (ed in 2:(Nedge(tr) - 1)){ # start at 2 because there isn't a parent edge to edge 1, end at Nedge- 1 because no parent to the last edge either
#       parent_edge <- find_parent_edge(tr, ed)
#       if (all_confidence_by_edge[[ge]][ed] == 1 & all_confidence_by_edge[[ge]][parent_edge] == 1){
#         edge_and_parent_both_confident[[ge]][ed] <- 1
#       }
#     }
#     edge_and_parent_both_confident[[ge]][1] <- all_confidence_by_edge[[ge]][1] # have to accoutn for the fact that there isn't a parent edge to edge 1
#     edge_and_parent_both_confident[[ge]][Nedge(tr)] <- all_confidence_by_edge[[ge]][Nedge(tr)] # have to accoutn for the fact that there isn't a parent edge to last edge
#   }
#
#
#   # Identify high confidence transition edges by overlapping the above and transitions
#   for (k in 1:ncol(geno)){
#     edge_and_parent_confident_and_trans_edge[[k]] <- as.numeric((edge_and_parent_both_confident[[k]] + genotype_transition_by_edges[[k]]$transition) == 2)
#   }
#
#   # Return that overlap as high confidence transitions
#   return(edge_and_parent_confident_and_trans_edge)
# } # end assign_high_confidence_to_transition_edges()
