library(microbeGWAS)
context("Transition edges") #----------------------------------------------#
# TODO write tests for all functions

# TODO test identify_transition_edges

test_that("identify_transition_edges does X given Y", {


})

# identify_transition_edges <- function(tr, mat, num, node_recon, disc_cont){

#   # FUNCTION ------------------------------------------------------------------#
#   transition <- transition_direction <- parent_node <- child_node <- integer(Nedge(tr)) # initialize all as zeroes
#   older <- 1 # older node is 1st column in tr$edge
#   younger <- 2 # younger node is 2nd column in tr$edge
#   parent_0_child_1 <- 1
#   parent_1_child_0 <- -1
#   parent_equals_child <- 0
#   both_parent_and_child_are_one <- 2
#
#
#   for (i in 1:Nedge(tr)){
#     if (is_tip(tr$edge[i, older], tr)){
#       stop("tree invalid")
#     }
#     parent_node[i] <- node_recon[tr$edge[i, older] - Ntip(tr)]   # Assign node value
#     if (is_tip(tr$edge[i, younger], tr)){ # child is a tip
#       child_node[i]  <- mat[ , num][tr$edge[i, younger]]           # Assign tip value
#     } else {  # child is internal nodes
#       child_node[i]  <- node_recon[tr$edge[i, younger] - Ntip(tr)] # Assign node value
#     }
#
#     transition[i] <- sum(parent_node[i] + child_node[i])
#     # transition[i] is either 0, 1, or 2 for discrete traits
#     # transition[i] is not to be used when the trait is continuous because all, or very nearly all edges are transition edges.
#
#     if (parent_node[i] > child_node[i]){
#       transition_direction[i] <- parent_1_child_0
#     } else if (parent_node[i] < child_node[i]){
#       transition_direction[i] <- parent_0_child_1
#     }
#
#     if (disc_cont == "discrete"){
#       transition[transition == both_parent_and_child_are_one] <- parent_equals_child # If parent_node and child_node had same value (1) then no transition occured.
#     } else if (disc_cont == "continuous"){
#       transition <- NA # transition[i] is not to be used when the trait is continuous because all, or very nearly all edges are transition edges.
#     } else {
#       stop("disc_cont must be discrete or continuous")
#     }
#
#   }
#   # Check and return output ----------------------------------------------------
#   results <- list("transition" = transition, "trans_dir" = transition_direction)
#   return(results)
# } # end identify_transition_edges()
#
# TODO test keep_at_least_two_high_conf_trans_edges
# keep_at_least_two_high_conf_trans_edges <- function(genotype_transition, genotype_confidence){
#   # TODO
#   # 1) Update description of inputs
#   # 2) check inputs.
#   # 3) Update description of output.
#   # 4) Check output.
#   #
#   # Function description -------------------------------------------------------
#   # Since we're looking for convergence of transitions we need a second quality
#   # control step where we remove genotypes that have only 1 transition-edge or
#   # where the transition edges are identical!
#   #
#   # Inputs:
#   # genotype_transition.
#   # genotype_confidence.
#   #
#   # Output:
#   # has_at_least_two_high_confidence_transition_edges.
#   #
#   # Check inputs ---------------------------------------------------------------
#
#   # Function -------------------------------------------------------------------
#   has_at_least_two_high_confidence_transition_edges <- rep(FALSE, length(genotype_transition))
#   for (p in 1:length(genotype_transition)){
#     if (sum(genotype_transition[[p]]$transition * genotype_confidence[[p]]) > 1){
#       has_at_least_two_high_confidence_transition_edges[p] <- TRUE
#     }
#   }
#
#   # Check and return output ----------------------------------------------------
#   return(has_at_least_two_high_confidence_transition_edges)
# } # end keep_at_least_two_high_conf_trans_edges()
#
# TODO test keep_hits_with_more_change_on_trans_edges
# keep_hits_with_more_change_on_trans_edges <- function(results, pvals, a){
#   # TODO
#   # 1) Update description of inputs
#   # 2) check inputs.
#   # 3) Update description of output.
#   # 4) Check output.
#   #
#   # Function description -------------------------------------------------------
#   # SUBSET SIGNIFICANT HITS WHERE THE MEDIAN(DELTA PHENOTYPE) ON TRANSITION EDGES IS > MEDIAN(DELTA PHENOTYPE) ON ALL EDGES
#   #
#   # Inputs:
#   # results.
#   # pvals.
#   # a.      Number. Alpha (significance threshold).
#   # Output:
#   # has_at_least_two_high_confidence_transition_edges.
#   #
#   # Check inputs ---------------------------------------------------------------
#   check_if_alpha_valid(a)
#
#   # Function -------------------------------------------------------------------
#   temp <- pvals$hit_pvals[(results$trans_median > results$all_edges_median), , drop = FALSE]
#   hits <- temp[temp[ , 1] < a, , drop = FALSE]
#
#   # Check and return output ----------------------------------------------------
#   return(hits)
# } # end keep_hits_with_more_change_on_trans_edges()
