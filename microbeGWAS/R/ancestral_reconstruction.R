# All of the functions necessary to perform the ancestral reconstruction of traits.

# General note on ape::ace() setting:
# Marginal = FALSE means that the marginal is in fact calculated. Do not set
# marginal to TRUE, because this only does the condition (only downstream edges
# considered). ace never calculates the joint.

pick_recon_model <- function(mat, tr, disc_cont, num, recon_method){
  # Function description -------------------------------------------------------
  # Given a discrete genotype or phenotype, build an ancestral reconstruction
  # using an equal rates model ("ER") and then repeat with an all rates
  # different model ("ARD"). Choose the best model based on p-value from
  # likelihood ratio test and difference in AIC.
  #
  # Inputs:
  # mat. Matrix. Genotype or phenotype binary matrix. (Only discrete phenotype)
  # tr. Phylo.
  # num. Number. Column index for matrix.
  # disc_cont. String. "discrete" or "continuous".
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  check_if_binary_matrix(mat)
  check_is_string(disc_cont)
  check_is_number(num)
  check_is_string(recon_method)
  check_if_ancestral_reconstruction_method_compatible_with_ape(recon_method)
  if (disc_cont != "discrete"){
    stop("Only pick recon model for discrete characters. Continuous characters must be BM.")
  }
  if (num > ncol(mat)){
    stop("Index must be 1 <= index <= ncol(matrix)")
  }

  # Function ------------------------------------------------------------------#
  # Note, SYMreconstruction removed because SYM == ER for binary inputs.
  # Use this function to choose the best model for reconstruction.

  # Cutoffs for comparing the ER and ARD:
  alpha <- 0.05 # For likelihood test
  significant_difference_in_AIC <- 2

  # Reference for model testing: https://www.r-phylo.org/wiki/HowTo/Ancestral_State_Reconstruction & http://blog.phytools.org/2017/07/comparing-fitted-discrete-character.html
  # Test ER vs ARD
  set.seed(1)
  ERreconstruction  <- ace(mat[ , num, drop = TRUE], tr, type = disc_cont, method = recon_method, marginal = FALSE, model = "ER")

  # Some ARD models don't work well with the data and given a warning message like:  "In sqrt(diag(solve(h))) : NaNs produced"
  # To ensure the ER model is prefered in this case use the following warning catching:
  error_msg <- "ARD_bad_fit"
  set.seed(1)
  ARDreconstruction <- tryCatch(ace(mat[ , num, drop = TRUE], tr, type = disc_cont, method = recon_method, marginal = FALSE, model = "ARD"), warning = function(x) {error_msg})
  # If ARD gave a warning, pick ER
  best_model <- "ER"
  if (length(ARDreconstruction) == 1){
    best_model <- "ER"
  } else { # Pick best of ER and ARD
    p_ER_ARD  <- 1 - pchisq(2*abs(ERreconstruction$loglik - ARDreconstruction$loglik), 1)
    if (p_ER_ARD < alpha & AIC(ERreconstruction) > (significant_difference_in_AIC + AIC(ARDreconstruction))){
      best_model <- "ARD"
    }
  }
  return(best_model)
} # end pick_recon_model()

continuous_ancestral_reconstruction <- function(tr, mat, num, disc_cont, recon_method){
  # Function Description -------------------------------------------------------
  # Perform ancestral state reconstruction using ape::ace() function on
  # continuous input data.
  #
  # Inputs:
  # mat. Matrix. Genotype or phenotype binary matrix. (Only discrete phenotype)
  # tr. Phylo.
  # num. Number. Column index for matrix.
  # disc_cont. String. "discrete" or "continuous".
  # recon_method. String. "ML", "REML", "pic", or "GLS."
  #
  # Outputs:
  # "ML_anc_rec" Vector. Reconstruction (numeric) at node values only.
  # "tip_and_node_recon" = tip_and_node_recon. A vector with the tip values followed by the node ancestral reconstruction.
  #
  # Check inputs ---------------------------------------------------------------
  check_if_ancestral_reconstruction_method_compatible_with_ape(recon_method)
  check_is_string(disc_cont)
  check_str_is_discrete_or_continuous(disc_cont)
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  check_is_number(num)
  if (class(mat) != "matrix"){stop("Ancestral reconstruction input must be a matrix")}

  # Function -------------------------------------------------------------------
  set.seed(1)
  reconstruction <- ace(mat[ , num, drop = TRUE],
                        tr,
                        model = "BM",
                        type = disc_cont,
                        method = recon_method,
                        marginal = FALSE)
  ML_anc_rec <- reconstruction$ace # vector containing reconstruction data for all internal nodes (#51 -99) where tips are 1-50.
  tip_and_node_recon <- c(mat[ , num, drop = TRUE], ML_anc_rec)
  names(tip_and_node_recon) <- c(1:sum(Ntip(tr), Nnode(tr)))

  # Return output --------------------------------------------------------------
  return(list("ML_anc_rec" = ML_anc_rec,
              "tip_and_node_recon" = tip_and_node_recon))
} # end continuous_ancestral_reconstruction()

continuous_get_recon_confidence <- function(recon_vector){
  # Function description -------------------------------------------------------
  # Given a vector of a continuous trait reconstruction, generate a dummy
  # onfidence vector. Ancestral reconstruction for continuous values only gives
  # back a 95% CI. We can't use any of this information to decide which nodes
  # are low confidence so treat all reconstructed values as high confidence,
  # which we stores as a value 1.
  #
  # Input:
  # recon_vector. Vector. Numeric. Length of Ntip(tr) + Nnode(tr).
  #
  # Output:
  # tip_and_node_anc_rec_confidence. Vector of 1s; length of Ntip(tr) + Nnode(tr).
  #
  # Check input ----------------------------------------------------------------
  if(!is.vector(recon_vector)){stop("input must be a vector")}

  # Function -------------------------------------------------------------------
  tip_and_node_anc_rec_confidence <- rep(1, length(recon_vector))

  # Return output --------------------------------------------------------------
  return(tip_and_node_anc_rec_confidence)
} # end continuous_get_recon_confidence()

discrete_ancestral_reconstruction <- function(tr, mat, num, disc_cont, recon_method){
  # Function Description -------------------------------------------------------
  # Perform ancestral state reconstruction using ape::ace() function on
  # continuous input data.
  #
  # Inputs:
  # mat. Matrix. Genotype or phenotype binary matrix. (Only discrete phenotype)
  # tr. Phylo.
  # num. Number. Column index for matrix.
  # disc_cont. String. "discrete".
  # recon_method. String. "ML", "REML", "pic", or "GLS."
  #
  # Outputs:
  # "ML_anc_rec" Vector. Reconstruction (numeric) at node values only.
  # "tip_and_node_recon" = tip_and_node_recon. A vector with the tip values followed by the node ancestral reconstruction.
  # "reconstruction" = Ape::ace() object. Contains multiple piecs of information.

  # Check inputs ---------------------------------------------------------------
  check_if_ancestral_reconstruction_method_compatible_with_ape(recon_method)
  check_is_string(disc_cont)
  check_str_is_discrete_or_continuous(disc_cont)
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  check_is_number(num)
  if (class(mat) != "matrix"){stop("Ancestral reconstruction input must be a matrix")}
  if(disc_cont != "discrete"){stop("Discrete ancestral reconstruction only")}

  # Function -------------------------------------------------------------------
  set.seed(1)
  recon_model <- pick_recon_model(mat, tr, disc_cont, num, recon_method)
  reconstruction <- ace(mat[ , num, drop = TRUE],
                        tr,
                        model = recon_model,
                        type = disc_cont,
                        method = recon_method,
                        marginal = FALSE)
  ML_anc_rec <- as.numeric(colnames(reconstruction$lik.anc)[apply(reconstruction$lik.anc, 1, which.max)]) # Extract the mostly likely character state using which.max
  names(ML_anc_rec) <- c((Ntip(tr) + 1):(Ntip(tr) + Nnode(tr)))
  tip_and_node_recon <- c(mat[ , num, drop = TRUE], ML_anc_rec)
  names(tip_and_node_recon) <- c(1:sum(Ntip(tr), Nnode(tr)))

  # Return outputs -------------------------------------------------------------
  return(list("tip_and_node_recon" = tip_and_node_recon,
              "ML_anc_rec" = ML_anc_rec,
              "reconstruction" = reconstruction))
} # end discrete_ancestral_reconstruction()

discrete_get_recon_confidence <- function(recon, tr, ML_cutoff){
  # Function description -------------------------------------------------------
  # Compute ancestral reconstruction from a continuous or discrete input (genotype or phenotype matrix).
  #
  # Inputs:
  # tre:       Phylo.     Rooted phylogenetic tree.
  # ML_cutoff: Numeric.   Indicating the row of the matrix from which the ancestral reconstruction is to be built.
  # recon: Ape reconstruction object
  #
  # Outputs:
  # "tip_and_node_rec_conf" = tip_and_node_anc_rec_confidence. A vector with a 1 for each tip with the node ancestral reconstruction confidence appended.

  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  check_if_alpha_valid(ML_cutoff)
  # Function -------------------------------------------------------------------

  anc_rec_confidence <- apply(recon$lik.anc, 1, max) # Get the highest confidence value at each node
  tip_and_node_anc_rec_confidence <- c(rep(1, Ntip(tr)), anc_rec_confidence) # count all tips as high confidence
  tip_and_node_anc_rec_confidence <- discretize_confidence_using_threshold(tip_and_node_anc_rec_confidence, ML_cutoff) # Count anything lower than threshold as low confidence

  # Return output --------------------------------------------------------------
  return(tip_and_node_anc_rec_confidence)
} # end discrete_get_recon_confidence()


ancestral_reconstruction_by_ML <- function(tr, mat, num, disc_cont){
  # Function description -------------------------------------------------------
  # Compute ancestral reconstruction from a continuous or discrete input (genotype or phenotype matrix).
  #
  # Inputs:
  # tree:      Phylo.     Rooted phylogenetic tree.
  # mat:       Matrix.    Either the phenotype matrix or the genotype matrix.
  # num:       Numeric.   Indicating the row of the matrix from which the ancestral reconstruction is to be built.
  # disc_cont: Character. Either "discrete" or "continuous".
  #
  # Outputs:
  # "node_anc_rec" = ML_anc_rec. Vector. Ancestral reconstruction for each node.
  # "tip_and_node_rec_conf" = tip_and_node_anc_rec_confidence. A vector with a 1 for each tip with the node ancestral reconstruction confidence appended.
  # "recon_edge_mat" = reconstruction_as_edge_mat. Matrix. Ancestral reconstruction configured into the edge matrix.
  # "tip_and_node_recon" = tip_and_node_recon. A vector with the tip values followed by the node ancestral reconstruction.
  #
  # Check input ----------------------------------------------------------------
  check_is_string(disc_cont)
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  check_is_number(num)
  if (class(mat) != "matrix"){stop("Ancestral reconstruction input must be a matrix")}
  check_str_is_discrete_or_continuous(disc_cont)
  # Function -------------------------------------------------------------------
  # Compute ancestral reconstruction
  recon_method <- "ML" # ML == Maximum Likelihood.
  ML_significance_threshold <- .875 # ML literature suggests that a ratio of 7:1 suggests a high confidence ancestral reconstruction per node .875/.125 = 7

  if (disc_cont == "continuous"){
    cont_results <- continuous_ancestral_reconstruction(tr, mat, num, disc_cont, recon_method) # RECONSTRUCTION
    ML_anc_rec <- cont_results$ML_anc_rec
    tip_and_node_recon <- cont_results$tip_and_node_recon
    tip_and_node_anc_rec_confidence <- continuous_get_recon_confidence(tip_and_node_recon) # CONFIDENCE IN RECONSTRUCTION

  } else if (disc_cont == "discrete"){
    discrete_results <- discrete_ancestral_reconstruction(tr, mat, num, disc_cont, recon_method) # RECONSTRUCTION
    ML_anc_rec <- discrete_results$ML_anc_rec
    tip_and_node_recon <- discrete_results$tip_and_node_recon
    tip_and_node_anc_rec_confidence <- discrete_get_recon_confidence(discrete_results$reconstruction, tr, ML_significance_threshold) # CONFIDENCE IN RECONSTRUCTION
  }
  reconstruction_as_edge_mat <- convert_to_edge_mat(tr, tip_and_node_recon)

  # Check and return outputs ---------------------------------------------------
  results <- list("node_anc_rec" = ML_anc_rec,
                  "tip_and_node_rec_conf" = tip_and_node_anc_rec_confidence,
                  "recon_edge_mat" = reconstruction_as_edge_mat,
                  "tip_and_node_recon" = tip_and_node_recon)
  if (class(reconstruction_as_edge_mat) != "matrix"){stop("reconstruction_as_edge_mat must be a matrix")}
  if (ncol(reconstruction_as_edge_mat) != ncol(tr$edge)){stop("reconstruction_as_edge_mat must have same dim as tr$edge")}
  if (nrow(reconstruction_as_edge_mat) != nrow(tr$edge)){stop("reconstruction_as_edge_mat must have same dim as tr$edge")}
  if (Nnode(tr) != length(ML_anc_rec)){stop("Ancestral reconstruction of nodes is the wrong length.")}
  if (Nnode(tr) + Ntip(tr) != length(tip_and_node_anc_rec_confidence)){stop("Ancestral reconstruction confidence is the wrong length.")}
  if (Nnode(tr) + Ntip(tr) != length(tip_and_node_recon)){stop("Ancestral reconstruction is the wrong length.")}
  return(results)
} # end ancestral_reconstruction_by_ML()

convert_to_edge_mat <- function(tr, tip_and_node_reconstruction){
  # Function description -------------------------------------------------------
  # Convert the reconstruction to be in the same format at tree$edge, where each
  # row is an edge. This format will make it much easier to calculate the
  # phenotype change on each edge for the continuous phenotype GWAS.
  #
  # Inputs:
  # tr. Phylo.
  # tip_and_node_reconstruction. Vector. Ordered by tips then by nodes.
  #
  # Outputs:
  # reconstruction_as_edge_mat. Matrix.
  #
  # Check inputs ---------------------------------------------------------------
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)
  if (length(tip_and_node_reconstruction) != Nnode(tr) + Ntip(tr)){
    stop("Reconstruction incorrectly formatted.")
  }

  # Function -------------------------------------------------------------------
  reconstruction_as_edge_mat <- tr$edge
  for (k in 1:nrow(reconstruction_as_edge_mat)){
    reconstruction_as_edge_mat[k, 1] <- tip_and_node_reconstruction[tr$edge[k, 1]]
    reconstruction_as_edge_mat[k, 2] <- tip_and_node_reconstruction[tr$edge[k, 2]]
  }

  # Return output --------------------------------------------------------------
  return(reconstruction_as_edge_mat)
} # end convert_to_edge_mat()
