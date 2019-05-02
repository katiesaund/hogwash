pick_recon_model <- function(mat, tr, disc_cont, num, recon_method){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # mat. Matrix. Genotypoe or phenotype binary matrix. (Only discrete phenotype)
  # tr. Phylo.
  # num. Number. index for matrix.
  # disc_cont. String. "discrete" or "continuous".

  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------

  # Note, SYMreconstruction removed because SYM == ER for binary inputs.
  # Use this function to choose the best model for reconstruction.

  # Check inputs --------------------------------------------------------------#
  check_is_string(disc_cont)
  if (disc_cont != "discrete"){
    stop("Only pick recon model for discrete characters. Continuous characters must be BM.")
  }

  check_if_binary_matrix(mat)

  # Function ------------------------------------------------------------------#
  alpha <- 0.05
  significant_difference_in_AIC <- 2

  # Reference for model testing: https://www.r-phylo.org/wiki/HowTo/Ancestral_State_Reconstruction & http://blog.phytools.org/2017/07/comparing-fitted-discrete-character.html
  # Test ER vs ARD
  ERreconstruction  <- ace(mat[ , num, drop = TRUE], tr, type = disc_cont, method = recon_method, marginal = FALSE, model = "ER")

  # Some ARD models don't work well with the data and given a warning message like:  "In sqrt(diag(solve(h))) : NaNs produced"
  # To ensure the ER model is prefered in this case use the following warning catching:
  error_msg <- "ARD_bad_fit"
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
    kER  <- length(ERreconstruction$rates)
    kARD <- length(ARDreconstruction$rates)
    new_p_ER_ARD  <- pchisq(-2*(logLik(ERreconstruction)-logLik(ARDreconstruction)),  df = kARD - kER,  lower.tail = FALSE)

    num_digits <- 4
    if (round(p_ER_ARD, num_digits) != round(new_p_ER_ARD, num_digits)){
      stop("ER_ARD loglikelihood test should be the same from both calculations")
    }
  }

  return(best_model)
} # end pick_recon_model()


ancestral_reconstruction_by_ML <- function(tr, mat, num, disc_cont){
  # TODO:
  # 1) Break this function into two subfuctions: one for continuous, one for discrete
  # 2) For those two subfunctions:
  #    a) subfunctionalize
  #    b) comment
  #    c) add output descriptions
  #    d) finish checking all inputs are correctly formatted.

  # Function description -------------------------------------------------------
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # tree:      Phylo.     Rooted phylogenetic tree.
  # mat:       Matrix.    Either the phenotype matrix or the genotype matrix.
  # num:       Numeric.   Indicating the row of the matrix from which the ancestral reconstruction is to be built.
  # disc_cont: Character. Either "discrete" or "continuous".
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec
  # "anc_rec_confidence" = tip_and_node_anc_rec_confidence
  #
  # Check input ----------------------------------------------------------------
  check_is_string(disc_cont)
  check_for_root_and_bootstrap(tr)

  # Function -------------------------------------------------------------------

  # Compute ancestral reconstruction
  recon_method <- "ML" # Method: ML.         Maximum Likelihood.
  ML_significance_threshold <- .875 # ML literature suggests that a ratio of 7:1 suggests a high confidence ancestral reconstruction per node .875/.125 = 7

  if (disc_cont == "continuous"){

    # RECONSTRUCTION -----------------------------------------------------------
    set.seed(3)
    reconstruction <- ace(mat[ , num, drop = TRUE], tr, type = disc_cont, method = recon_method, model = "BM", marginal = FALSE)
    ML_anc_rec <- reconstruction$ace # vector containing reconstruction data for all internal nodes (#51 -99) where tips are 1-50.
    tip_and_node_recon <- c(mat[ , num, drop = TRUE], ML_anc_rec)
    names(tip_and_node_recon) <- c(1:sum(Ntip(tr), Nnode(tr)))

    # CONFIDENCE IN RECONSTRUCTION ---------------------------------------------
    tip_and_node_anc_rec_confidence <- rep(1, length(tip_and_node_recon)) # Ancestral reconstruction for continuous values only gives back a 95% CI. We can't use any of this information to decide which nodes are low confidence so treat all reconstructed values as high confidence.

  } else if (disc_cont == "discrete"){

    # RECONSTRUCTION -----------------------------------------------------------
    set.seed(4)
    recon_model <- pick_recon_model(mat, tr, disc_cont, num, recon_method)
    reconstruction <- ace(mat[ , num, drop = TRUE], tr,
                          model = recon_model,
                          type = disc_cont,
                          method = recon_method,
                          marginal = FALSE) # Marginal = FALSE means that the marginal is in fact calculated. Do not set marginal to TRUE, because this only does the condition (only downstream edges considered). ace never calculates the joint.
    ML_anc_rec <- as.numeric(colnames(reconstruction$lik.anc)[apply(reconstruction$lik.anc, 1, which.max)]) # Extract the mostly likely character state using which.max
    names(ML_anc_rec) <- c((Ntip(tr) + 1):(Ntip(tr) + Nnode(tr)))
    tip_and_node_recon <- c(mat[ , num, drop = TRUE], ML_anc_rec)
    names(tip_and_node_recon) <- c(1:sum(Ntip(tr), Nnode(tr)))

    # CONFIDENCE IN RECONSTRUCTION ---------------------------------------------
    anc_rec_confidence <- apply(reconstruction$lik.anc, 1, max) # Get the highest confidence value at each node
    tip_and_node_anc_rec_confidence <- c(rep(1, Ntip(tr)), anc_rec_confidence) # count all tips as high confidence
    tip_and_node_anc_rec_confidence <- discretize_confidence_using_threshold(tip_and_node_anc_rec_confidence, ML_significance_threshold) # Count anything lower than threshold as low confidence

  } else {
    stop("Ancestral reconstruction cannot run.")
  }

  # AS A MATRIX, IN THE SAME FORMAT AS TREE$EDGE. SO NODE VALUES WILL APPEAR MULTIPLE TIMES IN THE tr.
  # THIS FORMAT WILL MAKE IT MUCH EASIER TO CALCULATE PHENOTYPE CHANGE ON EDGES.
  reconstruction_as_edge_mat <- tr$edge
  for (k in 1:nrow(reconstruction_as_edge_mat)){
    reconstruction_as_edge_mat[k, 1] <- tip_and_node_recon[tr$edge[k, 1]]
    reconstruction_as_edge_mat[k, 2] <- tip_and_node_recon[tr$edge[k, 2]]
  }


  if (Nnode(tr) != length(ML_anc_rec)){
    stop("Ancestral reconstruction of nodes is the wrong length.")
  }

  results <- list("node_anc_rec" = ML_anc_rec,
                  "tip_and_node_rec_conf" = tip_and_node_anc_rec_confidence,
                  "recon_edge_mat" = reconstruction_as_edge_mat,
                  "tip_and_node_recon" = tip_and_node_recon
  )
  return(results)
} # end ancestral_reconstruction_by_ML()
