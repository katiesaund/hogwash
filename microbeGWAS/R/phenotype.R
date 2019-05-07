assign_pheno_type <- function(mat){
  # Function description -------------------------------------------------------
  # Determine if the matrix is discrete or continuous.
  #
  # Input:
  # mat: Matrix. Should be a phenotype with 1 column.
  #
  # Outputs:
  # type. Character. Either "discrete" or "continuous".

  # Check input ----------------------------------------------------------------
  check_dimensions(mat, NULL, 1, 1, 1)

  # Function -------------------------------------------------------------------
  type <- "discrete"
  if (sum(!(mat %in% c(0, 1))) > 0){
    type <- "continuous"
  }

  # Check and return output ----------------------------------------------------
  check_is_string(type)
  return(type)
} # end assign_pheno_type()

calculate_phenotype_change_on_edge <- function(edge_list, phenotype_by_edges){
  # Function description -------------------------------------------------------
  # Quantify absoluate value of phenotype change on each tree edge.
  #
  # Input:
  # edge_list.          Numeric vector. Each number is the index of the tree edge to be used
  # phenotype_by_edges. Mat.            Dimensions: Nedge x 2 matrix. Entries are the phenotype value at the node, where row is the edge, 1st column is the parent node and 2nd column is the child node.
  #
  # Outputs:
  # delta.              Numeric vector.   Length = length(edge_list).

  # Check input ----------------------------------------------------------------
  if (max(edge_list) > nrow(phenotype_by_edges)){
    stop("Cannot calculate phenotype change on edges.")
  }
  check_dimensions(phenotype_by_edges, exact_rows = NULL, min_rows = max(edge_list), exact_cols = NULL, min_cols = 2)

  # Function -------------------------------------------------------------------
  delta <- rep(NA, length(unique(edge_list)))
  for (j in 1:length(edge_list)){
    delta[j] <- abs(phenotype_by_edges[edge_list[j], 1] - phenotype_by_edges[edge_list[j], 2])
  }

  # Check and return output ----------------------------------------------------
  return(delta)
} # end calculate_phenotype_change_on_edge()

calc_raw_diff <- function(edge_list, ph_edges){
  # Function description ------------------------------------------------------#
  # Subtract child node phenotype value from parent node phenotype value.
  #
  # Inputs:
  # TODO
  # edge_list. ?
  # ph_edges. ?
  #
  # Output:
  # delta. ?
  #
  # Check inputs ---------------------------------------------------------------
  # TODO
  # Function -------------------------------------------------------------------
  delta <- rep(NA, length(unique(edge_list)))
  for (j in 1:length(edge_list)){
    delta[j] <- ph_edges[edge_list[j], 1] - ph_edges[edge_list[j], 2]
  }
  # Return outputs -------------------------------------------------------------
  return(delta)
} # end calc_raw_diff()


check_if_phenotype_normal <- function(pheno, continuous_or_discrete){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  if (continuous_or_discrete == "continuous"){
    result <- shapiro.test(unlist(pheno))
    alpha <- 0.05
    if (result$p < alpha){
      print("Consider normalizing your phenotype")
    }
  }
} # end check_if_phenotype_normal

check_if_convergence_occurs <- function(pheno, tree, continuous_or_discrete){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  if (continuous_or_discrete == "continuous"){
    set.seed(1)
    geiger_BM <- fitContinuous(tree, pheno, model = "BM")
    geiger_white <- fitContinuous(tree, pheno, model = "white")

    if (geiger_white$opt$aicc < geiger_BM$opt$aicc){
      print("WN better than BM")
    }

    # TODO Add this as a plot to output?
    #pdf(paste(test_dir, "/evolutionary_model_barplot.pdf", sep =""))
    #par(mfrow = c(1,1))
    #barplot(c(geiger_BM$opt$aicc, geiger_OU$opt$aicc, geiger_white$opt$aicc),
    #        names.arg = c("Brownian Motion", "Ornstein-Uhlenbeck", "White Noise"),
    #        col = c("black", "black", "black"),ylim = c(0, 600),
    #        ylab = "AICc", xlab = "Evolutionary model",
    #        main = "Model values for toxin")
    #dev.off()
  }
} # end check_if_convergence_occurs()

