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
  if (!typeof(mat) %in% c("double", "integer", "numeric")){stop("Must be a numeric matrix")}

  # Function -------------------------------------------------------------------
  type <- "discrete"
  if (sum(!(mat %in% c(0, 1))) > 0){
    type <- "continuous"
  }

  # Check and return output ----------------------------------------------------
  check_is_string(type)
  check_str_is_discrete_or_continuous(type)

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
  # Given a continuous phenotype, check if the phenotype follows a normal
  # distribution. The brownian motional model of ancestral reconstruction, which
  # this library uses, assumes a normal distribution of the phenotype.
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # pheno. Matrix. A one column numeric matrix.
  # Continuous_or_discrete. String. Either 'continuous' or 'discrete.'
  #
  # Outputs:
  # None.
  #
  # Check input ----------------------------------------------------------------
  check_dimensions(pheno, NULL, 1, 1, 1)
  if (!typeof(pheno) %in% c("double", "integer", "numeric")){stop("Must be a numeric matrix")}
  check_is_string(continuous_or_discrete)
  check_str_is_discrete_or_continuous(continuous_or_discrete)

  # Function -------------------------------------------------------------------
  if (continuous_or_discrete == "continuous"){
    if (length(unique(pheno)) == 1){stop("phenotype must have some variability")}

    result <- shapiro.test(unlist(pheno))
    alpha <- 0.05
    if (result$p < alpha){
      print("Please consider normalizing your phenotype so as to not violate
            assumptions used in the ancestral reconstruction.")
    }
  }
} # end check_if_phenotype_normal

check_if_convergence_occurs <- function(pheno, tr, continuous_or_discrete){
  # Function description -------------------------------------------------------
  # Given the phenotype and tr, test if a Brownian motion model or a white
  # noise model better fit the data. If a white noise model better fits the data
  # it could be because convergent or parallel evolution in happening on
  # disparate parts of the tr It's just one more clue.
  #
  # Inputs:
  # pheno. Numeric matrix. One column.
  # tr. Phylo.
  # continuous_or_discrete. String. Either 'continuous' or 'discrete'
  #
  # Outputs:
  # TODO? Plot?
  #
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)
  check_dimensions(pheno, NULL, 1, 1, 1)
  check_str_is_discrete_or_continuous(continuous_or_discrete)

  # Function -------------------------------------------------------------------
  if (continuous_or_discrete == "continuous"){
    set.seed(1)
    geiger_BM <- geiger::fitContinuous(tr, pheno, model = "BM")
    geiger_white <- geiger::fitContinuous(tr, pheno, model = "white")

    enough_of_a_difference_in_AIC <- 2
    if (geiger_white$opt$aicc + enough_of_a_difference_in_AIC < geiger_BM$opt$aicc){
      print("A white noise model better describes phenotype than does a Brownina motion model.")
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
  # Return output --------------------------------------------------------------
  # TODO? Add the plot? If I do, then I need to add unit tests
  # TODO where is best place to put this function in the large script?
} # end check_if_convergence_occurs()

