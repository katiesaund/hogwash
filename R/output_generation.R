create_contingency_table <- function(genotype_by_edges, phenotype_by_edges, geno, test_type){
  # TODO
  # add names to contingency table results
  # add validations
  # add unit tests
  # add function description
  # genotype_by_edges should be a list
  # phenotype_by_edges should be a numeric vector
  # length of genotype_by_edges should be ncol of geno
  # length of _phenotype_by_edges should be same length as every entry in genotype_by_edges
  # all genotype_by_edges and phenotype_by_edges should be 0s or 1
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # test_type. Character string. Either "convergence" or "synchronous"
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------
  check_is_string(test_type)

  # Function -------------------------------------------------------------------



  all_tables <- rep(list(NULL), length(genotype_by_edges))
  contingency_table <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)

  if (test_type == "synchronous") {
    row.names(contingency_table) <- c("geno_transition", "geno_not_trans")
    colnames(contingency_table) <- c("pheno_transition", "pheno_not_trans")
  } else if (test_type == "convergence") {
    row.names(contingency_table) <- c("geno_transition", "geno_not_trans")
    colnames(contingency_table) <- c("pheno_present", "pheno_absent")
  } else {
    stop("Test must be either synchronous or convergence")
  }


  for (i in 1:length(genotype_by_edges)) {
    temp_table <- contingency_table
    temp_table[1, 1] <- sum(genotype_by_edges[[i]]  + phenotype_by_edges == 2,  na.rm = TRUE)
    temp_table[2, 2] <- sum(genotype_by_edges[[i]]  + phenotype_by_edges == 0,  na.rm = TRUE)
    temp_table[1, 2] <- sum(-genotype_by_edges[[i]] + phenotype_by_edges == -1, na.rm = TRUE)
    temp_table[2, 1] <- sum(-genotype_by_edges[[i]] + phenotype_by_edges == 1,  na.rm = TRUE)
    all_tables[[i]] <- temp_table
  }

  names(all_tables) <- colnames(geno)

  # Return output --------------------------------------------------------------
  return(all_tables)
} # end create_contigency_table()
