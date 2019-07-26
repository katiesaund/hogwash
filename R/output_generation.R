#' create_contingency_table
#'
#' @description This function generates a contingency table for the genotype
#'  and phenotype. Each entry in the table correpsonds to a tree edge. For
#'  synchronous it counts the number of tree edges that are +/- genotype
#'  transitions and +/- phenotype transitions. For phyc it counts the number of
#'  tree edges that are +/- genotype transitions and +/- phenotype presence.
#'
#' @param genotype_by_edges List.
#' @param phenotype_by_edges Numeric vector.
#' @param geno Matrix.
#' @param test_type Character. Either "synchronous" or "phyc."
#'
#' @return all_tables. List of matrices.
#'
#' @noRd
#'
create_contingency_table <- function(genotype_by_edges,
                                     phenotype_by_edges,
                                     geno,
                                     test_type){
  # Check input ----------------------------------------------------------------
  check_str_is_test_name(test_type)
  check_equal(length(genotype_by_edges), ncol(geno))
  check_equal(length(phenotype_by_edges), length(genotype_by_edges[[1]]))
  check_if_binary_vector(phenotype_by_edges)
  check_if_binary_vector(genotype_by_edges[[1]])

  # Function -------------------------------------------------------------------
  all_tables <- rep(list(NULL), length(genotype_by_edges))
  contingency_table <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)

  if (test_type == "synchronous") {
    row.names(contingency_table) <- c("geno_transition", "geno_not_trans")
    colnames(contingency_table) <- c("pheno_transition", "pheno_not_trans")
  } else if (test_type == "phyc") {
    row.names(contingency_table) <- c("geno_transition", "geno_not_trans")
    colnames(contingency_table) <- c("pheno_present", "pheno_absent")
  } else {
    stop("Test must be either synchronous or phyc")
  }

  for (i in 1:length(genotype_by_edges)) {
    temp_table <- contingency_table
    temp_table[1, 1] <-
      sum(genotype_by_edges[[i]]  + phenotype_by_edges == 2,  na.rm = TRUE)
    temp_table[2, 2] <-
      sum(genotype_by_edges[[i]]  + phenotype_by_edges == 0,  na.rm = TRUE)
    temp_table[1, 2] <-
      sum(-genotype_by_edges[[i]] + phenotype_by_edges == -1, na.rm = TRUE)
    temp_table[2, 1] <-
      sum(-genotype_by_edges[[i]] + phenotype_by_edges == 1,  na.rm = TRUE)
    all_tables[[i]] <- temp_table
  }
  names(all_tables) <- colnames(geno)

  # Return output --------------------------------------------------------------
  return(all_tables)
} # end create_contigency_table()
