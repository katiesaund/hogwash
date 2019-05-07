# FUNCTIONS FOR PHYC ----------------------------------------------------------#
create_contingency_table <- function(genotype_by_edges, phenotype_by_edges, geno){
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
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------

  all_tables <- rep(list(NULL), length(genotype_by_edges))
  contingency_table <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
  row.names(contingency_table) <- c("geno_present", "geno_absent")
  colnames(contingency_table) <- c("pheno_present", "pheno_absent")
  for (i in 1:length(genotype_by_edges)){
    temp_table <- contingency_table
    temp_table[1, 1] <- sum(genotype_by_edges[[i]]  + phenotype_by_edges == 2,  na.rm = TRUE)
    temp_table[2, 2] <- sum(genotype_by_edges[[i]]  + phenotype_by_edges == 0,  na.rm = TRUE)
    temp_table[1, 2] <- sum(-genotype_by_edges[[i]] + phenotype_by_edges == -1, na.rm = TRUE)
    temp_table[2, 1] <- sum(-genotype_by_edges[[i]] + phenotype_by_edges == 1,  na.rm = TRUE)
    all_tables[[i]] <- temp_table
  }

  names(all_tables) <- colnames(geno)
  return(all_tables)
} # end create_contigency_table()


create_file_name <- function(output_dir, output_name, other_info){
  # Function description -------------------------------------------------------
  # Create a file name string that includes the path to the output directory.
  #
  # Input:
  # output_dir.  Character.
  # output_name. Character.
  # other_info.  Character.
  #
  # Output:
  # file_name. Character.
  #
  # Check input ----------------------------------------------------------------
  check_if_dir_exists(output_dir)
  check_is_string(output_name)
  check_is_string(other_info)

  # Function -------------------------------------------------------------------
  file_name <- paste(output_dir, "/", "phyc_", output_name, "_", other_info, sep = "")

  # Check and return output ----------------------------------------------------
  check_is_string(file_name)
  return(file_name)
} # end create_file_name()
