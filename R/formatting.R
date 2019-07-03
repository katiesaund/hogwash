format_tree <- function(tr){
  check_class(tr, "phylo")
  if (is.null(tr$node.label)) {
    stop("trees must have support values at the nodes")
  }

  if (!ape::is.rooted(tr)) {
    tr <- phytools::midpoint.root(tr)
  }
  for (i in 1:ape::Nnode(tr)) {
    if (tr$node.label[i] == "") {
      tr$node.label[i] <- 0
    }
  }
  tr$node.label <- as.numeric(tr$node.label)

  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)

  return(tr)
}

convert_matrix_to_vector <- function(mat){
  # Function description -------------------------------------------------------
  # Convert a single column matrix into a vector, retain row names as names of vector.
  #
  # Input:
  # mat. Matrix. Matrix should have only 1 column.
  #
  # Output:
  # vec. Vector.
  #
  # Check input ----------------------------------------------------------------
  check_dimensions(mat, NULL, 1, 1, 1)

  # Function -------------------------------------------------------------------
  vec <- as.vector(unlist(mat[ , 1]))
  names(vec) <- row.names(mat)

  # Check and return output ----------------------------------------------------
  check_if_vector(vec)
  return(vec)
} # end convert_matrix_to_vector()


create_test_data <- function(){
  # Function description -------------------------------------------------------
  # Create a set of reproducible test data.
  #
  # Input:
  # None.
  #
  # Output:
  # tree.             Phylo.
  # phenotype_matrix. Matrix.
  # genotype_matrix.  Matrix.
  #
  # Function -------------------------------------------------------------------
  # Create tree
  set.seed(1)
  tips            <- 50
  tree            <- ape::rtree(n = tips, rooted = TRUE)
  tree$node.label <- truncnorm::rtruncnorm(n = ape::Nnode(tree), sd = 10, mean = 85, a = 0, b = 100) # dummy tree bootstrap values

  # Create continous phenotype
  phenotype_matrix <- as.matrix(phytools::fastBM(tree))

  # Create genotypes
  genotype_matrix <- matrix(NA, nrow = tips, ncol = 100)
  for (i in 1:ncol(genotype_matrix)) {
    genotype_matrix[ , i] <- stats::rbinom(tips, 1, 0.5)
  }
  row.names(genotype_matrix)  <- tree$tip.label
  colnames(genotype_matrix) <- paste("snp", c(1:100), sep = "_")

  # Check and return output ----------------------------------------------------
  check_for_root_and_bootstrap(tree)
  check_dimensions(phenotype_matrix, ape::Ntip(tree), 2, 1, 1)
  check_dimensions(genotype_matrix, ape::Ntip(tree), 2, NULL, 1)
  check_if_binary_matrix(genotype_matrix)
  check_rownames(phenotype_matrix, tree)
  check_rownames(genotype_matrix, tree)

  results <- list("tree" = tree, "phenotype" = phenotype_matrix, "genotype" = genotype_matrix)
  return(results)
} # end create_test_data()

prepare_phenotype <- function(pheno, disc_cont, tr){
  # Function description -------------------------------------------------------
  # Prepare phenotype for downstream analysese.
  # Convert phenotype matrix to phenotype vector.
  # Check if the continuous phenotype is normally distributed.
  # Check if convergence is possible in the phenotype given the tree.
  #
  # Inputs:
  # pheno. Matrix with 1 column.
  # disc_cont. String. Either "discrete" or "continuous"
  # tr. Phylo. Ntip = nrow(pheno).
  #
  # Outputs:
  # pheno_vector. Vector. Length = Ntip(tr).
  #
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_dimensions(pheno, exact_rows = ape::Ntip(tr), exact_cols = 1, min_rows = 1, min_cols = 1)
  check_str_is_discrete_or_continuous(disc_cont)

  # Function -------------------------------------------------------------------
  check_if_phenotype_normal(pheno, disc_cont)
  # check_if_convergence_occurs(pheno, tr, disc_cont)
  pheno_vector <- convert_matrix_to_vector(pheno)
  check_convergence_possible(pheno_vector, disc_cont)

  # Check and return output --------------------------------------------------------------
  check_equal(length(pheno_vector), ape::Ntip(tr))
  return(pheno_vector)
} # end prepare_phenotype

# END OF SCRIPT ---------------------------------------------------------------#
