#' Determine if phenotype is discrete or continuous
#'
#' @description Determine if the matrix is discrete or continuous.
#'
#' @param mat Matrix. Should be a phenotype with 1 column.
#'
#' @return type. Character. Either "discrete" or "continuous".
#'
#' @noRd
assign_pheno_type <- function(mat){
  # Check input ----------------------------------------------------------------
  check_dimensions(mat, NULL, 1, 1, 1)
  if (!typeof(mat) %in% c("double", "integer", "numeric")) {
    stop("Must be a numeric matrix")
  }

  # Function -------------------------------------------------------------------
  type <- "discrete"
  if (sum(!(mat %in% c(0, 1))) > 0) {
    type <- "continuous"
  }

  # Check and return output ----------------------------------------------------
  check_is_string(type)
  check_str_is_discrete_or_continuous(type)
  return(type)
}

#' Calculate phenotype change per tree edge
#'
#' @description Quantify absoluate value of phenotype change on each tree edge.
#'
#' @param edge_list Numeric vector. Each number is the index of the tree edge to
#'   be used
#' @param phenotype_by_edges Mat. Dimensions: Nedge x 2 matrix. Entries are the
#'   phenotype value at the node, where row is the edge, 1st column is the
#'   parent node and 2nd column is the child node.
#'
#' @return Numeric vector.   Length = length(edge_list).
#' @noRd
calculate_phenotype_change_on_edge <- function(edge_list, phenotype_by_edges){
  # Check input ----------------------------------------------------------------
  if (max(edge_list) > nrow(phenotype_by_edges)) {
    stop("Cannot calculate phenotype change on edges.")
  }
  check_dimensions(phenotype_by_edges,
                   exact_rows = NULL,
                   min_rows = max(edge_list),
                   exact_cols = NULL,
                   min_cols = 2)
  if (!is.vector(edge_list)) {
    stop("edge_list must be a vector of indices of edges")
  }
  check_is_number(edge_list[1])

  # Function -------------------------------------------------------------------
  delta <- rep(NA, length(unique(edge_list)))
  for (j in 1:length(edge_list)) {
    delta[j] <-
      abs(phenotype_by_edges[edge_list[j], 1] -
            phenotype_by_edges[edge_list[j], 2])
  }

  # Check and return output ----------------------------------------------------
  if (sum(delta < 0) > 0) {
    stop("Delta phenotype should always be recorded as an absolute value.")
  }

  return(delta)
}

#' Calculate the phenotype difference on each tree edge
#'
#' @description Subtract child node phenotype value from parent node phenotype
#'   value. Identical to calculate_phenotype_change_on_edge EXCEPT that
#'   calculating the raw difference rather than the absolute value of the
#'   difference.
#'
#' @param edge_list Numeric vector. Each number is the index of the tree edge to
#'   be used.
#' @param phenotype_by_edges Mat. Dimensions: Nedge x 2 matrix. Entries are the
#'   phenotype value at the node, where row is the edge, 1st column is the
#'   parent node and 2nd column is the child node.
#'
#' @return delta. Numeric vector. Length = length(edge_list).
#' @noRd
calc_raw_diff <- function(edge_list, phenotype_by_edges){
  # Check input ----------------------------------------------------------------
  if (max(edge_list) > nrow(phenotype_by_edges)) {
    stop("Cannot calculate phenotype change on edges.")
  }
  check_dimensions(phenotype_by_edges,
                   exact_rows = NULL,
                   min_rows = max(edge_list),
                   exact_cols = NULL,
                   min_cols = 2)
  if (!is.vector(edge_list)) {
    stop("edge_list must be a vector of indices of edges")
  }
  check_is_number(edge_list[1])

  # Function -------------------------------------------------------------------
  delta <- rep(NA, length(unique(edge_list)))
  for (j in 1:length(edge_list)) {
    delta[j] <-
      phenotype_by_edges[edge_list[j], 1] - phenotype_by_edges[edge_list[j], 2]
  }

  # Return output --------------------------------------------------------------
  return(delta)
}

#' Report best fitting model for the phenotype and tree
#'
#' @param pheno Matrix. One column matrix. Row.names correspond to tree$tip.labels.
#' @param tree Phylo.
#' @param root Logical. Default = TRUE If TRUE, midpoint roots the tree.
#'
#' @export
report_phylogenetic_signal <- function(pheno, tree, root = TRUE){

  # Checks here are written out explicitly so that user sees specific, helpful
  # error messages.
  if (!methods::is(pheno, "matrix")) {
    stop("Provide a phenotype matrix.")
  }
  if (!methods::is(tree, "phylo")) {
    stop("Provide a tree in a phylo object")
  }
  if (nrow(pheno) != ape::Ntip(tree)) {
    stop("Provide a tree with tips that correspond to the phenotype rows.")
  }
  if (sum(row.names(pheno) != tree$tip.labels) > 0) {
    stop("Provide a phenotype matrix where row.names match tree$tip.labels.")
  }
  if (ncol(pheno) != 1) {
    stop("Provide a phenotype matrix with exactly 1 column of data.")
  }
#TODO unit tests
  cont_or_disc <- assign_pheno_type(pheno)

  if (root) {
    tree <- phytools::midpoint.root(tree)
    # TODO reorder phenotype by tree tips
  }

  if (cont_or_disc == "continuous") {
    lambda <- calculate_lambda(pheno, tree)
    report_lambda(lambda)
  } else {
    d_stat <- calculate_d(pheno, tree)
    report_d(d_stat)
  }
}

# TODO add documentation and unit tests
calculate_lambda <- function(pheno, tree) {
  # Checks here are written out explicitly so that user sees specific, helpful
  # error messages.
  if (!methods::is(pheno, "matrix")) {
    stop("Provide a phenotype matrix.")
  }
  if (!methods::is(tree, "phylo")) {
    stop("Provide a tree in a phylo object")
  }
  if (nrow(pheno) != ape::Ntip(tree)) {
    stop("Provide a tree with tips that correspond to the phenotype rows.")
  }
  if (sum(row.names(pheno) != tree$tip.labels) > 0) {
    stop("Provide a phenotype matrix where row.names match tree$tip.labels.")
  }
  if (ncol(pheno) != 1) {
    stop("Provide a phenotype matrix with exactly 1 column of data.")
  }

  # Function
  lambda_out <-
    phytools::phylosig(tree, pheno, method = "lambda", test = TRUE)
  lambda <- lambda_out$lambda

  return(lambda)
}

report_lambda <- function(lambda) {
  msg <- "The phenotype is more conserved (clumpy) than expected by Brownian Motion; Pagel's lambda = "
  if (lambda < 1.05 & lambda > 0.95) {
    msg <- "The phenotype is modeled well by Brownian Motion; Pagel's lambda = "
  } else if (lambda < 0.05 & lambda > -0.05) {
    msg <- "The phenotype is modeled well by White Noise (Random); Pagel's lambda = "
  } else if (lambda < 0.95 & lambda > 0.05) {
    msg <- "The phenotype is in-between a Brownian Motion and White Noise model; Pagel's lambda = "
  } else if (lambda < -0.05) {
    msg <- "The phenotype is negatively correlated among closely related species; Pagel's lambda = "
  }
  print(paste0(msg, round(lambda, 5)))
}

calculate_d <- function(pheno, tree) {
  #TODO add unit tests and documentation
  # Checks here are written out explicitly so that user sees specific, helpful
  # error messages.
  if (!methods::is(pheno, "matrix")) {
    stop("Provide a phenotype matrix.")
  }
  if (!methods::is(tree, "phylo")) {
    stop("Provide a tree in a phylo object")
  }
  if (nrow(pheno) != ape::Ntip(tree)) {
    stop("Provide a tree with tips that correspond to the phenotype rows.")
  }
  if (sum(row.names(pheno) != tree$tip.labels) > 0) {
    stop("Provide a phenotype matrix where row.names match tree$tip.labels.")
  }
  if (ncol(pheno) != 1) {
    stop("Provide a phenotype matrix with exactly 1 column of data.")
  }
  if (sum(sort(unique(pheno)) != c(0, 1)) > 0) {
    stop("Provide a binary phenotype (0/1)")
  }
  # Function
  tree$node.label <- NULL
  temp_trait <- pheno[, 1, drop = FALSE]
  temp_trait <- convert_trait_vec_to_df(temp_trait, tree)
  compar_data_obj <-
    caper::comparative.data(data = temp_trait,
                            phy = tree,
                            names.col =  ID)
  d_stat_results <- caper::phylo.d(data = compar_data_obj,
                                   names.col = ID,
                                   binvar = trait,
                                   phy = tree)
  d_stat <- d_stat_results$DEstimate
  return(d_stat)
}

report_d <- function(d_stat){
  msg <- "The phenotype is negatively correlated among closely related species; D = "
  if (d_stat < 1.05 & d_stat > 0.95) { # 1
    msg <- "The phenotype is modeled well by White Noise (Random); D = "
  } else if (d_stat < 0.05 & d_stat > -0.05) { #0
    msg <- "The phenotype is modeled well by Brownian Motion; D = "
  } else if (d_stat < 0.95 & d_stat > 0.05) {
    msg <- "The phenotype is in-between a Brownian Motion and White Noise model; D = "
  } else if (d_stat < -0.05) { # <0
    msg <- "The phenotype is more clumpy than expected by Brownian Motion; D = "
  }
  print(paste0(msg, d_stat))
}

convert_trait_vec_to_df <- function(trait_vec, tree){
  #TODO Add unit tests
  trait_df <- as.data.frame(trait_vec)
  trait_df <- cbind(tree$tip.label, trait_df)
  colnames(trait_df) <- c("ID", "trait")
  row.names(trait_df) <- NULL
  return(trait_df)
}

 #' Report best fitting model for the phenotype and tree
#'
#' @param pheno Matrix. One column matrix. Row.names correspond to tree$tip.labels.
#' @param tree Phylo.
#'
#' @export
internal_report_phylogenetic_signal <- function(pheno, tree){
# TODO incorporate this function within hogwash carefully!
  check_class(pheno, "matrix")
  check_class(tree, "phylo")
  check_equal(nrow(pheno), ape::Ntip(tree))
  if (sum(row.names(pheno) != tree$tip.labels) > 0) {
    stop("Provide a phenotype matrix where row.names match tree$tip.labels.")
  }
  check_dimensions(pheno, exact_cols = 1, min_rows = 1, min_cols = 1)
  #TODO unit tests
  cont_or_disc <- assign_pheno_type(pheno)

  # if (root) { ?? include this?
  #   tree <- phytools::midpoint.root(tree)
  #   # TODO reorder phenotype by tree tips
  # }

  if (cont_or_disc == "continuous") {
    lambda <- calculate_lambda(pheno, tree)
    report_lambda(lambda)
  } else {
    d_stat <- calculate_d(pheno, tree)
    report_d(d_stat)
  }
}
