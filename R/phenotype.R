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

#' Check if the phenotype is normally distributed
#'
#' @description Given a continuous phenotype, check if the phenotype follows a
#'   normal distribution. The Brownian Motional model of ancestral
#'   reconstruction, which this library uses, assumes a normal distribution of
#'   the phenotype. If it doesn't print a statement, but continue with
#'   hogwash().
#'
#' @param pheno Matrix. A one column numeric matrix.
#' @param continuous_or_discrete String. Either 'continuous' or 'discrete.'
#'
#' @noRd
check_if_phenotype_normal <- function(pheno, continuous_or_discrete){
  # Check input ----------------------------------------------------------------
  check_dimensions(pheno, NULL, 1, 1, 1)
  if (!typeof(pheno) %in% c("double", "integer", "numeric")) {
    stop("Must be a numeric matrix")
  }
  check_is_string(continuous_or_discrete)
  check_str_is_discrete_or_continuous(continuous_or_discrete)

  # Function -------------------------------------------------------------------
  msg <- "Consider making phenotype normal so ancestral reconstruction
  assumptions are not violated."
  if (continuous_or_discrete == "continuous") {
    if (length(unique(pheno)) == 1) {
      stop("phenotype must have some variability")
    }
    result <- stats::shapiro.test(unlist(pheno))
    alpha <- 0.05
    if (result$p < alpha) {
      print(msg)
    }
  }
}
