#' Calculate convergence (epsilon, beta genotype, and beta phenotype) within
#'   PhyC test
#' @description Given phenotype and genotype information, calculate summary
#'   statistics that describe the number of edges on the tree where the
#'   phenotype is present and the genotype transitions as well as their overlap.
#' @param geno_trans_edge_list A list of genotypes. Length == number of
#'   genotypes. Length of individual vectors within == Nedge(tree). Individual
#'   vectors are binary.
#' @param pheno_recon_vec A binary vector with length == Nedge(tree).
#' @param high_conf
#'   \describe{
#'    \item{tr_and_pheno_hi_conf}{A vector of high confidence edges. Only takes
#'    into acount the tree edge lengths, bootstrap values, and phenotype
#'    ancestral reconstruction ML values. Binary vector.}
#'     \item{high_conf_ordered_by_edges}{A list of high confidence edges. Length
#'     of list == number of genotypes. Length of individual vectors within ==
#'     Nedge(tree). Individual vectors are binary.}
#'   }
#' @return results. List of five objects:
#'   \describe{
#'     \item{intersection}{Numeric vector. Intersection of the geno_beta and
#'      pheno_beta  for each genotype. Length == number of genotypes.}
#'     \item{num_hi_conf_edges}{Numeric vector. Number of high confidence
#'     edges per genotype. Length == number of genotypes.}
#'     \item{pheno_beta}{Number. Count of how many tree edges are phenotype
#'     transitions and the phenotype ancestral reconstruction and tree edge are
#'     high confidence. Length == 1.}
#'     \item{geno_beta}{Numeric vector. count of how many tree edges are
#'     gentoype transitions and the genotype ancestral reconstruction and tree
#'     edge are high confidence. Length == number of genotypes.}
#'     \item{epsilon}{Numeric vector. 2 x (edges with both high confidence
#'      genotype transition and phenotype presence) / sum(edges with high
#'      confidence gentoype transition and/or phenotype presence). Length ==
#'      number of genotypes.}
#'   }
#' @noRd
calculate_phyc_gamma <- function(geno_trans_edge_list,
                                 pheno_recon_vec,
                                 high_conf){
  high_conf_edge_list <- high_conf$high_conf_ordered_by_edges

  check_equal(length(geno_trans_edge_list), length(high_conf_edge_list))
  check_equal(length(geno_trans_edge_list[[1]]), length(high_conf_edge_list[[1]]))
  check_equal(length(geno_trans_edge_list[[1]]), length(pheno_recon_vec))
  epsilon <- geno_beta <- intersection <- rep(0, length(geno_trans_edge_list))
  pheno_beta <- sum(pheno_recon_vec == 1 & high_conf$tr_and_pheno_hi_conf == 1)

  for (i in 1:length(geno_trans_edge_list)) {
    pheno_1_geno_0_to_1 <-
      sum(pheno_recon_vec == 1 &
            geno_trans_edge_list[[i]] == 1 &
            high_conf_edge_list[[i]] == 1)
    intersection[i] <- pheno_1_geno_0_to_1
    geno_beta[i] <- sum(geno_trans_edge_list[[i]] == 1 &
                          high_conf_edge_list[[i]] == 1)
    epsilon[i] <- (2 * intersection[i]) / (pheno_beta + geno_beta[i])
  }
  num_hi_conf_edges <- unlist(lapply(high_conf_edge_list, sum))
  results <- list("intersection" = intersection,
                  "num_hi_conf_edges" = num_hi_conf_edges,
                  "pheno_beta" = pheno_beta,
                  "geno_beta" = geno_beta,
                  "epsilon" = epsilon)
  return(results)
}

#' Calculate convergence (epsilon, beta genotype, and beta phenotype) within the
#'   Synchronous test
#'
#' @description Given phenotype and genotype information, calculate summary
#'   statistics that describe the number of edges on the tree where the
#'   phenotype and the genotype transitions (0 to 1 or 1 to 0) as well as their
#'   overlap.
#' @param geno_trans_edge_list A list of genotypes. Length == number of
#'   genotypes. Length of individual vectors within == Nedge(tree). Individual
#'   vectors are binary.
#' @param pheno_recon_vec A list of lists. Each sublist has two names vectors.
#'   Individual vectors are binary with length == Nedge(tree). Named vectors are
#'   $transition and $trans_dir.
#' @param high_conf
#'   \describe{
#'    \item{tr_and_pheno_hi_conf}{A vector of high confidence edges. Only takes
#'    into acount the tree edge lengths, bootstrap values, and phenotype
#'    ancestral reconstruction ML values. Binary vector.}
#'     \item{high_conf_ordered_by_edges}{A list of high confidence edges. Length
#'     of list == number of genotypes. Length of individual vectors within ==
#'     Nedge(tree). Individual vectors are binary.}
#'   }
#'
#' @return results. List of five objects:
#'   \describe{
#'     \item{intersection}{Numeric vector. Intersection of the geno_beta and
#'      pheno_beta for each genotype. Length == number of genotypes.}
#'     \item{num_hi_conf_edges}{Numeric vector. Number of high confidence
#'     edges per genotype. Length == number of genotypes.}
#'     \item{pheno_beta}{Number. Count of how many tree edges are phenotype
#'     transitions and the phenotype ancestral reconstruction and tree edge are
#'     high confidence. Length == 1.}
#'     \item{geno_beta}{Numeric vector. count of how many tree edges are
#'     gentoype transitions and the genotype ancestral reconstruction and tree
#'     edge are high confidence. Length == number of genotypes.}
#'     \item{epsilon}{Numeric vector. 2 x (edges with both high confidence
#'     genotype AND phenotype transition) / sum(edges with high confidence
#'     gentoype and/or phenotype transitions). Length == number of genotypes.}
#'   }
#' @noRd
calculate_synchronous_gamma <- function(geno_trans_edge_list,
                                        pheno_trans_vec,
                                        high_conf){
  high_conf_edge_list <- high_conf$high_conf_ordered_by_edges
  check_equal(length(geno_trans_edge_list), length(high_conf_edge_list))
  check_equal(length(geno_trans_edge_list[[1]]), length(high_conf_edge_list[[1]]))
  check_equal(length(geno_trans_edge_list[[1]]), length(pheno_trans_vec$transition))
  epsilon <- geno_beta <- intersection <- rep(0, length(geno_trans_edge_list))
  pheno_beta <-
    sum(pheno_trans_vec$transition == 1 & high_conf$tr_and_pheno_hi_conf == 1)
  for (i in 1:length(geno_trans_edge_list)) {
    pheno_trans_and_geno_trans <-
      sum(pheno_trans_vec$transition == 1 &
            geno_trans_edge_list[[i]] == 1 &
            high_conf_edge_list[[i]] == 1)
    intersection[i] <- pheno_trans_and_geno_trans
    geno_beta[i] <- sum(geno_trans_edge_list[[i]] == 1 &
                          high_conf_edge_list[[i]] == 1)
    epsilon[i] <- (2 * intersection[i]) / (pheno_beta + geno_beta[i])
  }
  num_hi_conf_edges <- unlist(lapply(high_conf_edge_list, sum))
  results <- list("intersection" = intersection,
                  "num_hi_conf_edges" = num_hi_conf_edges,
                  "pheno_beta" = pheno_beta,
                  "geno_beta" = geno_beta,
                  "epsilon" = epsilon)
  return(results)
}

#' Calculate convergence (epsilon, beta genotype, and beta phenotype) within the
#'   Continuous test
#'
#' @description Given phenotype and genotype information, calculate summary
#'   statistics that describes the elementwise intersection of genotype
#'   transition edges and phenotype delta divided by their union. Delta
#'   phenotype is scaled from 0 to 1.
#' @param geno_trans_edge_list List of lists. Each entry has two named lists:
#'   $transition and $transition_dir. Number of sub-lists = number of genotypes.
#'   Length(each sublist) == Nedge(tree). Binary.
#' @param pheno_recon_mat Matrix. Structured like tree$edge. Numeric.
#' @param high_conf
#'   \describe{
#'    \item{tr_and_pheno_hi_conf}{A vector of high confidence edges. Only takes
#'    into acount the tree edge lengths, bootstrap values, and phenotype
#'    ancestral reconstruction ML values. Binary vector.}
#'    \item{high_conf_ordered_by_edges}{A list of high confidence edges. Length
#'    of list == number of genotypes. Length of individual vectors within ==
#'    Nedge(tree). Individual vectors are binary.}
#'   }
#' @return results. List of five objects:
#'   \describe{
#'     \item{intersection}{Numeric vector. Intersection of the geno_beta and
#'      pheno_beta for each genotype. Length == number of genotypes.}
#'     \item{num_hi_conf_edges}{Number of high confidence edges.}
#'     \item{pheno_beta}{Numeric vector. Beta(phenotype). Length == number of
#'     genotypes.}
#'     \item{geno_beta}{Numeric vector. Beta(genotype). Length == number of
#'     genotypes.}
#'     \item{epsilon}{Numeric vector. Epsilon value for each genotype. Length ==
#'     number of genotypes.}
#'   }
#' @noRd
calculate_continuous_convergence <- function(pheno_recon_mat,
                                       high_conf){
  high_conf_edge_list <- high_conf$high_conf_ordered_by_edges
  geno_trans_edge_list <- high_conf$genotype_transition
  check_equal(length(geno_trans_edge_list), length(high_conf_edge_list))
  check_equal(length(geno_trans_edge_list[[1]]$transition),
              length(high_conf_edge_list[[1]]))
  check_equal(nrow(pheno_recon_mat), length(high_conf_edge_list[[1]]))
  check_dimensions(pheno_recon_mat,
                   exact_rows = length(high_conf_edge_list[[1]]),
                   min_rows = 1,
                   exact_cols = 2,
                   min_cols = 2)

  # Get indices for all high confidence genotype transition edges
  num_geno <- length(geno_trans_edge_list)
  geno_non_trans_index_list <- high_conf_index_list <- rep(list(0), num_geno)
  for (i in 1:num_geno) {
    geno_non_trans_index_list[[i]] <-
      which(geno_trans_edge_list[[i]]$transition == 0)
    high_conf_index_list[[i]] <-
      which(high_conf_edge_list[[i]] == 1 & high_conf$tr_and_pheno_hi_conf == 1)
  }

  # Get phenotype delta for all edges and specifically for genotype transition
  #   edges

  pheno_delta <- pheno_delta_geno_non_trans <- rep(list(0), num_geno)
  num_tree_edge <- nrow(pheno_recon_mat)

  for (i in 1:num_geno) {
    pheno_delta[[i]] <-
      calculate_phenotype_change_on_edge(1:num_tree_edge,
                                         pheno_recon_mat)
    pheno_delta_geno_non_trans[[i]] <-
      calculate_phenotype_change_on_edge(geno_non_trans_index_list[[i]],
                                         pheno_recon_mat)
  }

  pheno_beta <- epsilon <- geno_beta <- intersection <- rep(0, num_geno)

  for (i in 1:num_geno) {
    scaled_pheno <- scales::rescale(pheno_delta[[i]], to = c(0, 1))
    pheno_beta[i] <- sum(scaled_pheno * (1 * (high_conf_edge_list[[i]] == 1)))
    geno_beta[i] <- sum(geno_trans_edge_list[[i]]$transition == 1 &
                          high_conf_edge_list[[i]] == 1)
    intersection[i] <-
      sum(
        (scaled_pheno * (1 * (high_conf_edge_list[[i]] == 1))) *
          (geno_trans_edge_list[[i]]$transition == 1 &
             high_conf_edge_list[[i]] == 1))
    epsilon[i] <-
      intersection[i] / (geno_beta[i] + pheno_beta[i] - intersection[i])
  }
  num_hi_conf_edges <- unlist(lapply(high_conf_edge_list, sum))
  results <- list("intersection" = intersection,
                  "num_hi_conf_edges" = num_hi_conf_edges,
                  "pheno_beta" = pheno_beta,
                  "geno_beta" = geno_beta,
                  "epsilon" = epsilon)
  return(results)
}
