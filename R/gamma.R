#' Calculate gamma within PhyC test
#'
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
#' @return results. List of three objects:
#'   \describe{
#'     \item{gamma_avg}{Numeric. Average gamma value of all genotypes. A single
#'     value.}
#'     \item{gamma_percent}{Numeric vector. Gamma value for each genotype.
#'     Length == number of genotypes.}
#'    \item{gamma_count}{Numeric vector. Raw gamma value for each genotype.
#'    Length == number of genotypes.}
#'   }
#' @noRd
calculate_phyc_gamma <- function(geno_trans_edge_list,
                                 pheno_recon_vec,
                                 high_conf){
  high_conf_edge_list <- high_conf$high_conf_ordered_by_edges

  check_equal(length(geno_trans_edge_list), length(high_conf_edge_list))
  check_equal(length(geno_trans_edge_list[[1]]), length(high_conf_edge_list[[1]]))
  check_equal(length(geno_trans_edge_list[[1]]), length(pheno_recon_vec))
  epsilon <- geno_beta <- gamma_count <- gamma_percent <-
    rep(0, length(geno_trans_edge_list))
  pheno_beta <- sum(pheno_recon_vec == 1 & high_conf$tr_and_pheno_hi_conf == 1)

  for (i in 1:length(geno_trans_edge_list)) {
    pheno_1_geno_0_to_1 <-
      sum(pheno_recon_vec == 1 &
            geno_trans_edge_list[[i]] == 1 &
            high_conf_edge_list[[i]] == 1)
    gamma_count[i] <- pheno_1_geno_0_to_1
    gamma_percent[i] <- gamma_count[i] / sum(high_conf_edge_list[[i]])
    geno_beta[i] <- sum(geno_trans_edge_list[[i]] == 1 &
                          high_conf_edge_list[[i]] == 1)
    epsilon[i] <- (2 * gamma_count[i]) / (pheno_beta + geno_beta[i])
  }
  gamma_avg <- mean(gamma_percent)
  num_hi_conf_edges <- unlist(lapply(high_conf_edge_list, sum))
  results <- list("gamma_avg" = gamma_avg,
                  "gamma_percent" = gamma_percent,
                  "gamma_count" = gamma_count,
                  "num_hi_conf_edges" = num_hi_conf_edges,
                  "pheno_beta" = pheno_beta,
                  "geno_beta" = geno_beta,
                  "epsilon" = epsilon)
  return(results)
}

#' Calculate gamma within synchronous test
#'
#' @description Given phenotype and genotype information, calculate a summary
#'   statistic that describes the number of edges on the tree where the
#'   phenotype and the genotype transitions (0 to 1 or 1 to 0). This summary
#'   stat will be used to evaluate the appropriateness / effectiveness of
#'   hogwash on this dataset.
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
#' @return results. List of three objects:
#'   \describe{
#'     \item{gamma_avg}{Numeric. Average gamma value of all genotypes. A single
#'     value.}
#'     \item{gamma_percent}{Numeric vector. Gamma value for each genotype.
#'     Length == number of genotypes.}
#'     \item{gamma_count}{Numeric vector. Raw gamma value for each genotype.
#'     Length == number of genotypes.}
#'   }
#' @noRd
calculate_synchronous_gamma <- function(geno_trans_edge_list,
                                        pheno_trans_vec,
                                        high_conf){
  high_conf_edge_list <- high_conf$high_conf_ordered_by_edges
  check_equal(length(geno_trans_edge_list), length(high_conf_edge_list))
  check_equal(length(geno_trans_edge_list[[1]]), length(high_conf_edge_list[[1]]))
  check_equal(length(geno_trans_edge_list[[1]]), length(pheno_trans_vec$transition))
  epsilon <- geno_beta <- gamma_count <- gamma_percent <-
    rep(0, length(geno_trans_edge_list))
  pheno_beta <-
    sum(pheno_trans_vec$transition == 1 & high_conf$tr_and_pheno_hi_conf == 1)
  for (i in 1:length(geno_trans_edge_list)) {
    pheno_trans_and_geno_trans <-
      sum(pheno_trans_vec$transition == 1 &
            geno_trans_edge_list[[i]] == 1 &
            high_conf_edge_list[[i]] == 1)
    gamma_count[i] <- pheno_trans_and_geno_trans
    gamma_percent[i] <- gamma_count[i] / sum(high_conf_edge_list[[i]])
    geno_beta[i] <- sum(geno_trans_edge_list[[i]] == 1 &
                          high_conf_edge_list[[i]] == 1)
    epsilon[i] <- (2 * gamma_count[i]) / (pheno_beta + geno_beta[i])
  }
  gamma_avg <- mean(gamma_percent)
  num_hi_conf_edges <- unlist(lapply(high_conf_edge_list, sum))
  results <- list("gamma_avg" = gamma_avg,
                  "gamma_percent" = gamma_percent,
                  "gamma_count" = gamma_count,
                  "num_hi_conf_edges" = num_hi_conf_edges,
                  "pheno_beta" = pheno_beta,
                  "geno_beta" = geno_beta,
                  "epsilon" = epsilon)
  return(results)
}

#' Calculate gamma within continuous test
#'
#' @description Given phenotype and genotype information, calculate a summary
#'   statistic that describes the number of edges on the tree where the
#'   phenotype change is large (delta phenotype > median(delta phenotype)) and
#'   the genotype transitions (0 to 1 or 1 to 0). This summary stat will be used
#'   to evaluate the appropriateness / effectiveness of hogwash on this dataset.
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
#' @return results. List of three objects:
#'   \describe{
#'     \item{gamma_avg}{Numeric. Average gamma value of all genotypes. A single
#'     value.}
#'     \item{gamma_percent}{Numeric vector. Gamma value for each genotype.
#'     Length == number of genotypes.}
#'     \item{gamma_count}{Numeric vector. Raw gamma value for each genotype.
#'     Length == number of genotypes.}
#'   }
#' @noRd
calculate_continuous_gamma <- function(pheno_recon_mat,
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
  pheno_delta <- rep(0, nrow(pheno_recon_mat))
  for (i in 1:nrow(pheno_recon_mat)) {
    pheno_delta[i] <- calculate_phenotype_change_on_edge(i, pheno_recon_mat)
  }
  pheno_beta <- sum(pheno_delta > stats::median(pheno_delta) &
                      high_conf$tr_and_pheno_hi_conf == 1)
  epsilon <- geno_beta <- gamma_count <- gamma_percent <-
    rep(0, length(geno_trans_edge_list$transition))
  for (i in 1:length(geno_trans_edge_list)) {
    pheno_large_geno_trans <-
      sum(pheno_delta > stats::median(pheno_delta) &
            geno_trans_edge_list[[i]]$transition == 1 &
            high_conf_edge_list[[i]] == 1)
    gamma_count[i] <- pheno_large_geno_trans
    gamma_percent[i] <- gamma_count[i] / sum(high_conf_edge_list[[i]])
    geno_beta[i] <- sum(geno_trans_edge_list[[i]]$transition == 1 &
                          high_conf_edge_list[[i]] == 1)
    epsilon[i] <- (2 * gamma_count[i]) / (geno_beta[i] + pheno_beta)
  }
  gamma_avg <- mean(gamma_percent)
  num_hi_conf_edges <- unlist(lapply(high_conf_edge_list, sum))
  results <- list("gamma_avg" = gamma_avg,
                  "gamma_percent" = gamma_percent,
                  "gamma_count" = gamma_count,
                  "num_hi_conf_edges" = num_hi_conf_edges,
                  "pheno_beta" = pheno_beta,
                  "geno_beta" = geno_beta,
                  "epsilon" = epsilon)
  return(results)
}
