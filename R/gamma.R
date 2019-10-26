#' Calculate gamma within phyC test
#'
#' @param geno_trans_edge_list A list of genotypes. Length == number of
#'   genotypes. Length of individual vectors within == Nedge(tree). Individual
#'   vectors are binary.
#' @param pheno_recon_vec A binary vector with length == Nedge(tree).
#' @param high_conf_edge_list A list of high confidence edges. Length of list
#'  == number of genotypes. Length of individual vectors within == Nedge(tree).
#'  Individual vectors are binary.
#' @return results. List of three objects:
#'   $gamma_avg Numeric. Average gamma value of all genotypes. A single value.
#'   $gamma_percent Numeric vector. Gamma value for each genotype. Length ==
#'     number of genotypes.
#'  $gamma_count Numeric vector. Raw gamma value for each genotype. Length ==
#'    number of genotypes.
#' @noRd
#'
calculate_phyc_gamma <- function(geno_trans_edge_list,
                                 pheno_recon_vec,
                                 high_conf_edge_list){
  gamma_list <- gamma_percent <- rep(0, length(geno_trans_edge_list))
  for (i in 1:length(geno_trans_edge_list)) {
    pheno_1_geno_0_to_1 <-
      sum(pheno_recon_vec == 1 &
            geno_trans_edge_list[[i]] == 1 &
            high_conf_edge_list[[i]] == 1)
    gamma_list[[i]] <- pheno_1_geno_0_to_1
    gamma_percent[[i]] <- gamma_list[[i]] / sum(high_conf_edge_list[[i]])
  }
  gamma_avg <- mean(gamma_percent)

  results <- list("gamma_avg" = gamma_avg,
                  "gamma_percent" = gamma_percent,
                  "gamma_count" = gamma_list)
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
#'   Individual vectors are binary with length == Nedge(tree). Named vectors
#'   are $transition and $trans_dir.
#' @param high_conf_edge_list A list of high confidence edges. Length of list
#'  == number of genotypes. Length of individual vectors within == Nedge(tree).
#'  Individual vectors are binary.
#'
#' @return results. List of three objects:
#'   $gamma_avg Numeric. Average gamma value of all genotypes. A single value.
#'   $gamma_percent Numeric vector. Gamma value for each genotype. Length ==
#'     number of genotypes.
#'  $gamma_count Numeric vector. Raw gamma value for each genotype. Length ==
#'    number of genotypes.
#' @noRd
#'
calculate_synchronous_gamma <- function(geno_trans_edge_list,
                                        pheno_trans_vec,
                                        high_conf_edge_list){
  gamma_list <- gamma_percent <- rep(0, length(geno_trans_edge_list))
  for (i in 1:length(geno_trans_edge_list)) {
    pheno_1_geno_0_to_1 <-
      sum(pheno_trans_vec$transition == 1 &
            geno_trans_edge_list[[i]] == 1 &
            high_conf_edge_list[[i]] == 1)
    gamma_list[[i]] <- pheno_1_geno_0_to_1
    gamma_percent[[i]] <- gamma_list[[i]] / sum(high_conf_edge_list[[i]])
  }
  gamma_avg <- mean(gamma_percent)

  results <- list("gamma_avg" = gamma_avg,
                  "gamma_percent" = gamma_percent,
                  "gamma_count" = gamma_list)
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
#' @param high_conf_edge_list List of vectors. Length(list) = number of
#'   genotypes. Length(vectors) == nedge(tree). Binary.
#'
#' @return results. List of three objects:
#'   $gamma_avg Numeric. Average gamma value of all genotypes. A single value.
#'   $gamma_percent Numeric vector. Gamma value for each genotype. Length ==
#'     number of genotypes.
#'  $gamma_count Numeric vector. Raw gamma value for each genotype. Length ==
#'    number of genotypes.
#' @noRd
#'
calculate_continuous_gamma <- function(geno_trans_edge_list,
                                       pheno_recon_mat,
                                       high_conf_edge_list){

  pheno_delta <- rep(0, nrow(pheno_recon_mat))
  for (i in 1:nrow(pheno_recon_mat)) {
    pheno_delta[i] <- calculate_phenotype_change_on_edge(i, pheno_recon_mat)
  }

  gamma_list <- gamma_percent <- rep(0, length(geno_trans_edge_list$transition))
  for (i in 1:length(geno_trans_edge_list)) {
    pheno_large_geno_trans <-
      sum(pheno_delta > median(pheno_delta) &
            geno_trans_edge_list[[i]]$transition == 1 &
            high_conf_edge_list[[i]] == 1)
    gamma_list[[i]] <- pheno_large_geno_trans
    gamma_percent[[i]] <- gamma_list[[i]] / sum(high_conf_edge_list[[i]])
  }
  gamma_avg <- mean(gamma_percent)

  results <- list("gamma_avg" = gamma_avg,
                  "gamma_percent" = gamma_percent,
                  "gamma_count" = gamma_list)
  return(results)
}
