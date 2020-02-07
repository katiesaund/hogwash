#' Check format of user provided data
#'
#' @description Check that all of the inputs into hogwash are in the valid
#'   format.
#' @param pheno Matrix. Phenotype matrix row.names should be the same as the
#'   tr$tip.label and exactly 1 column. Phenotype rownames should be in the same
#'   order as the tr$tip.label.
#' @param tr Phylo. Tree.
#' @param geno Matrix. Genotype matrix should have same rows as tr$tip.label, at
#'   least 2 genotypes in the columns. Genotype rownames should be in the same
#'   order as the tr$tip.label.
#' @param name Character. Output name (prefix).
#' @param dir Character. Output path.
#' @param perm Number. Times to shuffle the data on the tree to create a null
#'   distribution for the permutation test.
#' @param fdr Number. False discovery rate. Between 0 and 1.
#' @param bootstrap Number. Confidence threshold for tree bootstrap values.
#' @param group_genotype_key Either NULL or a matrix. Dimenions: nrow = number
#'   of unique genotypes, ncol = 2.
#'
#' @return discrete_or_continuous. Character. Either "discrete" or "continuous".
#'   Describes the input phenotype.
#' @noRd
check_input_format <- function(pheno,
                               tr,
                               geno,
                               name,
                               dir,
                               perm,
                               fdr,
                               bootstrap,
                               group_genotype_key){
  # Check input ----------------------------------------------------------------
  check_dimensions(geno, ape::Ntip(tr), 2, NULL, 2)
  check_dimensions(pheno, ape::Ntip(tr), 2, 1, 1)
  check_rownames(geno, tr)
  check_rownames(pheno, tr)
  check_for_NA_and_inf(geno)
  check_for_NA_and_inf(pheno)
  check_if_binary_matrix(geno)
  check_is_string(name)
  check_if_dir_exists(dir)
  check_if_permutation_num_valid(perm)
  check_num_between_0_and_1(fdr)
  check_num_between_0_and_1(bootstrap)
  if (!is.null(group_genotype_key)) {
    check_class(group_genotype_key, "matrix")
    check_dimensions(group_genotype_key,
                     min_rows = 1,
                     exact_cols = 2,
                     min_cols = 2)
    check_for_NA_and_inf(group_genotype_key)
    if (nrow(group_genotype_key) != ncol(geno)) {
      stop("If providing a genotype grouping key there must be a 1-to-1 mapping between each variant in the genotype key's first column and the variants in genotype matrix columns.")
    }
  }
  if (ape::Ntip(tr) < 7) {
    stop("Tree must have at least 7 tips")
  }


  # Function -------------------------------------------------------------------
  discrete_or_continuous <- assign_pheno_type(pheno)

  # Check and return output ----------------------------------------------------
  check_is_string(discrete_or_continuous)
  return(discrete_or_continuous)
}
