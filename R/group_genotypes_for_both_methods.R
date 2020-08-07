# For functions that are used by both grouping methods -- either pre- or post-AR
# In practice, these are the functions that were in the code back at the time
# of the first arxiv submission and were unchanged when I rashly updated the
# grouping algorithm in July.

#' Remove common and rare ungrouped genotypes
#'
#' @description Remove genotypes that are too common or rare for ancestral
#'   reconstruction to work. Given that this genotype is not grouped return NULL
#'   for the variable snps_per_gene. Keep track of which genotypes got removed
#'   in no_convergence_genotypes.
#'
#' @param geno Matrix. Binary. Ncol = number of genotypes. Nrow = Ntip(tr).
#' @param tr Phylo.
#'
#' @return List of three objects:
#'   \describe{
#'     \item{snp_per_gene}{NULL.}
#'     \item{genotype}{Matrix.}
#'     \item{convergene_not_possible_genotypes}{Character vector. Vector of
#'     genotype names.}
#'   }
#' @noRd
prepare_ungrouped_genotype <- function(geno, tr){
  # Check input ----------------------------------------------------------------
  check_dimensions(geno, exact_rows = ape::Ntip(tr), min_rows = 1, min_cols = 1)
  check_if_binary_matrix(geno)
  check_for_root_and_bootstrap(tr)

  # Function -------------------------------------------------------------------
  simple_geno <- remove_rare_or_common_geno(geno, tr)
  snps_per_gene <- NULL

  # Check and return output --------------------------------------------------
  results <-
    list("snps_per_gene" = snps_per_gene,
         "genotype" = simple_geno$mat,
         "no_convergence_genotypes" = simple_geno$dropped_genotype_names)
  return(results)
}


#' Prepare genotype, paying attention to (not) grouping
#'
#' @description Funnel the genotype to be prepared for downstream use. The
#'   preparation depends on if the genotype is going to be grouped or not.
#' @param group_logical Logical.
#' @param geno Genotype matrix. Binary. Nrow = Ntip(tr). Ncol = number of
#'   ungrouped genotypes.
#' @param tr Phylo.
#' @param lookup Either NULL or a Matrix. Ncol = 2.
#'
#' @return  prepped_geno. List of multiple objects. Content depends on value of
#'   group_logical. For grouped genotypes output includes:
#'   \describe{
#'     \item{snp_per_gene}{Named table. Names are genotypes. Values are number
#'     of not-yet-grouped-genotypes that go into the grouped genotype.}
#'     \item{unique_genes}{Character vector. Vector of genotype names.}
#'     \item{gene_snp_lookup}{Character matrix. Ncol = 2. Nrow = number of
#'     genotypes that are neither omni-present or completely absent.}
#'     \item{genotype}{Matrix.}
#'
#'   }
#'   For not grouped genotypes output includes:
#'   \describe{
#'     \item{snp_per_gene}{NULL}
#'     \item{genotype}{Matrix}
#'     \item{convergene_not_possible_genotypes}{Character vector. Vector of
#'     genotype names.}
#'
#'   }
#' @noRd
prepare_genotype <- function(group_logical, geno, tr, lookup, group_method){
  # Check input ----------------------------------------------------------------
  check_class(group_logical, "logical")
  check_for_root_and_bootstrap(tr)
  if (!is.null(lookup)) {
    check_dimensions(lookup, exact_cols = 2, min_cols = 2, min_rows = 1)
  }
  check_dimensions(geno, exact_rows = ape::Ntip(tr), min_rows = 1, min_cols = 1)
  check_is_string(group_method)
  if (!group_method %in% c("post-ar", "pre-ar")) {
    stop("If grouping the genotype please select either: post-ar or pre-ar")
  }

  # Function -------------------------------------------------------------------
  if (group_logical) {
    if (group_method == "post-ar") {
      prepped_geno <- prepare_grouped_genotype_post_ar(geno, lookup)
    } else if (group_method == "pre-ar") {
      prepped_geno <- prepare_grouped_genotype_pre_ar(geno, lookup, tr)
    } else {
      stop("Group method invalid")
    }
  } else {
    prepped_geno <- prepare_ungrouped_genotype(geno, tr)
  }
  return(prepped_geno)
}
