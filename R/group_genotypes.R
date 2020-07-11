#' Build grouped genotype
#'
#' @description Build presence/absence of the grouped genotypes (e.g. gene) from
#'   the presence/absence of the ungrouped genotypes (e.g. snps).
#'
#' @param geno Matrix. Columns are genotypes. Rows are isolates.
#' @param gene_to_snp_lookup_table Matrix. 2 columns. 1st column are the
#'   ungrouped genotypes 2nd column are the grouped genotypes. The table's 1st
#'   column must contain only genotypes that are found in the row.names(geno).
#'
#' @return samples_by_genes. Matrix.
#'
#' @noRd
build_gene_genotype_from_snps <- function(geno, gene_to_snp_lookup_table){
  # Check input ----------------------------------------------------------------
  check_class(geno, "matrix")
  check_class(gene_to_snp_lookup_table, "matrix")
  if (sum(!gene_to_snp_lookup_table[ , 1] %in% colnames(geno)) > 0) {
    stop("gene_to_snp_lookup_table must only contain genotypes in geno")
  }
  check_dimensions(gene_to_snp_lookup_table, NULL, 1, 2, 2)
  check_dimensions(geno, NULL, 1, NULL, 1)

  # Function -------------------------------------------------------------------
  unique_genes <- unique(gene_to_snp_lookup_table[, 2])
  samples_by_genes <- matrix(0, nrow = nrow(geno), ncol = length(unique_genes))
  colnames(samples_by_genes) <- unique_genes
  row.names(samples_by_genes) <- row.names(geno)

  for (j in 1:length(unique_genes)) {
    temp_mat <-
      geno[, colnames(geno) %in% gene_to_snp_lookup_table[ , 1][gene_to_snp_lookup_table[, 2] == unique_genes[j]],
           drop = FALSE]
    class(temp_mat) <- "numeric"
    temp_column <- rowSums(temp_mat)
    samples_by_genes[, j] <- temp_column
  }

  samples_by_genes <- samples_by_genes > 0
  class(samples_by_genes) <- "numeric"
  # Return output --------------------------------------------------------------
  return(samples_by_genes)
}

#' Remove rare and common genotypes from grouped genotypes
#'
#' @description Remove genotypes that are too common (omnipresent) or rare
#'   (completely absent) for ancestral reconstruction to work from both the
#'   genotype and the lookup key.
#'
#' @param geno Matrix. Binary. Nrow = Ntip(tr). Ncol = number of original
#'   genotypes.
#' @param lookup Matrix. Ncol = 2. Nrow = genotypes with group assignments.
#'
#' @return List of four objects:
#'   \describe{
#'     \item{snp_per_gene.}{Named table. Names are genotypes. Values are number
#'     of not-yet-grouped-genotypes that go into the grouped genotype.}
#'     \item{unique_genes.}{Character vector. Vector of genotype names.}
#'     \item{gene_snp_lookup.}{Character matrix. Ncol = 2. Nrow = number of
#'     genotypes that are neither omni-present or completely absent.}
#'     \item{genotype.}{Matrix.}
#'   }
#' @noRd
prepare_grouped_genotype <- function(geno, lookup){
  # Check input ----------------------------------------------------------------
  check_dimensions(lookup, min_rows = 1, exact_cols = 2, min_cols = 2)
  check_dimensions(geno, min_cols = 1, min_rows = 1)
  check_if_binary_matrix(geno)

  # Function -------------------------------------------------------------------
  geno_to_drop_bc_not_present <-
    c(colnames(geno)[colSums(geno) == 0],
      colnames(geno)[colSums(geno) == nrow(geno)])

  # we don't want to remove snps that are too rare or too common until snps are
  # grouped into genes, then run on the grouped genes. But we need to remove
  # SNPs that don't occur for ape::ace to work.
  genotype <- geno[, colSums(geno) > 0]
  genotype <- genotype[, colSums(genotype) < nrow(genotype)]

  gene_snp_lookup <-
    lookup[!(lookup[, 1] %in% geno_to_drop_bc_not_present),
           , drop = FALSE]
  gene_snp_lookup <-
    gene_snp_lookup[gene_snp_lookup[, 1] %in% colnames(genotype),
                    , drop = FALSE]
  unique_genes <- unique(gene_snp_lookup[, 2])
  snps_per_gene <- table(gene_snp_lookup[, 2])

  genotype <-
    genotype[, colnames(genotype) %in% gene_snp_lookup[, 1], drop = FALSE]
  # Check intermediate data ----------------------------------------------------
  if (ncol(genotype) < 2) {
    stop("There are fewer than 2 genotypes that have variable presence/absence and are named in the grouping key")
  }
  if (nrow(gene_snp_lookup) < 2) {
    stop("There are fewer than 2 genotypes that are named in the grouping key and found in the genotype matrix")
  }

  # Do the grouping step ----
  grouped_genotype <- build_gene_genotype_from_snps(genotype, gene_snp_lookup)

  # Return output ----
  results <- list("snps_per_gene" = snps_per_gene,
                  "unique_genes" = unique_genes,
                  "gene_snp_lookup" = gene_snp_lookup,
                  "genotype" = grouped_genotype)
  return(results)
}

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
prepare_genotype <- function(group_logical, geno, tr, lookup){
  # Check input ----------------------------------------------------------------
  check_class(group_logical, "logical")
  check_for_root_and_bootstrap(tr)
  if (!is.null(lookup)) {
    check_dimensions(lookup, exact_cols = 2, min_cols = 2, min_rows = 1)
  }
  check_dimensions(geno, exact_rows = ape::Ntip(tr), min_rows = 1, min_cols = 1)
  #
  # Function -------------------------------------------------------------------
  if (group_logical) {
    prepped_geno <- prepare_grouped_genotype(geno, lookup)
  } else {
    prepped_geno <- prepare_ungrouped_genotype(geno, tr)
  }
  return(prepped_geno)
}
