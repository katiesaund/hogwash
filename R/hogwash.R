#' Bacterial genome-wide association study.
#'
#' @description This function runs a bacterial genome-wide association test. It
#'   runs either the Continuous Test (when given continuous phenotype data) or
#'   both the Synchronous Test and PhyC (when given binary phenotype data).
#'

#' @param pheno Matrix. Dimensions: nrow = number of samples, ncol = 1. Either
#'   continuous or binary (0/1). Row.names() must match tree$tip.label.
#' @param geno Matrix. Dimensions: nrow = number of samples, ncol = number of
#'   genotypes. Binary (0/1). Row.names() must match tree$tip.label.
#' @param tree Phylo object. If unrooted, will be rooted using
#'   phytools::midpoint.root() method.
#' @param file_name Character. Name of output files.
#' @param dir Character. Path to output directory.
#' @param perm Integer. Number of permutations to run.
#' @param fdr Numeric. False discovery rate. Between 0 and 1.
#' @param bootstrap Numeric. Confidence threshold for tree bootstrap values.
#' @param group_genotype_key Matrix. Dimenions: nrow = number of unique genotypes, ncol = 2.
#'
#' @export
#'
#' @author Katie Saund
#' @references Farhat MR, Shapiro BJ, Kieser KJ, et al. Genomic analysis identifies targets of convergent positive selection in drug-resistant Mycobacterium tuberculosis. Nat Genet. 2013;45(10):1183–1189. doi:10.1038/ng.2747

hogwash <- function(pheno,
                    geno,
                    tree,
                    file_name = paste(Sys.Date(), "_hogwash", sep = ""),
                    dir = ".",
                    perm = 10000,
                    fdr = 0.15,
                    bootstrap = 0.70,
                    group_genotype_key = NULL){
  args <- NULL
  args$phenotype <- pheno
  args$genotype <- geno
  args$tree <- tree
  args$output_name <- file_name
  args$output_dir <- dir
  args$perm <- perm
  args$fdr <- fdr
  args$bootstrap_cutoff <- bootstrap
  args$gene_snp_lookup <- group_genotype_key
  args$group_genotype <- FALSE
  if (!is.null(group_genotype_key)) {
    args$group_genotype <- TRUE
  }
  args$discrete_or_continuous <- check_input_format(args$phenotype, args$tree, args$genotype, args$output_name, args$output_dir, args$perm, args$fdr)
  select_test_type(args)
}

# End of script ----------------------------------------------------------------
