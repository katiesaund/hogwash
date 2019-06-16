#' Main hogwash function
#'
#' @export

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
  args$annot <- NULL
  args$gene_snp_lookup <- group_genotype_key
  group_genotype <- FALSE
  if (!is.null(group_genotype_key)) {
    group_genotype <- TRUE
  }
  args$discrete_or_continuous <- check_input_format(args$phenotype, args$tree, args$genotype, args$output_name, args$output_dir, args$perm, args$fdr, args$annot)
  select_test_type(args)
}

# End of script ----------------------------------------------------------------
