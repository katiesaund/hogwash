hogwash <- function(pheno,
                    geno,
                    tree,
                    file_name,
                    dir,
                    perm = 10000,
                    fdr = 0.15,
                    bootstrap = 0.70,
                    annotation = NULL,
                    group_genotype_key = NULL){

  args                        <- NULL
  args$test                   <- FALSE
  args$phenotype              <- read_in_tsv_matrix(pheno)
  args$genotype               <- read_in_tsv_matrix(geno)
  args$tree                   <- read.tree(tree)
  args$tree$node.label[1]     <- 0
  args$tree$node.label        <- as.numeric(args$tree$node.label)
  args$output_name            <- file_name
  args$output_dir             <- dir
  args$perm                   <- as.numeric(perm) #typically 10,000
  args$fdr                    <- as.numeric(fdr)
  args$annot                  <- NULL
  if (!is.null(annotation)){
   args$annot                 <- read_in_tsv_matrix(annotation)
  }
  args$discrete_or_continuous <- check_input_format(args$phenotype, args$tree, args$genotype, args$output_name, args$output_dir, args$perm, args$fdr, args$annot)
  args$bootstrap_cutoff       <- as.numeric(bootstrap)
  args$group_genotype         <- NULL
  if (!is.null(args$group_genotype)){
    args$gene_snp_lookup      <- read_in_tsv_matrix(group_genotype_key)
  }
  select_test_type(args)
}

# End of script ----------------------------------------------------------------
