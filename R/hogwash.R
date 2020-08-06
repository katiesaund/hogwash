#' Run bacterial GWAS
#'
#' @description This function runs a bacterial genome-wide association test. It
#'   runs either the Continuous Test when given continuous phenotype data. When
#'   given binary data the user may run either the Synchronous Test or PhyC or
#'   both.
#'
#' @details Overview: hogwash reads in one phenotype (either continuous or
#'   binary), a matrix of binary genotypes, and a phylogenetic tree. Given these
#'   inputs it performs an ancestral reconstruction of that phenotype and each
#'   genotype. The ancestral reconstructions are used to perform one of several
#'   tests to associate the the genotypes with the phenotype:
#'   \enumerate{
#'    \item Continuous Test
#'    \item Synchronous Test
#'    \item PhyC Test (Farhat et al.)
#'   }
#'   Once a test finishes running it returns (i) p-values for all genotypes
#'   tested, (ii) a manhattan plot of those p-values; if any of the genotypes
#'   tested were significant associated with the phenotype after FDR correction
#'   it also returns (iii) a list of significant hits and (iv) figures of the
#'   genotype & phenotype reconstructions on the tree.
#'
#'   Grouping: A feature of hogwash is the ability to organize genotypes into
#'   biologically meaningful groups. Testing for an association between an
#'   individual SNP and a phenotype is quite stringent, but patterns may emerge
#'   when grouping together biologically related genotypes. For example,
#'   grouping together all variants (insertions, deletions and SNPs) within a
#'   gene or promoter region could allow the user to identify a particular gene
#'   as being associated with a phenotype while any individual variant within
#'   that gene may not have deep penetrance in the isolates being tested.
#'   Grouped genotypes could have increased power to identify convergent
#'   evolution because they captures larger trends in functional impact at the
#'   group level and reduce the multiple testing correction burden. Use cases
#'   for this method could be to group SNPs into genes, kmers or genes into
#'   pathways, etc... Each of the three tests can be run on disaggregated data
#'   or aggregated data with the inclusion of a grouping key.
#'
#' @param pheno Matrix. Dimensions: nrow = number of samples, ncol = 1. Either
#'   continuous or binary (0/1). Row.names() must match tree$tip.label. Required
#'   input.
#' @param geno Matrix. Dimensions: nrow = number of samples, ncol = number of
#'   genotypes. Binary (0/1). Row.names() must match tree$tip.label. Required
#'   input.
#' @param tree Phylo object. If unrooted, will be rooted using
#'   phytools::midpoint.root() method. Required input.
#' @param file_name Character. Suffix for output files. Default value is the
#'   current date: YYYY-MM-DD.
#' @param dir Character. Path to output directory. Default value is current
#'   directory: "."
#' @param perm Integer. Number of permutations to run. Default value is: 10,000.
#' @param fdr Numeric. False discovery rate. Between 0 and 1. Default value is:
#'   0.15.
#' @param bootstrap Numeric. Confidence threshold for tree bootstrap values.
#'   Default value is: 0.70.
#' @param group_genotype_key Matrix. Dimenions: nrow = number of unique
#'   genotypes, ncol = 2. Optional input.
#' @param test Character. Default = "both". User can supply three options:
#'   "both", "phyc", or "synchronous". Determines which test is run for binary
#'   data.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Both Synchronous Test & PhyC for discrete phenotype
#' phenotype <- hogwash::antibiotic_resistance
#' genotype <- hogwash::snp_genotype
#' tree <- hogwash::tree
#' hogwash(pheno = phenotype,
#'         geno = genotype,
#'         tree = tree)
#'
#' # Continuous Test for continuous phenotype
#' phenotype <- hogwash::growth
#' genotype <- hogwash::snp_genotype
#' tree <- hogwash::tree
#' hogwash(pheno = phenotype,
#'         geno = genotype,
#'         tree = tree)
#'
#' # Continuous Test while grouping SNPs into genes
#' phenotype <- hogwash::growth
#' genotype <- hogwash::snp_genotype
#' tree <- hogwash::tree
#' key <- hogwash::snp_gene_key
#' hogwash(pheno = phenotype,
#'         geno = genotype,
#'         tree = tree,
#'         group_genotype_key = key)
#'
#' # Both Synchronous Test & PhyC while grouping SNPs into genes
#' phenotype <- hogwash::antibiotic_resistance
#' genotype <- hogwash::snp_genotype
#' tree <- hogwash::tree
#' key <- hogwash::snp_gene_key
#' hogwash(pheno = phenotype,
#'         geno = genotype,
#'         tree = tree,
#'         group_genotype_key = key)
#' }
#'
#' @author Katie Saund
#'
#' @references Farhat MR, Shapiro BJ, Kieser KJ, et al. Genomic analysis
#'  identifies targets of convergent positive selection in drug-resistant
#'  \emph{Mycobacterium tuberculosis}. Nat Genet. 2013;45(10):1183–1189.
#'  doi:10.1038/ng.2747
#'
hogwash <- function(pheno,
                    geno,
                    tree,
                    file_name = "hogwash",
                    dir = ".",
                    perm = 10000,
                    fdr = 0.15,
                    bootstrap = 0.70,
                    group_genotype_key = NULL,
                    test = "both"){
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
  args$test <- test
  args$discrete_or_continuous <-
    check_input_format(args$phenotype,
                       args$tree,
                       args$genotype,
                       args$output_name,
                       args$output_dir,
                       args$perm,
                       args$fdr,
                       args$bootstrap_cutoff,
                       args$gene_snp_lookup,
                       args$test)
  select_test_type(args)
}
