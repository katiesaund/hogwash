context("Run continuous") #----------------------------------------------------#

test_that("run_continuous() doesn't give any errors when given a continuous
          phenotype and snps (do not group into genes)" ,{
  args                        <- NULL
  args$test                   <- FALSE
  args$phenotype              <- hogwash::growth
  args$tree                   <- hogwash::tree
  args$tree$node.label[1]     <- 0
  args$tree$node.label        <- as.numeric(args$tree$node.label)
  args$genotype               <- hogwash::snp_genotype
  args$output_name            <- "dummy_grouped"
  args$output_dir             <- "."
  args$perm                   <- 1000
  args$fdr                    <- 0.15
  args$discrete_or_continuous <- "continuous"
  args$bootstrap_cutoff       <- 0.7
  args$group_genotype         <- FALSE
  args$gene_snp_lookup        <- NULL
  expect_error(run_continuous(args), NA)
})

test_that("run_continuous() doesn't give any errors when given a continuous
          phenotype and snps are grouped into genes" ,{
  args                        <- NULL
  args$test                   <- FALSE
  args$phenotype              <- hogwash::growth
  args$tree                   <- hogwash::tree
  args$tree$node.label[1]     <- 0
  args$tree$node.label        <- as.numeric(args$tree$node.label)
  args$genotype               <- hogwash::snp_genotype
  args$output_name            <- "dummy_grouped"
  args$output_dir             <- "."
  args$perm                   <- 1000
  args$fdr                    <- 0.15
  args$discrete_or_continuous <- "continuous"
  args$bootstrap_cutoff       <- 0.7
  args$group_genotype         <- TRUE
  if (args$group_genotype ) {
    args$gene_snp_lookup      <- hogwash::snp_gene_key
  } else {
    args$gene_snp_lookup <- NULL
  }
  expect_error(run_continuous(args), NA)
})
