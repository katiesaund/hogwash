context("Run synchronous") #---------------------------------------------------#

test_that("run_synchronous() doesn't give any errors when given a discrete
          phenotype and snps that group into genes", {
  args                        <- NULL
  args$test                   <- FALSE
  args$phenotype              <- hogwash::antibiotic_resistance
  args$tree                   <- hogwash::tree
  args$tree$node.label[1]     <- 0
  args$tree$node.label        <- as.numeric(args$tree$node.label)
  args$genotype               <- hogwash::snp_genotype
  args$output_name            <- "dummy"
  args$output_dir             <- "."
  args$perm                   <- 1000
  args$fdr                    <- 0.15
  args$discrete_or_continuous <- "discrete"
  args$bootstrap_cutoff       <- 0.7
  args$group_genotype         <- TRUE
  if (args$group_genotype ) {
    args$gene_snp_lookup      <- hogwash::snp_gene_key
  } else {
    args$gene_snp_lookup <- NULL
  }
  expect_error(run_synchronous(args), NA)
})

test_that("run_synchronous() doesn't give any errors when given a discrete
          phenotype and snps (do not group into genes)", {
  args                        <- NULL
  args$test                   <- FALSE
  args$phenotype              <- hogwash::antibiotic_resistance
  args$tree                   <- hogwash::tree
  args$tree$node.label[1]     <- 0
  args$tree$node.label        <- as.numeric(args$tree$node.label)
  args$genotype               <- hogwash::snp_genotype
  args$output_name            <- "dummy"
  args$output_dir             <- "."
  args$perm                   <- 1000
  args$fdr                    <- 0.15
  args$discrete_or_continuous <- "discrete"
  args$bootstrap_cutoff       <- 0.7
  args$group_genotype         <- FALSE
  args$gene_snp_lookup <- NULL
  expect_error(run_synchronous(args), NA)
})
