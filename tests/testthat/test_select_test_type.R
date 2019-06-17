context("Select test type") #--------------------------------------------------#

# test select_test_type() ------------------------------------------------------
test_that("select_test_type doesn't throw error when given valid continuous input", {
  # Set up
  args                        <- NULL
  args$test                   <- FALSE
  args$phenotype              <- hogwash::growth
  args$tree                   <- hogwash::tree
  args$tree$node.label[1]     <- 0
  args$tree$node.label        <- as.numeric(args$tree$node.label)
  args$genotype               <- hogwash::snp_genotype
  args$output_name            <- "dummy_grouped"
  args$output_dir             <- "."
  args$perm                   <- 10
  args$fdr                    <- 0.15
  args$annotation             <- NULL
  args$discrete_or_continuous <- "continuous"
  args$bootstrap_cutoff       <- 0.7
  args$group_genotype         <- TRUE
  if (args$group_genotype ) {
    args$gene_snp_lookup      <- hogwash::snp_gene_key
  } else {
    args$gene_snp_lookup <- NULL
  }

  # Test
  expect_error(select_test_type(args), NA)
})

test_that("select_test_type doesn't throw error given valid discrete input", {
  # Set up
  args                        <- NULL
  args$test                   <- FALSE
  args$phenotype              <- hogwash::antibiotic_resistance
  args$tree                   <- hogwash::tree
  args$tree$node.label[1]     <- 0
  args$tree$node.label        <- as.numeric(args$tree$node.label)
  args$genotype               <- hogwash::snp_genotype
  args$output_name            <- "dummy_grouped"
  args$output_dir             <- "."
  args$perm                   <- 10
  args$fdr                    <- 0.15
  args$annotation             <- NULL
  args$discrete_or_continuous <- "discrete"
  args$bootstrap_cutoff       <- 0.7
  args$group_genotype         <- FALSE
  args$gene_snp_lookup <- NULL
  # Test
  expect_error(select_test_type(args), NA)
})



test_that("select_test_type throws error when given invalid input", {
  # Set up
  input <- NULL
  input$discrete_or_continuous <- "fake"

  # Test
  expect_error(select_test_type(input))
  expect_error(select_test_type("fake"))
})
