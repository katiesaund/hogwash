context("Read in and format inputs") #-----------------------------------------#

# test check_input_format
test_that("check_input_format works for valid inputs", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key),
               NA)
})

test_that("check_input_format gives error for invalid input", {
  temp_pheno <- letters[1:10]
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key))
})

test_that("check_input_format gives error for invalid input", {
  temp_pheno <- hogwash::growth
  temp_tree <- "tree"
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key))
})

test_that("check_input_format gives error for invalid input", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- 1:10
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key))
})

test_that("check_input_format gives error for invalid input", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- matrix(NA, 10, 10)
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key))
})

test_that("check_input_format gives error for invalid input", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "/this/is/a/fake/dir/"
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key))
})

test_that("check_input_format gives error for invalid input", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 0.10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key))
})

test_that("check_input_format gives error for invalid input", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key))
})

test_that("check_input_format gives error for invalid input", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 10
  temp_group_genotype_key <- hogwash::snp_gene_key

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key))
})

test_that("check_input_format gives error for invalid input", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- letters[1:10]

  # Test
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key))
})
