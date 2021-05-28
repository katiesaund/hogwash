# read_format_inputs ----------------------------------------------------------#

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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"
  temp_strain_key <- NULL

  # Test
  expect_equal(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type,
                                  temp_strain_key),
               "continuous")
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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"


  # Test - Bad phenotype
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - bad tree
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - bad genotype
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - bad name
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - bad directory
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - bad permutations
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - bad FDR
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - Bad bootstrap
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
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
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - Bad group genotype key
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
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
  temp_test <- "both"
  temp_group_method <- "AR"
  temp_tr_type <- "phylogram"

  # Test - Bad group method
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
})


test_that("check_input_format gives error when phenotype has NA", {
  temp_pheno <- hogwash::growth
  temp_pheno[1, 1] <- NA
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - NA in phenotype
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
})

test_that("check_input_format gives error when genotype has NA", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_geno[1, 1] <- NA
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - NA in genotype
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
})

test_that("check_input_format gives error when group_genotype_key has NA", {
  temp_pheno <- hogwash::growth
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key
  temp_group_genotype_key[1, 1] <- NA
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - NA in group_genotype_key
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
})

test_that("check_input_format gives error when tree has too few tips", {
  temp_tree <- hogwash::tree
  temp_geno <- hogwash::snp_genotype
  temp_tree <- ape::drop.tip(temp_tree, c("t2", "t3", "t4", "t5"))
  temp_pheno <- hogwash::growth
  temp_pheno <- temp_pheno[row.names(temp_pheno) %in% temp_tree$tip.label, , drop = FALSE]
  temp_geno <- temp_geno[row.names(temp_geno) %in% temp_tree$tip.label, , drop = FALSE]
  temp_name <- "temp_name"
  temp_dir <- "."
  temp_perm <- 10
  temp_fdr <- 0.10
  temp_bootstrap <- 0.875
  temp_group_genotype_key <- hogwash::snp_gene_key
  temp_test <- "both"
  temp_group_method <- "pre-ar"
  temp_tr_type <- "phylogram"

  # Test - NA in group_genotype_key
  expect_error(check_input_format(temp_pheno,
                                  temp_tree,
                                  temp_geno,
                                  temp_name,
                                  temp_dir,
                                  temp_perm,
                                  temp_fdr,
                                  temp_bootstrap,
                                  temp_group_genotype_key,
                                  temp_group_method,
                                  temp_test,
                                  temp_tr_type))
})
