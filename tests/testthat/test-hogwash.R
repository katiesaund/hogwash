# hogwash ----------------------------------------------------------------------
test_that("hogwash doesn't give errors when run on discrete data (grouping)", {
  expect_error(hogwash(pheno = hogwash::antibiotic_resistance,
                       geno = hogwash::snp_genotype,
                       tree = hogwash::tree,
                       file_name = "dummy",
                       dir = ".",
                       perm = 10,
                       fdr = 0.15,
                       bootstrap = 0.70,
                       group_genotype_key = hogwash::snp_gene_key,
                       strain_key = NULL), NA)
})

test_that("hogwash doesn't give errors when run on discrete data (no grouping)", {
  expect_error(hogwash(pheno = hogwash::antibiotic_resistance,
                       geno = hogwash::snp_genotype,
                       tree = hogwash::tree,
                       file_name = "dummy",
                       dir = ".",
                       perm = 10,
                       fdr = 0.15,
                       bootstrap = 0.70,
                       group_genotype_key = NULL,
                       strain_key = NULL), NA)
})

test_that("hogwash doesn't give errors when run on continuous data (grouping)", {
  expect_error(hogwash(pheno = hogwash::growth,
                       geno = hogwash::snp_genotype,
                       tree = hogwash::tree,
                       file_name = "dummy",
                       dir = ".",
                       perm = 10,
                       fdr = 0.15,
                       bootstrap = 0.70,
                       group_genotype_key = hogwash::snp_gene_key,
                       strain_key = NULL), NA)
})

test_that("hogwash doesn't give errors when run on continuous data (no grouping)", {
  expect_error(hogwash(pheno = hogwash::growth,
                       geno = hogwash::snp_genotype,
                       tree = hogwash::tree,
                       file_name = "dummy",
                       dir = ".",
                       perm = 10,
                       fdr = 0.15,
                       bootstrap = 0.70,
                       group_genotype_key = NULL,
                       strain_key = NULL), NA)
})


test_that("hogwash gives errors when missing critical data", {
  # No phenotype
  expect_error(hogwash(pheno = NULL,
                       geno = hogwash::snp_genotype,
                       tree = hogwash::tree,
                       file_name = "dummy"))
  # No genotype
  expect_error(hogwash(pheno = hogwash::growth,
                       geno = NULL,
                       tree = hogwash::tree,
                       file_name = "dummy"))
  # No tree
  expect_error(hogwash(pheno = hogwash::growth,
                       geno = hogwash::snp_genotype,
                       tree = NULL,
                       file_name = "dummy"))
})
