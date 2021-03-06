context("Run phyc") #----------------------------------------------------------#

test_that("run_phyc() doesn't give any errors when given a discrete phenotype
          and snps that group into genes", {
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
  args$grouping_method        <- "post-ar"
  args$tree_type              <- "fan"

  expect_error(run_phyc(args), NA)
})

test_that("run_phyc() doesn't give any errors when given a discrete phenotype
          and snps (do not group into genes)", {
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
  args$gene_snp_lookup        <- NULL
  args$grouping_method        <- "post-ar"
  args$tree_type              <- "fan"

  expect_error(run_phyc(args), NA)
})

test_that("run_phyc() doesn't give any errors when given a discrete phenotype
          and snps (do not group into genes) but tree is a different order", {
            args                        <- NULL
            args$test                   <- FALSE
            args$phenotype              <- rbind(hogwash::antibiotic_resistance,
                                                 hogwash::antibiotic_resistance,
                                                 hogwash::antibiotic_resistance)
            args$genotype               <- rbind(hogwash::snp_genotype,
                                                 hogwash::snp_genotype,
                                                 hogwash::snp_genotype)
            row.names(args$genotype) <-
              row.names(args$phenotype) <-
              paste0("t", 1:21)
            set.seed(1)
            args$tree                   <- ape::rtree(21)
            args$tree                   <- ape::unroot(args$tree)
            args$tree$node.label        <- rep(100, 21)
            args$tree$node.label[1]     <- 0
            args$tree$node.label        <- as.numeric(args$tree$node.label)
            args$output_name            <- "dummy"
            args$output_dir             <- "."
            args$perm                   <- 1000
            args$fdr                    <- 0.15
            args$discrete_or_continuous <- "discrete"
            args$bootstrap_cutoff       <- 0.7
            args$group_genotype         <- FALSE
            args$gene_snp_lookup        <- NULL
            args$grouping_method        <- "post-ar"
            args$tree_type              <- "fan"

            expect_error(run_phyc(args), NA)
})


test_that("run_phyc() doesn't give any errors when genotype matrix not originally
          in the same order at the tree and phenotype matrix (do not group
          into genes)", {
            args                        <- NULL
            args$test                   <- FALSE
            args$phenotype              <- hogwash::antibiotic_resistance
            args$tree                   <- hogwash::tree
            args$tree$node.label[1]     <- 0
            args$tree$node.label        <- as.numeric(args$tree$node.label)
            args$genotype               <- hogwash::snp_genotype
            row.names(args$genotype)[1:3] <- c("t5", "t2", "t7")
            args$output_name            <- "dummy"
            args$output_dir             <- "."
            args$perm                   <- 1000
            args$fdr                    <- 0.15
            args$discrete_or_continuous <- "discrete"
            args$bootstrap_cutoff       <- 0.7
            args$group_genotype         <- FALSE
            args$gene_snp_lookup        <- NULL
            args$grouping_method        <- "post-ar"
            args$tree_type              <- "fan"

            expect_error(run_phyc(args), NA)
})

test_that("run_phyc() doesn't give any errors when phenotype matrix not originally
          in the same order at the tree and genotype matrix (do not group
          into genes)", {
            args                        <- NULL
            args$test                   <- FALSE
            args$phenotype              <- hogwash::antibiotic_resistance
            row.names(args$phenotype)[1:3] <- c("t5", "t2", "t7")
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
            args$grouping_method        <- "post-ar"
            args$gene_snp_lookup        <- NULL
            args$tree_type              <- "fan"

            expect_error(run_phyc(args), NA)
})

test_that("run_phyc() results are the same if the genotype or phenotype rownames
            do or do not match the tree tip label order", {
            args                        <- NULL
            args$test                   <- FALSE
            args$phenotype              <- hogwash::antibiotic_resistance
            row.names(args$phenotype)[1:3] <- c("t5", "t2", "t7")
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
            args$gene_snp_lookup        <- NULL
            args$grouping_method        <- "post-ar"
            args$tree_type              <- "fan"

            pheno_order <- run_phyc(args)


            args                        <- NULL
            args$test                   <- FALSE
            args$phenotype              <- hogwash::antibiotic_resistance
            args$tree                   <- hogwash::tree
            args$tree$node.label[1]     <- 0
            args$tree$node.label        <- as.numeric(args$tree$node.label)
            args$genotype               <- hogwash::snp_genotype
            row.names(args$genotype)[1:3] <- c("t5", "t2", "t7")
            args$output_name            <- "dummy"
            args$output_dir             <- "."
            args$perm                   <- 1000
            args$fdr                    <- 0.15
            args$discrete_or_continuous <- "discrete"
            args$bootstrap_cutoff       <- 0.7
            args$group_genotype         <- FALSE
            args$gene_snp_lookup        <- NULL
            args$grouping_method        <- "post-ar"
            args$tree_type              <- "fan"

            geno_order <- run_phyc(args)

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
            args$gene_snp_lookup        <- NULL
            args$grouping_method        <- "post-ar"
            args$tree_type              <- "fan"

            orders_match <- run_phyc(args)

            expect_identical(geno_order, orders_match)
            expect_identical(pheno_order, orders_match)
            expect_identical(pheno_order, geno_order)
})


test_that("run_phyc() gives sames results for the tree, when the tree is
          unrooted prior to running hogwash and when the tree has been midpoint
          rooted before submission.", {
            args                        <- NULL
            args$test                   <- FALSE
            args$phenotype              <- rbind(hogwash::antibiotic_resistance,
                                                 hogwash::antibiotic_resistance,
                                                 hogwash::antibiotic_resistance)
            args$genotype               <- rbind(hogwash::snp_genotype,
                                                 hogwash::snp_genotype,
                                                 hogwash::snp_genotype)
            row.names(args$genotype) <-
              row.names(args$phenotype) <-
              paste0("t", 1:21)
            set.seed(1)
            args$tree                   <- ape::rtree(21)
            args$tree                   <- ape::unroot(args$tree)
            args$tree$node.label        <- rep(100, 21)
            args$tree$node.label[1]     <- 0
            args$tree$node.label        <- as.numeric(args$tree$node.label)
            args$output_name            <- "dummy"
            args$output_dir             <- "."
            args$perm                   <- 1000
            args$fdr                    <- 0.15
            args$discrete_or_continuous <- "discrete"
            args$bootstrap_cutoff       <- 0.7
            args$group_genotype         <- FALSE
            args$gene_snp_lookup        <- NULL
            args$grouping_method        <- "post-ar"
            args$tree_type              <- "fan"

            unrooted_out <- run_phyc(args)

            args$tree <- phytools::midpoint.root(args$tree)
            # tree$tip.label != tree$edge
            midpoint_rooted_out <- run_phyc(args)


            tip_log <- args$tree$edge[, 2] <= ape::Ntip(args$tree)
            ordered_tips <- args$tree$edge[tip_log, 2]
            args$tree$tip.label <- args$tree$tip.label[ordered_tips]
            # Now tree$tip.label == tree$edge
            fixed_midpoint_rooted_out <- run_phyc(args)

            expect_identical(unrooted_out, midpoint_rooted_out)
            expect_identical(fixed_midpoint_rooted_out, midpoint_rooted_out)
})
