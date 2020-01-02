run_plots <- function() {

# Continuous grouped
args                        <- NULL
args$test                   <- FALSE
args$phenotype              <- hogwash::growth
args$tree                   <- hogwash::tree
args$tree$node.label[1]     <- 0
args$tree$node.label        <- as.numeric(args$tree$node.label)
args$genotype               <- hogwash::snp_genotype
args$output_name            <- "dummy"
args$output_dir             <- "."
args$perm                   <- 1000
args$fdr                    <- 0.95
args$discrete_or_continuous <- "continuous"
args$bootstrap_cutoff       <- 0.7
args$group_genotype         <- TRUE
args$gene_snp_lookup        <- hogwash::snp_gene_key
run_continuous(args)

# continuous not grouped
args                        <- NULL
args$test                   <- FALSE
args$phenotype              <- hogwash::growth
args$tree                   <- hogwash::tree
args$tree$node.label[1]     <- 0
args$tree$node.label        <- as.numeric(args$tree$node.label)
args$genotype               <- hogwash::snp_genotype
args$output_name            <- "dummy"
args$output_dir             <- "."
args$perm                   <- 1000
args$fdr                    <- 1.0
args$discrete_or_continuous <- "continuous"
args$bootstrap_cutoff       <- 0.7
args$group_genotype         <- FALSE
args$gene_snp_lookup        <- NULL
run_continuous(args)


# synchronous
# Grouped
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
args$fdr                    <- 0.95
args$discrete_or_continuous <- "discrete"
args$bootstrap_cutoff       <- 0.7
args$group_genotype         <- TRUE
args$gene_snp_lookup        <- hogwash::snp_gene_key
run_synchronous(args)

# not grouped
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
args$fdr                    <- 0.95
args$discrete_or_continuous <- "discrete"
args$bootstrap_cutoff       <- 0.7
args$group_genotype         <- FALSE
args$gene_snp_lookup        <- NULL
run_synchronous(args)

# PhyC ungrouped
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
args$fdr                    <- 0.95
args$discrete_or_continuous <- "discrete"
args$bootstrap_cutoff       <- 0.7
args$group_genotype         <- FALSE
args$gene_snp_lookup        <- NULL
run_phyc(args)

# phyc grouped
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
args$fdr                    <- 0.95
args$discrete_or_continuous <- "discrete"
args$bootstrap_cutoff       <- 0.7
args$group_genotype         <- TRUE
args$gene_snp_lookup      <- hogwash::snp_gene_key
run_phyc(args)
}
