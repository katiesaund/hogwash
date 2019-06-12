# Katie Saund

# This is an implementation of phyC, originally described by Farhat et al. Nature Genetics 2013.
# This script will:
# 1. Read in data:
#    a. 1 phenotype (either continuous or discrete),
#    b. multiple binary genotypes,
#    c. and a phylogenetic tree.
# 2. Perform ancestral reconstructions of that phenotype and those genotypes over the tree.
# 3. Perform a test to associate the presence/absence of the genotypes with the phenotype:
#    a. co-occurrence of transitions of reconstructed genotype/phenotype:
#       i. Significantly more (or fewer) edges when both genotype and phenotype transition from one character state to another.
# 4. Return p-values for all genotypes tested and a list of significant hits, if any, after FDR correction.
# 5. Return figures of all genotype & phenotype reconstructions and transition trees.

# LOAD PHYC LIBRARY OF FUNCTIONS ----------------------------------------------#
library(microbeGWAS)

# READ IN ARGUMENTS -----------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE) # Grab arguments from the PBS script
args <- read_in_arguments(args)
select_test_type(args)

# END OF PHYC -----------------------------------------------------------------#
