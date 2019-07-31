## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(hogwash)

## ---- echo = FALSE, out.width = "700px"----------------------------------
knitr::include_graphics("tree_with_nomenclature_for_wiki.png")

## ---- echo = FALSE, out.width = "700px"----------------------------------
knitr::include_graphics("ancestral_reconstruction.png")

## ---- echo = FALSE, out.width = "700px"----------------------------------
knitr::include_graphics("both_transition_edges.png")

## ---- echo = FALSE, out.width = "700px"----------------------------------
knitr::include_graphics("continuous_algorithm_demo.png")

## ---- echo = FALSE, out.width = "700px"----------------------------------
knitr::include_graphics("both_transition_edges.png")

## ------------------------------------------------------------------------
growth_phenotype <- hogwash::growth
snp_genotype <- hogwash::snp_genotype
tree <- hogwash::tree
key <- hogwash::snp_gene_key


## ------------------------------------------------------------------------
antibiotic_phenotype <- hogwash::antibiotic_resistance

## ---- echo = FALSE, out.width = "700px"----------------------------------
knitr::include_graphics("tree_with_nomenclature_for_wiki.png")

## ---- eval=FALSE---------------------------------------------------------
#  hogwash(pheno = growth_phenotype,
#          geno = snp_genotype,
#          tree = tree)

## ---- eval=FALSE---------------------------------------------------------
#  hogwash(pheno = growth_phenotype,
#          geno = snp_genotype,
#          tree = tree,
#          group_genotype_key = key)

## ---- eval=FALSE---------------------------------------------------------
#  hogwash(pheno = antibiotic_phenotype,
#          geno = snp_genotype,
#          tree = tree)

## ---- eval=FALSE---------------------------------------------------------
#  hogwash(pheno = antibiotic_phenotype,
#          geno = snp_genotype,
#          tree = tree,
#          group_genotype_key = key)
#  

## ---- echo = FALSE, out.width = "700px"----------------------------------
knitr::include_graphics("2019-06-16_fqR_roary_convergence_test.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("continuous_phenotype_ancestral_reconstruction.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("Genotype_transitions.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("delta_phenotype_by_edge_type.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("KS_null_distribution.png")

## ---- echo = FALSE, out.width = "400px"----------------------------------
knitr::include_graphics("ordered_by_tree_edge.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("synchronous_grouped_genotype_summary_heatmap.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("synchronous_phenotype_transitio.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("synchronous_genotype_transition.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("synchronous_null_distribution.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("phenotype_binary_reconstruction.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("convergence_genotype_reconstruction.png")

## ---- echo = FALSE, out.width = "500px"----------------------------------
knitr::include_graphics("convergence_null_distribution.png")

