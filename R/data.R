#' Binary phenotype
#'
#' An example binary phenotype matrix.
#' @format A matrix with 7 rows and 1 columns. Each row name corresponds to a
#'   tree tip.
#' \describe{
#'   \item{resistance}{Antibiotic resistance encoded as 1; antibiotic
#'   susceptible encoded as 0}
#'  }
"antibiotic_resistance"

#' Continuous phenotype
#'
#' An example continous phenotype matrix.
#' @format A matrix with 7 rows and 1 columns. Each row name corresponds to a
#'   tree tip.
#' \describe{
#'   \item{growth}{Dummy growth rates}
#'  }
"growth"

#' Key
#'
#' An example mapping SNPs to genes.
#' @format A matrix with 8 rows and 2 columns.
#' \describe{
#'   \item{SNP}{SNP identifiers}
#'   \item{GENE}{Gene names}
#'  }
"snp_gene_key"

#' Binary genotype
#'
#' An example binary genotype matrix.
#' @format A matrix with 7 rows and 8 columns. Rows correspond to isolates and
#'   tree tips. Columns correspond to genotypes.
#' \describe{
#'   \item{SNP1}{SNP1 present encoded as 1, absent encoded as 0}
#'  }
"snp_genotype"

#' Phylogenetic tree
#'
#' An example tree with 7 tips.
#'
#' @format A rooted phylogenetic tree with 7 tips and 6 internal nodes.
#' \describe{
#'   \item{tip.labels}{Names of each tree tip}
#'   \item{node.lables}{Bootstrap values for each node}
#'  }
"tree"
