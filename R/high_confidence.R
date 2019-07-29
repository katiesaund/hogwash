#' discretize_conf_with_cutoff
#'
#' @description Given a vector with values that describe confidence, binarize
#'  the vector a accoriding to a cutoff value.
#' @param confidence_vector Numeric vector.
#' @param threshold Number.
#'
#' @return Confidence vector. Binary vector.
#'
#' @noRd
#'
discretize_conf_with_cutoff <- function(confidence_vector, threshold){
  # Check inputs ---------------------------------------------------------------
  check_is_number(threshold)
  if (!is.vector(confidence_vector)) {
    stop("Input must be a numeric vector")
  }
  check_is_number(confidence_vector[1])

  # Function -------------------------------------------------------------------
  confidence_vector[confidence_vector  < threshold] <- 0
  confidence_vector[confidence_vector >= threshold] <- 1

  # Check and return output ----------------------------------------------------
  check_if_binary_vector(confidence_vector)
  return(confidence_vector)
} # end discretize_conf_with_cutoff()


#' report_num_high_confidence_trans_edge
#'
#' @description Given a genotype for which you have: a list of vectors that
#'  stores if there is a genotype transition or not for each edge
#'  (genotype_transition), a list of vectors that stores if that edge is high
#'  confidence or not (high_conf_edges), and a character vector of the genotype
#'  names -- create an object that stores the number of high confidence
#'  transition edges per genotype.
#'
#' @param genotype_transition List of numeric vectors. Number of lists == number
#'  of genotypes. Length of vector == Nedge(tr).
#' @param high_conf_edges Binary vector. List of numeric vectors. Number of
#'   lists == number of genotypes. Length of vector == Nedge(tr).
#' @param geno_names Character vector. Length == ncol(genotype_matrix).
#'
#' @return num_high_conf_trans_edges. Numeric vector. Count of number
#'   of high confidence transitions per genotype. Vector is named with genotype
#'   names.
#' @noRd
#'
report_num_high_confidence_trans_edge <- function(genotype_transition,
                                                  high_conf_edges,
                                                  geno_names){
  # Check input ----------------------------------------------------------------
  check_equal(length(genotype_transition), length(geno_names))
  check_equal(length(high_conf_edges), length(geno_names))
  if (!is.vector(genotype_transition[[1]]$transition)) {
    stop("Input must have a vector.")
  }
  if (!is.vector(high_conf_edges[[1]])) {
    stop("Input must have a vector.")
  }
  if (!is.character(geno_names[1])) {
    stop("Input must be a character vector.")
  }

  # Function -------------------------------------------------------------------
  num_high_conf_trans_edges <- rep(0, length(high_conf_edges))
  for (p in 1:length(high_conf_edges)) {
    num_high_conf_trans_edges[p] <-
      sum(genotype_transition[[p]]$transition * high_conf_edges[[p]])
  }
  names(num_high_conf_trans_edges) <- geno_names

  # Return output --------------------------------------------------------------
  return(num_high_conf_trans_edges)
} # end report_num_high_confidence_trans_edge

#' assign_high_confidence_to_transition_edges
#'
#' @description Identify all edges for which the edge is high confidence and a
#'  transition edge.
#'
#' @param tr Phylo.
#' @param all_confidence_by_edge List of vectors. Each vector is binary.
#'  Length(list) == number of genomes.
#' @param geno_trans_by_edge List of vectors. Each vector is binary.
#'  Length(list) == number of genomes.
#' @param geno Matrix. Binary.
#'
#' @return edge_confident_and_trans_edge. List of vector. Each vector is binary.
#'  Length(list) == number of genomes.
#' @noRd
#'
assign_high_confidence_to_transition_edges <- function(tr,
                                                       all_confidence_by_edge,
                                                       geno_trans_by_edge,
                                                       geno){
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  check_equal(length(geno_trans_by_edge[[1]]$transition),
              ape::Nedge(tr))
  check_for_root_and_bootstrap(tr)
  check_equal(length(all_confidence_by_edge[[1]]), ape::Nedge(tr))

  # Function -------------------------------------------------------------------
  edge_confident_and_trans_edge <- rep(list(NULL), ncol(geno))
  for (k in 1:ncol(geno)) {
    edge_confident_and_trans_edge[[k]] <-
      as.numeric( (all_confidence_by_edge[[k]] +
                    geno_trans_by_edge[[k]]$transition) == 2)
  }

  # Return output --------------------------------------------------------------
  return(edge_confident_and_trans_edge)
} # end assign_high_confidence_to_transition_edges()


#' prepare_high_confidence_objects
#'
#' @description Identify high confidence edges (considering: tree bootstrap
#'  values, phenotype reconstruction, tree edge lengths, and ancestral
#'  reconstruction of genotype).
#'
#' @param genotype_transition List of lists. Number of lists = number of
#'  genotypes. Each list is made of a $transition and $trans_dir list.
#'  Length(transition) == Length(trans_dir) == Nedge(tree)
#' @param tr Phylo.
#' @param pheno_tip_node_recon_conf List of confidence values. Binary.
#'  Length(list) == Ntip() + Nedge()
#' @param boot_threshold Numeric. Between 0 and 1.
#' @param geno Matrix. Binary. Nrow = Ntip(tree). Ncol = Number of genotypes.
#' @param geno_conf_edge List of lists. Binary. Number of lists = number of
#'  genotypes. Length(each individual list) == Nedge(Tree)
#' @param geno_recon_edge List of lists. Binary. Number of lists = number of
#'  genotypes. Length(each individual list) == Nedge(Tree)
#' @param snps_in_each_gene Either Null or named table where names are genotypes
#'  and the values are number of not-yet-grouped-genotypes that go into the
#'  grouped genotype.
#'
#' @return List of objects.
#'  * dropped_genotypes Character vector. Names of the genotypes not being kept.
#'  * hi_confidence_transition_edge. List.
#'  * genotype. Matrix.
#'  * snps_per_gene. Either Null or named table where names are genotypes and
#'       the Values are number of not-yet-grouped-genotypes that go into the
#'       grouped genotype.
#'  * genotype_transition Object with two lists: $trans_dir and $transition.
#'  * geno_recon_edge. List of lists. Binary. Number of lists = number of
#'      genotypes. Length(each individual list) == Nedge(tree).
#'  * high_conf_ordered_by_edges. List.
#'  * num_high_conf_trans_edges. Numeric vector. Count of number of
#'      high confidence transitions per genotype. Vector is named with genotype
#'      names.
#'  * tr_and_pheno_hi_conf. Vector. Binary. Length = Nedge(tree).
#'
#' @noRd
#'
prepare_high_confidence_objects <- function(genotype_transition,
                                            tr,
                                            pheno_tip_node_recon_conf,
                                            boot_threshold,
                                            geno,
                                            geno_conf_edge,
                                            geno_recon_edge,
                                            snps_in_each_gene){
  # Check input ----------------------------------------------------------------
  check_equal(length(genotype_transition), ncol(geno))
  check_equal(length(genotype_transition[[1]]$transition), ape::Nedge(tr))
  check_for_root_and_bootstrap(tr)
  check_equal(length(pheno_tip_node_recon_conf),
              c(ape::Ntip(tr) + ape::Nnode(tr)))
  check_num_between_0_and_1(boot_threshold)
  check_dimensions(geno,
                   exact_rows = ape::Ntip(tr),
                   min_rows = ape::Ntip(tr),
                   exact_cols = NULL,
                   min_cols = 1)
  check_equal(length(geno_conf_edge), ncol(geno))
  check_equal(length(geno_conf_edge[[1]]), ape::Nedge(tr))
  check_equal(length(geno_recon_edge), ncol(geno))
  check_equal(length(geno_recon_edge[[1]]), ape::Nedge(tr))
  check_if_binary_vector(geno_conf_edge[[1]])
  check_if_binary_vector(geno_recon_edge[[1]])
  # TODO add checks for snps_in_each_gene

  # Function -------------------------------------------------------------------
  pheno_conf_ordered_by_edges <-
    reorder_tip_and_node_to_edge(pheno_tip_node_recon_conf, tr)
  tree_conf <- get_bootstrap_confidence(tr, boot_threshold)
  tree_conf_ordered_by_edges <- reorder_tip_and_node_to_edge(tree_conf, tr)
  short_edges <- identify_short_edges(tr)

  high_confidence_edges <-
    pheno_conf_ordered_by_edges + tree_conf_ordered_by_edges + short_edges == 3
  check_equal(length(high_confidence_edges), ape::Nedge(tr))
  all_high_confidence_edges <- rep(list(0), ncol(geno))

  # ADD IN GENO RECONSTRUCTION CONFIDENCE
  for (k in 1:ncol(geno)) {
    all_high_confidence_edges[[k]] <-
      as.numeric(geno_conf_edge[[k]] + high_confidence_edges == 2)
  }
  only_high_conf_geno_trans <-
    assign_high_confidence_to_transition_edges(tr,
                                               all_high_confidence_edges,
                                               genotype_transition,
                                               geno)
  for (i in 1:ncol(geno)) {
    genotype_transition[[i]]$transition <- only_high_conf_geno_trans[[i]]
    genotype_transition[[i]]$trans_dir <-
      only_high_conf_geno_trans[[i]] * genotype_transition[[i]]$trans_dir
  }
  names(only_high_conf_geno_trans) <- colnames(geno)
  geno_trans_by_edge <-
    report_num_high_confidence_trans_edge(genotype_transition,
                                          all_high_confidence_edges,
                                          colnames(geno))

  # KEEP ONLY GENOTYPES WITH AT LEAST TWO HIGH CONFIDENCE TRANSITION EDGES
  geno_to_keep <- keep_two_plus_hi_conf_tran_ed(genotype_transition,
                                                all_high_confidence_edges)
  genotype_transition <- genotype_transition[geno_to_keep]
  geno_recon_edge <- geno_recon_edge[geno_to_keep]
  high_conf_ordered_by_edges <- all_high_confidence_edges[geno_to_keep]
  dropped_genotypes <- get_dropped_genotypes(geno, geno_to_keep)
  geno <- geno[, geno_to_keep, drop = FALSE]
  snps_in_each_gene <-
    snps_in_each_gene[names(snps_in_each_gene) %in% colnames(geno)]

  # Return output --------------------------------------------------------------
  results <-
    list("dropped_genotypes" = dropped_genotypes,
          "hi_confidence_transition_edge" = only_high_conf_geno_trans,
          "genotype" = geno,
          "snps_per_gene" = snps_in_each_gene,
          "genotype_transition" = genotype_transition,
          "geno_recon_edge" = geno_recon_edge,
          "high_conf_ordered_by_edges" = high_conf_ordered_by_edges,
          "num_high_conf_trans_edges" = geno_trans_by_edge,
          "tr_and_pheno_hi_conf" = high_confidence_edges)
  return(results)
} # end prepare_high_confidence_objects()
