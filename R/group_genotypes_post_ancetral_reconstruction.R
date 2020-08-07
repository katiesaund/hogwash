#' Build grouped ancestral reconstruction and confidence
#'
#' @param tr Phylo.
#' @param geno Matrix. Binary. Genotypes in columns, isolates in rows.
#' @param g_recon_and_conf List of 2 objects:
#'   \describe{
#'     \item{tip_and_node_recon}{Length == number of genotypes. Length of
#'     vectors within each object is number of tips + number of nodes in the
#'     tree. Both objects are binary (1/0).}
#'     \item{tip_and_node_rec_conf}{Length == number of genotypes.Length of
#'     vectors within each object is number of tips + number of nodes in the
#'     tree. Both objects are binary (1/0).}
#'   }
#' @param lookup Matrix. Two columns. 1st column are the names
#'  of the genotypes to be grouped. 2nd column contains the appropriate groups.
#'
#' @return List of two objects:
#'   \describe{
#'     \item{tip_node_recon}{List of the genotypes built into groups. Length =
#'     number of unique grouped genotypes.}
#'     \item{tip_node_conf}{List of genotype confidences built into groups.
#'     Length = number of unique grouped genotypes.}
#'   }
#' @noRd
build_gene_anc_recon_and_conf_from_snp_post_ar <- function(tr,
                                                           geno,
                                                           g_recon_and_conf,
                                                           lookup){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  check_dimensions(geno, ape::Ntip(tr), 2, NULL, 2)
  check_equal(length(g_recon_and_conf), ncol(geno))
  check_dimensions(lookup, NULL, 1, 2, 2)
  check_if_binary_vector_numeric(g_recon_and_conf[[1]]$tip_and_node_recon)
  check_if_binary_vector_numeric(g_recon_and_conf[[1]]$tip_and_node_rec_conf)

  # Function -------------------------------------------------------------------
  unique_genes <- unique(lookup[, 2])

  # Create two matrices:
  # (1) tip_nodes_by_snp_mat_recon
  # (2) tip_nodes_by_snp_mat_confi
  # Where the rows are tips then nodes, names 1:sum(Ntip + Nnode)
  # Where the columns are the genotype names (SNPs not genes yet)
  tip_nodes_by_snp_mat_recon <-
    tip_nodes_by_snp_mat_confi <-
    matrix(0, nrow = (ape::Nnode(tr) + ape::Ntip(tr)), ncol = ncol(geno))
  check_equal(nrow(tip_nodes_by_snp_mat_recon),
              length(g_recon_and_conf[[1]]$tip_and_node_recon))
  for (k in 1:ncol(geno)) {
    tip_nodes_by_snp_mat_recon[, k] <-
      g_recon_and_conf[[k]]$tip_and_node_recon
    tip_nodes_by_snp_mat_confi[, k] <-
      g_recon_and_conf[[k]]$tip_and_node_rec_conf
  }

  row.names(tip_nodes_by_snp_mat_recon) <-
    row.names(tip_nodes_by_snp_mat_confi) <-
    c(1:nrow(tip_nodes_by_snp_mat_recon))
  colnames(tip_nodes_by_snp_mat_recon) <-
    colnames(tip_nodes_by_snp_mat_confi) <-
    colnames(geno)

  # Create a thrid matrix called recon_times_confi
  # This matrix records reconstruction times confidence values
  recon_times_confi <- tip_nodes_by_snp_mat_recon * tip_nodes_by_snp_mat_confi

  # Build the three gene matrices from the SNP matrices
  # (1) Gene reconstruction (gene_mat_built_from_snps)
  # (2) Gene present and confident (gene_presence_confidence)
  # (3) Gene absent and confidence (gene_absence_confidence)
  gene_mat_built_from_snps <-
    gene_presence_confidence <-
    gene_absence_confidence <-
    matrix(0,
           nrow = nrow(tip_nodes_by_snp_mat_recon),
           ncol = length(unique_genes))

  for (j in 1:length(unique_genes)) {
    # Generate temporary matrices- subset to just the SNPs found in gene "j"
    temp_mat <-
      tip_nodes_by_snp_mat_recon[,
                                 colnames(tip_nodes_by_snp_mat_recon) %in% lookup[ , 1][lookup[, 2] == unique_genes[j]],
                                 drop = FALSE]

    temp_conf <-
      tip_nodes_by_snp_mat_confi[,
                                 colnames(tip_nodes_by_snp_mat_confi) %in% lookup[ , 1][lookup[, 2] == unique_genes[j]],
                                 drop = FALSE]
    temp_recon_times_confi <-
      recon_times_confi[,
                        colnames(recon_times_confi) %in% lookup[ , 1][lookup[, 2] == unique_genes[j]],
                        drop = FALSE]
    class(temp_mat) <-
      class(temp_recon_times_confi) <-
      class(temp_conf) <-
      "numeric"
    for (r in 1:nrow(gene_presence_confidence)) {
      if (rowSums(temp_mat)[r] == 0 & rowSums(temp_conf)[r] > 0) {
        gene_absence_confidence[r, j] <- 1
      }
    }
    gene_mat_built_from_snps[, j] <- rowSums(temp_mat)
    gene_presence_confidence[, j] <- rowSums(temp_recon_times_confi)
  }

  gene_mat_built_from_snps <- gene_mat_built_from_snps > 0
  gene_presence_confidence <- gene_presence_confidence > 0
  class(gene_mat_built_from_snps) <- class(gene_presence_confidence) <- "numeric"
  gene_all_confidence <- gene_presence_confidence + gene_absence_confidence

  colnames(gene_mat_built_from_snps) <-
    colnames(gene_all_confidence) <-
    unique_genes
  row.names(gene_mat_built_from_snps) <-
    row.names(gene_all_confidence) <-
    c(1:nrow(gene_mat_built_from_snps))

  gene_list_built_from_snps <-
    gene_conf_list_built_from_snps <-
    rep(list(0), length(unique_genes))

  for (m in 1:length(unique_genes)) {
    gene_list_built_from_snps[[m]] <-
      gene_mat_built_from_snps[, m, drop = TRUE]
    gene_conf_list_built_from_snps[[m]] <-
      gene_all_confidence[, m, drop = TRUE]
  }
  names(gene_list_built_from_snps) <-
    names(gene_conf_list_built_from_snps) <-
    unique_genes

  # Return output --------------------------------------------------------------
  return(list("tip_node_recon" = gene_list_built_from_snps,
              "tip_node_conf" = gene_conf_list_built_from_snps ))
}

#' Build grouped transition information
#'
#' @description Calculate which edges of the tree are transitions for the
#'   grouped genotype from transition information of the ungrouped genotype. The
#'   function title reflects one such type of grouping: snps into genes.
#'
#' @param tr Phylo.
#' @param geno Matrix. Nrow = Ntip(tr). Ncol = Number of ungrouped genotypes.
#'   Binary.
#' @param geno_transition List of named vectors. Length == number of ungrouped
#'   genotypes. $transition. Length = Nedge(tr). Values or 0 or 1. $trans_dir.
#'   Length = Nedge(tr). Values -1, 0, or 1.
#' @param gene_to_snp_lookup_table Matrix. Ncol = 2. Nrow = number of unique
#'   ungrouped genotype - grouped genotype pairings.
#'
#' @return List of two objects:
#'   * transition
#'   * trans_dir
#'
#' @noRd
build_gene_trans_from_snp_trans_post_ar <- function(tr,
                                                    geno,
                                                    geno_transition,
                                                    gene_to_snp_lookup_table){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_dimensions(geno,
                   exact_rows = ape::Ntip(tr),
                   min_rows = 1,
                   min_cols = 1)
  check_if_binary_matrix(geno)
  check_equal(length(geno_transition), ncol(geno))
  check_equal(length(geno_transition[[1]]$transition), ape::Nedge(tr))
  check_equal(length(geno_transition[[1]]$trans_dir), ape::Nedge(tr))
  check_if_binary_vector_numeric(geno_transition[[1]]$transition)
  check_dimensions(gene_to_snp_lookup_table,
                   min_rows = 1,
                   exact_cols = 2,
                   min_cols = 1)

  # Function -------------------------------------------------------------------
  unique_genes <- unique(gene_to_snp_lookup_table[, 2])

  # Create two matrices:
  # (1) transition_edges_by_snp_mat (Transitions)
  # (2) trans_dir_edges_by_snp_mat (Transition Direction)
  # Where the matrices have one row for each tree edge and 1 column for each
  #   ungrouped genotype.
  transition_edges_by_snp_mat <-
    trans_dir_edges_by_snp_mat <-
    matrix(0, nrow = ape::Nedge(tr), ncol = ncol(geno))

  for (k in 1:ncol(geno)) {
    transition_edges_by_snp_mat[, k] <- geno_transition[[k]]$transition
    trans_dir_edges_by_snp_mat[, k] <- geno_transition[[k]]$trans_dir
  }
  row.names(transition_edges_by_snp_mat) <-
    row.names(trans_dir_edges_by_snp_mat) <-
    c(1:nrow(transition_edges_by_snp_mat))

  colnames(transition_edges_by_snp_mat) <-
    colnames(trans_dir_edges_by_snp_mat) <-
    colnames(geno)

  # Create two gene matrices from the snp matrices:
  # (1) gene_trans_mat_built_from_snp (Transitions)
  # (2) gene_trans_dir_mat_blt_frm_snp (Transition Direction)
  gene_trans_mat_built_from_snp <-
    gene_trans_dir_mat_blt_frm_snp <-
    matrix(0,
           nrow = nrow(transition_edges_by_snp_mat),
           ncol = length(unique_genes))

  for (j in 1:length(unique_genes)) {
    temp_mat <-
      transition_edges_by_snp_mat[,
                                  colnames(transition_edges_by_snp_mat) %in% gene_to_snp_lookup_table[ , 1][gene_to_snp_lookup_table[, 2] == unique_genes[j]],
                                  drop = FALSE]
    class(temp_mat) <- "numeric"
    temp_column <- rowSums(temp_mat)
    gene_trans_mat_built_from_snp[, j] <- temp_column

    temp_dir_mat <-
      trans_dir_edges_by_snp_mat[,
                                 colnames(trans_dir_edges_by_snp_mat) %in% gene_to_snp_lookup_table[ , 1][gene_to_snp_lookup_table[, 2] == unique_genes[j]],
                                 drop = FALSE]
    class(temp_dir_mat) <- "numeric"
    temp_dir_column <- rowSums(temp_dir_mat)
    gene_trans_dir_mat_blt_frm_snp[, j] <- temp_dir_column
  }

  gene_trans_mat_built_from_snp <- gene_trans_mat_built_from_snp > 0
  class(gene_trans_mat_built_from_snp) <- "numeric"
  colnames(gene_trans_mat_built_from_snp) <- unique_genes
  row.names(gene_trans_mat_built_from_snp) <-
    c(1:nrow(gene_trans_mat_built_from_snp))

  gene_trans_dir_mat_blt_frm_snp <- gene_trans_dir_mat_blt_frm_snp > 0
  class(gene_trans_dir_mat_blt_frm_snp) <- "numeric"
  colnames(gene_trans_dir_mat_blt_frm_snp) <- unique_genes
  row.names(gene_trans_dir_mat_blt_frm_snp) <-
    c(1:nrow(gene_trans_dir_mat_blt_frm_snp))

  gene_trans_list_from_snp <- rep(list(0), length(unique_genes))
  for (m in 1:length(unique_genes)) {
    gene_trans_list_from_snp[[m]] <-
      gene_trans_mat_built_from_snp[, m, drop = TRUE]
  }
  names(gene_trans_list_from_snp) <- unique_genes

  gene_trans_dir_list_from_snp <- rep(list(0), length(unique_genes))
  for (m in 1:length(unique_genes)) {
    gene_trans_dir_list_from_snp[[m]] <-
      gene_trans_dir_mat_blt_frm_snp[, m, drop = TRUE]
  }
  names(gene_trans_dir_list_from_snp) <- unique_genes

  temp_results <- rep(list(list()), length(unique_genes))
  for (i in 1:length(unique_genes)) {
    temp_results[[i]]$transition <-
      unname(unlist(gene_trans_list_from_snp[i]))
    temp_results[[i]]$trans_dir <-
      unname(unlist(gene_trans_dir_list_from_snp[i]))
  }
  # Return output --------------------------------------------------------------
  return(temp_results)
}

#' Build grouped genotype
#'
#' @description Build presence/absence of the grouped genotypes (e.g. gene) from
#'   the presence/absence of the ungrouped genotypes (e.g. snps).
#'
#' @param geno Matrix. Columns are genotypes. Rows are isolates.
#' @param gene_to_snp_lookup_table Matrix. 2 columns. 1st column are the
#'   ungrouped genotypes 2nd column are the grouped genotypes. The table's 1st
#'   column must contain only genotypes that are found in the row.names(geno).
#'
#' @return samples_by_genes. Matrix.
#'
#' @noRd
build_gene_genotype_from_snps_post_ar <- function(geno, gene_to_snp_lookup_table){
  # Check input ----------------------------------------------------------------
  check_class(geno, "matrix")
  check_class(gene_to_snp_lookup_table, "matrix")
  if (sum(!gene_to_snp_lookup_table[ , 1] %in% colnames(geno)) > 0) {
    stop("gene_to_snp_lookup_table must only contain genotypes in geno")
  }
  check_dimensions(gene_to_snp_lookup_table, NULL, 1, 2, 2)
  check_dimensions(geno, NULL, 1, NULL, 1)

  # Function -------------------------------------------------------------------
  unique_genes <- unique(gene_to_snp_lookup_table[, 2])
  samples_by_genes <- matrix(0, nrow = nrow(geno), ncol = length(unique_genes))
  colnames(samples_by_genes) <- unique_genes
  row.names(samples_by_genes) <- row.names(geno)

  for (j in 1:length(unique_genes)) {
    temp_mat <-
      geno[, colnames(geno) %in% gene_to_snp_lookup_table[ , 1][gene_to_snp_lookup_table[, 2] == unique_genes[j]],
           drop = FALSE]
    class(temp_mat) <- "numeric"
    temp_column <- rowSums(temp_mat)
    samples_by_genes[, j] <- temp_column
  }

  samples_by_genes <- samples_by_genes > 0
  class(samples_by_genes) <- "numeric"
  # Return output --------------------------------------------------------------
  return(samples_by_genes)
}

#' Remove rare and common genotypes from grouped genotypes
#'
#' @description Remove genotypes that are too common (omnipresent) or rare
#'   (completely absent) for ancestral reconstruction to work from both the
#'   genotype and the lookup key.
#'
#' @param geno Matrix. Binary. Nrow = Ntip(tr). Ncol = number of original
#'   genotypes.
#' @param lookup Matrix. Ncol = 2. Nrow = genotypes with group assignments.
#'
#' @return List of four objects:
#'   \describe{
#'     \item{snp_per_gene.}{Named table. Names are genotypes. Values are number
#'     of not-yet-grouped-genotypes that go into the grouped genotype.}
#'     \item{unique_genes.}{Character vector. Vector of genotype names.}
#'     \item{gene_snp_lookup.}{Character matrix. Ncol = 2. Nrow = number of
#'     genotypes that are neither omni-present or completely absent.}
#'     \item{genotype.}{Matrix.}
#'   }
#' @noRd
prepare_grouped_genotype_post_ar <- function(geno, lookup){
  # Check input ----------------------------------------------------------------
  check_dimensions(lookup, min_rows = 1, exact_cols = 2, min_cols = 2)
  check_dimensions(geno, min_cols = 1, min_rows = 1)
  check_if_binary_matrix(geno)

  # Function -------------------------------------------------------------------
  geno_to_drop_bc_not_present <-
    c(colnames(geno)[colSums(geno) == 0],
      colnames(geno)[colSums(geno) == nrow(geno)])

  # we don't want to remove snps that are too rare or too common until snps are
  # grouped into genes, then run on the grouped genes. But we need to remove
  # SNPs that don't occur for ape::ace to work.
  genotype <- geno[, colSums(geno) > 0]
  genotype <- genotype[, colSums(genotype) < nrow(genotype)]

  gene_snp_lookup <-
    lookup[!(lookup[, 1] %in% geno_to_drop_bc_not_present),
           , drop = FALSE]
  gene_snp_lookup <-
    gene_snp_lookup[gene_snp_lookup[, 1] %in% colnames(genotype),
                    , drop = FALSE]
  unique_genes <- unique(gene_snp_lookup[, 2])
  snps_per_gene <- table(gene_snp_lookup[, 2])

  genotype <-
    genotype[, colnames(genotype) %in% gene_snp_lookup[, 1], drop = FALSE]
  # Check and return output ----------------------------------------------------
  if (ncol(genotype) < 2) {
    stop("There are fewer than 2 genotypes that have variable presence/absence and are named in the grouping key")
  }
  if (nrow(gene_snp_lookup) < 2) {
    stop("There are fewer than 2 genotypes that are named in the grouping key and found in the genotype matrix")
  }
  results <- list("snps_per_gene" = snps_per_gene,
                  "unique_genes" = unique_genes,
                  "gene_snp_lookup" = gene_snp_lookup,
                  "genotype" = genotype)
  return(results)
}

#' Group genotypes
#'
#' @description Group genotypes into larger groups. Examples: SNPs into genes or
#'   genes into pathways.
#'
#' @param tr Phylo.
#' @param geno Binary matrix. Nrow = Ntip(tr). Ncol = Number of genotypes.
#' @param geno_recon_and_conf List of Lists of four objects. Length of list ==
#'   number of genotypes.
#'   \describe{
#'     \item{node_anc_rec.}{Reconstructed values correspond to each ancestral
#'     node.}
#'     \item{tip_and_node_recon.}{Reconstructed values corresponding to tips and
#'     nodes.}
#'     \item{tip_and_node_rec_conf.}{Ancestral reconstruction confidence values
#'     corresponding to tips and nodes.}
#'     \item{recon_edge_mat.}{Ancestral reconstruction structured as a matrix
#'     with each row correspoding to a tree edge. 1st column is parent node. 2nd
#'     column is child node.}
#'   }
#' @param genotype_transition_sync List of lists. Length of list = number of
#'   genotypes.
#'   \describe{
#'     \item{transition}{Length == Nedge(tr).}
#'     \item{trans_dir}{Length == Nedge(tr).}
#'   }
#' @param genotype_transition_phyc List of lists. Length of list = number of
#'  genotypes.
#'   \describe{
#'     \item{transition}{Length == Nedge(tr).}
#'     \item{trans_dir}{Length == Nedge(tr).}
#'   }
#' @param lookup Matrix. Characters. Nrow = number of unique genotype-group
#'  pairings.
#' @param uni_genes Character vector. Equivalent to unique(lookup[,2]).
#'
#' @return List of 6 objects:
#'   \describe{
#'     \item{geno_recon_ordered_by_edges.}{List.}
#'     \item{geno_conf_ordered_by_edges.}{List.}
#'     \item{geno_trans_synchronous.}{List of two objects. $transition and
#'     $trans_dir.}
#'     \item{geno_trans_phyc.}{List of two objects. $transition and $trans_dir.}
#'     \item{no_convergence_genotypes.}{Character vector.}
#'     \item{genotype.}{Matrix.}
#'   }
#'
#' @noRd
group_genotypes_post_ar <- function(tr,
                            geno,
                            geno_reconstruction_and_conf,
                            genotype_transition_sync,
                            genotype_transition_phyc,
                            lookup,
                            uni_genes){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_dimensions(geno,
                   exact_rows = ape::Ntip(tr),
                   min_rows = ape::Ntip(tr),
                   exact_cols = NULL,
                   min_cols = 1)
  check_equal(length(geno_reconstruction_and_conf), ncol(geno))
  check_equal(length(geno_reconstruction_and_conf[[1]]$node_anc_rec),
              ape::Nnode(tr))
  check_equal(length(geno_reconstruction_and_conf[[1]]$tip_and_node_recon),
              c(ape::Nnode(tr) + ape::Ntip(tr)))
  check_equal(length(geno_reconstruction_and_conf[[1]]$tip_and_node_rec_conf),
              c(ape::Nnode(tr) + ape::Ntip(tr)))
  check_equal(nrow(geno_reconstruction_and_conf[[1]]$recon_edge_mat),
              ape::Nedge(tr))
  check_equal(length(genotype_transition_sync), ncol(geno))
  check_equal(length(genotype_transition_sync[[1]]$transition), ape::Nedge(tr))
  check_equal(length(genotype_transition_phyc), ncol(geno))
  check_equal(length(genotype_transition_phyc[[1]]$transition), ape::Nedge(tr))
  check_dimensions(lookup,
                   min_rows = 2,
                   exact_cols = 2,
                   min_cols = 2)
  check_is_string(uni_genes[1])
  if (length(uni_genes) > length(unique(lookup[, 2]))) {
    stop("Too many genotypes in uni_genes")
  }

  # Function -------------------------------------------------------------------
  geno_recon_and_conf_tip_node <-
    build_gene_anc_recon_and_conf_from_snp_post_ar(tr,
                                                   geno,
                                                   geno_reconstruction_and_conf,
                                                   lookup)

  genotype_transition_sync <-
    build_gene_trans_from_snp_trans_post_ar(tr,
                                            geno,
                                            genotype_transition_sync,
                                            lookup)

  genotype_transition_phyc <-
    build_gene_trans_from_snp_trans_post_ar(tr,
                                            geno,
                                            genotype_transition_phyc,
                                            lookup)

  # make new geno (just at the tips, from the snps)
  geno <- build_gene_genotype_from_snps_post_ar(geno, lookup)
  simplified_genotype <- remove_rare_or_common_geno(geno, tr)
  geno <- simplified_genotype$mat
  genes_to_keep_in_consideration <-
    !(uni_genes %in% simplified_genotype$dropped_genotype_names)

  genotype_transition_sync <-
    genotype_transition_sync[genes_to_keep_in_consideration]
  genotype_transition_phyc <-
    genotype_transition_phyc[genes_to_keep_in_consideration]

  geno_recon_conf_tip_node_recon <-
    geno_recon_and_conf_tip_node$tip_node_recon[genes_to_keep_in_consideration]
  geno_rec_conf_tip_node_conf <-
    geno_recon_and_conf_tip_node$tip_node_conf[ genes_to_keep_in_consideration]

  geno_conf_ordered_by_edges <-
    geno_recon_ordered_by_edges <-
    rep(list(0), ncol(geno))
  for (k in 1:ncol(geno)) {
    geno_recon_ordered_by_edges[[k]] <-
      reorder_tip_and_node_to_edge(
        geno_recon_conf_tip_node_recon[[k]], tr)
    geno_conf_ordered_by_edges[[k]] <-
      reorder_tip_and_node_to_edge(
        geno_rec_conf_tip_node_conf[[k]], tr)
  }

  # Return output --------------------------------------------------------------
  results <-
    list(
      "geno_recon_ordered_by_edges" = geno_recon_ordered_by_edges,
      "geno_conf_ordered_by_edges" = geno_conf_ordered_by_edges,
      "geno_trans_synchronous" = genotype_transition_sync,
      "geno_trans_phyc" = genotype_transition_phyc,
      "no_convergence_genotypes" = simplified_genotype$dropped_genotype_names,
      "genotype" = geno)
  return(results)
}
