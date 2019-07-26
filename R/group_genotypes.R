#' build_gene_anc_recon_and_conf_from_snp
#'
#' @param tr Phylo.
#' @param geno Matrix. Binary. Genotypes in columns, isolates in rows.
#' @param g_recon_and_conf List of 2 objects: $tip_and_node_recon and
#'   $tip_and_node_rec_conf. Length of both objects is number of genotypes.
#'   Length of vectors within each object is number of tips + number of nodes
#'   in the tree. Both objects are binary (1/0).
#' @param gene_to_snp_lookup_table Matrix. Two columns. 1st column is the less
#'  condensed grouping, 2nd column is the condensed groups.
#'
#' @return List of two objects:
#'   * tip_node_recon. List of the genotypes built into groups. Length = number
#'        of unique grouped genotypes.
#'   * tip_node_conf. List of genotype confidences built into groups. Length =
#'       number of unique grouped genotypes.
#'
#' @noRd
#'
build_gene_anc_recon_and_conf_from_snp <- function(tr,
                                                   geno,
                                                   g_recon_and_conf,
                                                   gene_to_snp_lookup_table){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  check_dimensions(geno, ape::Ntip(tr), 2, NULL, 2)
  check_equal(length(g_recon_and_conf), ncol(geno))
  check_dimensions(gene_to_snp_lookup_table, NULL, 1, 2, 2)
  check_if_binary_vector_numeric(g_recon_and_conf[[1]]$tip_and_node_recon)
  check_if_binary_vector_numeric(g_recon_and_conf[[1]]$tip_and_node_rec_conf)

  # Function -------------------------------------------------------------------
  tip_nodes_by_snp_mat_recon <-
    tip_nodes_by_snp_mat_confi <-
    matrix(0, nrow = (ape::Nnode(tr) + ape::Ntip(tr)), ncol = ncol(geno))
  if (nrow(tip_nodes_by_snp_mat_recon) !=
      length(g_recon_and_conf[[1]]$tip_and_node_recon)) {
    stop("mismatch in size")
  }
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

  if (nrow(gene_to_snp_lookup_table) != ncol(tip_nodes_by_snp_mat_recon)) {
    stop("mismatch")
  }

  recon_times_confi <- tip_nodes_by_snp_mat_recon * tip_nodes_by_snp_mat_confi
  tip_nod_by_mat_recon_w_gene_id <-
    rbind(tip_nodes_by_snp_mat_recon,
          unlist(gene_to_snp_lookup_table[, 2, drop = TRUE]))
  tip_node_by_snp_mat_conf_w_id <-
    rbind(tip_nodes_by_snp_mat_confi,
          unlist(gene_to_snp_lookup_table[, 2, drop = TRUE]))
  recon_times_confi_with_gene_id <-
    rbind(recon_times_confi,
          unlist(gene_to_snp_lookup_table[, 2, drop = TRUE]))

  if (nrow(tip_nod_by_mat_recon_w_gene_id) !=
      (nrow(tip_nodes_by_snp_mat_recon) + 1)) {
    stop("rbind didn't work")
  }
  unique_genes <- unique(gene_to_snp_lookup_table[, 2])

  gene_mat_built_from_snps <-
    gene_presence_confidence <-
    gene_absence_confidence <-
    matrix(0,
           nrow = nrow(tip_nodes_by_snp_mat_recon),
           ncol = length(unique_genes))

  for (j in 1:length(unique_genes)) {
    # Matrix of just the SNPs found in gene "j"
    temp_mat <-
      tip_nod_by_mat_recon_w_gene_id[1:(
        nrow(tip_nod_by_mat_recon_w_gene_id) - 1),
        tip_nod_by_mat_recon_w_gene_id[nrow(
          tip_nod_by_mat_recon_w_gene_id), ] == unique_genes[j],
        drop = FALSE]
    temp_recon_times_confi <-
      recon_times_confi_with_gene_id[1:(
        nrow(recon_times_confi_with_gene_id) - 1),
        recon_times_confi_with_gene_id[nrow(
          recon_times_confi_with_gene_id), ] == unique_genes[j],
        drop = FALSE]
    temp_conf <-
      tip_node_by_snp_mat_conf_w_id[1:(
        nrow(tip_node_by_snp_mat_conf_w_id) - 1),
        tip_node_by_snp_mat_conf_w_id[nrow(
          tip_node_by_snp_mat_conf_w_id), ] == unique_genes[j],
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
  class(gene_mat_built_from_snps) <-
    class(gene_presence_confidence) <-
    "numeric"

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
} # end build_gene_anc_recon_from_snp()

#' build_gene_trans_from_snp_trans
#'
#' @description Calculate which edges of the tree are transitions for the
#'  grouped genotype from transition information of the ungrouped genotype. The
#'  function title reflects one such type of grouping: snps into genes.
#'
#' @param tr Phylo.
#' @param geno Matrix. Nrow = Ntip(tr). Ncol = Number of ungrouped genotypes.
#'   Binary.
#' @param geno_transition List of named vectors. Length == number of ungrouped
#'   genotypes. $transition. Length = Nedge(tr). Values or 0 or 1.
#'   $trans_dir. Length = Nedge(tr). Values -1, 0, or 1.
#' @param gene_to_snp_lookup_table Matrix. Ncol = 2. Nrow = number of ungrouped
#'   genotypes.
#'
#' @return List of two objects. $transition and $trans_dir.
#'
#' @noRd
#'
build_gene_trans_from_snp_trans <- function(tr,
                                            geno,
                                            geno_transition,
                                            gene_to_snp_lookup_table){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_dimensions(geno,
                   exact_rows = ape::Ntip(tr),
                   min_rows = 1,
                   exact_cols =  nrow(gene_to_snp_lookup_table),
                   min_cols = 1)
  check_if_binary_matrix(geno)
  check_equal(length(geno_transition), ncol(geno))
  check_equal(length(geno_transition[[1]]$transition), ape::Nedge(tr))
  check_equal(length(geno_transition[[1]]$trans_dir), ape::Nedge(tr))
  check_if_binary_vector_numeric(geno_transition[[1]]$transition)
  check_dimensions(gene_to_snp_lookup_table,
                   exact_rows = ncol(geno),
                   min_rows = 1,
                   exact_cols = 2,
                   min_cols = 1)

  # Function -------------------------------------------------------------------
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

  trans_ed_by_snp_mat_w_id <-
    rbind(transition_edges_by_snp_mat,
          unlist(gene_to_snp_lookup_table[, 2, drop = TRUE]))
  trans_dir_ed_by_snp_mat_w_id <-
    rbind(trans_dir_edges_by_snp_mat,
          unlist(gene_to_snp_lookup_table[, 2, drop = TRUE]))

  unique_genes <- unique(gene_to_snp_lookup_table[, 2])
  gene_trans_mat_built_from_snp <-
    gene_trans_dir_mat_blt_frm_snp <-
    matrix(0,
           nrow = nrow(transition_edges_by_snp_mat),
           ncol = length(unique_genes))

  for (j in 1:length(unique_genes)) {
    temp_mat <-
      trans_ed_by_snp_mat_w_id[1:(
        nrow(trans_ed_by_snp_mat_w_id) - 1),
        trans_ed_by_snp_mat_w_id[nrow(
          trans_ed_by_snp_mat_w_id), ] == unique_genes[j], drop = FALSE]
    class(temp_mat) <- "numeric"
    temp_column <- rowSums(temp_mat)
    gene_trans_mat_built_from_snp[, j] <- temp_column

    temp_dir_mat <-
      trans_dir_ed_by_snp_mat_w_id[1:(
        nrow(trans_dir_ed_by_snp_mat_w_id) - 1),
        trans_dir_ed_by_snp_mat_w_id[nrow(
          trans_dir_ed_by_snp_mat_w_id), ] == unique_genes[j], drop = FALSE]
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

  gene_transition_list_built_from_snps <- rep(list(0), length(unique_genes))
  for (m in 1:length(unique_genes)) {
    gene_transition_list_built_from_snps[[m]] <-
      gene_trans_mat_built_from_snp[, m, drop = TRUE]
  }
  names(gene_transition_list_built_from_snps) <- unique_genes

  gene_trans_dir_list_built_from_snps <- rep(list(0), length(unique_genes))
  for (m in 1:length(unique_genes)) {
    gene_trans_dir_list_built_from_snps[[m]] <-
      gene_trans_dir_mat_blt_frm_snp[, m, drop = TRUE]
  }
  names(gene_trans_dir_list_built_from_snps) <- unique_genes

  temp_results <- rep(list(list()), length(unique_genes))
  for (i in 1:length(unique_genes)) {
    temp_results[[i]]$transition <-
      unname(unlist(gene_transition_list_built_from_snps[i]))
    temp_results[[i]]$trans_dir <-
      unname(unlist(gene_trans_dir_list_built_from_snps[i]))
  }
  # Return output --------------------------------------------------------------
  return(temp_results)
} # end build_gene_trans_from_snp_trans()


#' build_gene_genotype_from_snps
#'
#' @description Build presence/absence of the grouped genotypes (e.g. gene) from
#'  the presence/absence of the ungrouped genotypes (e.g. snps).
#'
#' @param geno Matrix. Columns are genotypes. Rows are isolates.
#' @param gene_to_snp_lookup_table Matrix. 2 columns. 1st column are the
#'  ungrouped genotypes 2nd column are the grouped genotypes.
#'
#' @return samples_by_genes. Matrix.
#'
#' @noRd
#'
build_gene_genotype_from_snps <- function(geno, gene_to_snp_lookup_table){
  # Check input ----------------------------------------------------------------
  # TODO add input checks
  # Function -------------------------------------------------------------------
  unique_genes <- unique(gene_to_snp_lookup_table[, 2])
  samples_by_genes <- matrix(0, nrow = nrow(geno), ncol = length(unique_genes))
  colnames(samples_by_genes) <- unique_genes
  row.names(samples_by_genes) <- row.names(geno)

  snp_geno_with_gene_id <-
    rbind(geno, unlist(gene_to_snp_lookup_table[, 2, drop = TRUE]))
  if (nrow(snp_geno_with_gene_id) != (nrow(geno) + 1)) {
    stop("rbind didn't work")
  }

  for (j in 1:length(unique_genes)) {
    temp_mat <-
      snp_geno_with_gene_id[1:(
        nrow(snp_geno_with_gene_id) - 1),
        snp_geno_with_gene_id[nrow(
          snp_geno_with_gene_id), ] == unique_genes[j], drop = FALSE]
    class(temp_mat) <- "numeric"
    temp_column <- rowSums(temp_mat)
    samples_by_genes[, j] <- temp_column
  }

  samples_by_genes <- samples_by_genes > 0
  class(samples_by_genes) <- "numeric"
  # Return output --------------------------------------------------------------
  return(samples_by_genes)
} # end build_gene_genotype_from_snps()

#' prepare_grouped_genotype
#'
#' @description Remove genotypes that are too common (omnipresent) or rare
#'  (completely absent) for ancestral reconstruction to work from both the
#'  genotype and the lookup key.
#'
#' @param geno Matrix. Binary. Nrow = Ntip(tr). Ncol = number of original
#'  genotypes.
#' @param lookup Matrix. Ncol = 2. Nrow = number of genotypes.
#'
#' @return List of four objects:
#'  * $snp_per_gene. Named table. Names are genotypes. Values are number of
#'      not-yet-grouped-genotypes that go into the grouped genotype.
#'  * $unique_genes. Character vector. Vector of genotype names.
#'  * $gene_snp_lookup. Character matrix. Ncol = 2. Nrow = number of genotypes
#'      that are neither omni-present or completely absent.
#'  * $genotype. Matrix.
#' @noRd
#'
prepare_grouped_genotype <- function(geno, lookup){
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
  # SNPs that don't occur for ace to work.
  genotype <- geno[, colSums(geno) > 0]
  genotype <- genotype[, colSums(genotype) < nrow(genotype)]

  # TODO replace the magic numbers in the next four lines.
  # TODO Group into a function?
  gene_snp_lookup <-
    lookup[!(lookup[, 1] %in% geno_to_drop_bc_not_present),
           , drop = FALSE]
  gene_snp_lookup <-
    gene_snp_lookup[gene_snp_lookup[, 1] %in% colnames(genotype),
                    , drop = FALSE]
  unique_genes <- unique(gene_snp_lookup[, 2])
  snps_per_gene <- table(gene_snp_lookup[, 2])

  # Check and return output ----------------------------------------------------
  results <- list("snps_per_gene" = snps_per_gene,
                  "unique_genes" = unique_genes,
                  "gene_snp_lookup" = gene_snp_lookup,
                  "genotype" = genotype)
  return(results)
} # end prepare_grouped_genotype()

#' prepare_ungrouped_genotype
#'
#' @description Remove genotypes that are too common or rare for ancestral
#'  reconstruction to work. Given that this genotype is not grouped return NULL
#'  for the variable snps_per_gene. Keep track of which genotypes got removed in
#'  no_convergence_genotypes.
#'
#' @param geno Matrix. Binary. Ncol = number of genotypes. Nrow = Ntip(tr).
#' @param tr Phylo.
#'
#' @return List of three objects:
#'  * $snp_per_gene. Object of NULL.
#'  * $genotype. Matrix.
#'  * $convergene_not_possible_genotypes. Character vector. Vector of genotype
#'       names.
#' @noRd
prepare_ungrouped_genotype <- function(geno, tr){
  # Check input ----------------------------------------------------------------
  check_dimensions(geno, exact_rows = ape::Ntip(tr), min_rows = 1, min_cols = 1)
  check_if_binary_matrix(geno)
  check_for_root_and_bootstrap(tr)

  # Function -------------------------------------------------------------------
  simple_geno <- remove_rare_or_common_geno(geno, tr)
  snps_per_gene <- NULL

  # Check and return output --------------------------------------------------
  results <-
    list("snps_per_gene" = snps_per_gene,
         "genotype" = simple_geno$mat,
         "no_convergence_genotypes" = simple_geno$dropped_genotype_names)
  return(results)
} # end prepare_ungrouped_genotype()

#' prepare_genotype
#'
#' @description Funnel the genotype to be prepared for downstream use. The
#'  preparation depends on if the genotype is going to be grouped or not.
#' @param group_logical Logical.
#' @param geno Genotype matrix. Binary. Nrow = Ntip(tr). Ncol = number of
#'   ungrouped genotypes.
#' @param tr Phylo.
#' @param lookup Either NULL or a Matrix. Ncol = 2.
#'
#' @return  prepped_geno. List of multiple objects. Content depends on value of
#'  group_logical. For grouped genotypes output includes:
#'  * $snp_per_gene. Named table. Names are genotypes. Values are number of
#'      not-yet-grouped-genotypes that go into the grouped genotype.
#'  * $unique_genes. Character vector. Vector of genotype names.
#'  * $gene_snp_lookup. Character matrix. Ncol = 2. Nrow = number of genotypes
#'      that are neither omni-present or completely absent.
#'  * $genotype. Matrix.
#'
#'  For not grouped genotypes output includes:
#'  * $snp_per_gene. Object of NULL.
#'  * $genotype. Matrix.
#'  * $convergene_not_possible_genotypes. Character vector. Vector of genotype
#'      names.
#' @noRd
prepare_genotype <- function(group_logical, geno, tr, lookup){
  # Check input ----------------------------------------------------------------
  if (!is.logical(group_logical)) {
    stop("Input must be either TRUE or FALSE")
  }
  check_for_root_and_bootstrap(tr)
  if (!is.null(lookup)) {
    check_dimensions(lookup, exact_cols = 2, min_cols = 2, min_rows = 1)
  }
  check_dimensions(geno, exact_rows = ape::Ntip(tr), min_rows = 1, min_cols = 1)
  #
  # Function -------------------------------------------------------------------
  if (group_logical) {
    prepped_geno <- prepare_grouped_genotype(geno, lookup)
  } else {
    prepped_geno <- prepare_ungrouped_genotype(geno, tr)
  }
  return(prepped_geno)
} # end prepare_genotype()

#' group_genotypes
#'
#' @description Group genotypes into larger groups. Examples: SNPs into genes or
#'  genes into pathways.
#'
#' @param tr Phylo.
#' @param geno Binary matrix. Nrow = Ntip(tr). Ncol = Number of genotypes.
#' @param geno_recon_and_conf List of Lists of four objects.
#'   Length of list == number of genotypes.
#'   * $node_anc_rec. Reconstructed values correspond to each ancestral node.
#'   * $tip_and_node_recon. Reconstructed values corresponding to tips and
#'       nodes.
#'   * $tip_and_node_rec_conf. Ancestral reconstruction confidence values
#'       corresponding to tips and nodes.
#'   * $recon_edge_mat. Ancestral reconstruction structured as a matrix with
#'       each row correspoding to a tree edge. 1st column is parent node. 2nd
#'       column is child node.
#' @param genotype_transition_sync List of lists. Length of list = number of
#'   genotypes.
#'   * $transition. Length == Nedge(tr).
#'   * $trans_dir. Length == Nedge(tr).
#' @param genotype_transition_phyc List of lists. Length of list = number of
#'  genotypes.
#'  * $transition. Length == Nedge(tr).
#'  * $trans_dir. Length == Nedge(tr).
#' @param lookup Matrix. Characters. Nrow = number of genotypes.
#' @param uni_genes Character vector. Equivalent to unique(lookup[,2]).
#'
#' @return List of 6 objects.
#'  * geno_recon_ordered_by_edges. List.
#'  * geno_conf_ordered_by_edges. List.
#'  * geno_trans_synchronous. List of two objects. $transition and $trans_dir.
#'  * geno_trans_phyc. List of two objects. $transition and $trans_dir.
#'  * no_convergence_genotypes. Character vector.
#'  * genotype. Matrix.
#'
#' @noRd
#'
group_genotypes <- function(tr,
                            geno,
                            geno_reconstruction_and_conf,
                            genotype_transition_sync,
                            genotype_transition_phyc,
                            lookup,
                            uni_genes){
  # TODO
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
                   exact_rows = ncol(geno),
                   min_rows = ncol(geno),
                   exact_cols = 2,
                   min_cols = 2)

  # Function -------------------------------------------------------------------
  geno_recon_and_conf_tip_node <-
    build_gene_anc_recon_and_conf_from_snp(tr,
                                           geno,
                                           geno_reconstruction_and_conf,
                                           lookup)
  genotype_transition_sync <-
    build_gene_trans_from_snp_trans(tr, geno, genotype_transition_sync, lookup)
  genotype_transition_phyc <-
    build_gene_trans_from_snp_trans(tr, geno, genotype_transition_phyc, lookup)

  # make new geno (just at the tips, from the snps)
  geno <- build_gene_genotype_from_snps(geno, lookup)
  simplified_genotype <- remove_rare_or_common_geno(geno, tr)
  geno <- simplified_genotype$mat
  genes_to_keep_in_consideration <-
    !(uni_genes %in% simplified_genotype$dropped_genotype_names)

  genotype_transition_sync <-
    genotype_transition_sync[genes_to_keep_in_consideration]
  genotype_transition_phyc <-
    genotype_transition_phyc[genes_to_keep_in_consideration]

  geno_recon_and_confidence_tip_node_recon <-
    geno_recon_and_conf_tip_node$tip_node_recon[genes_to_keep_in_consideration]
  geno_recon_and_confidence_tip_node_confidence <-
    geno_recon_and_conf_tip_node$tip_node_conf[ genes_to_keep_in_consideration]

  geno_conf_ordered_by_edges <-
    geno_recon_ordered_by_edges <-
    rep(list(0), ncol(geno))
  for (k in 1:ncol(geno)) {
    geno_recon_ordered_by_edges[[k]] <-
      reorder_tip_and_node_to_edge(
        geno_recon_and_confidence_tip_node_recon[[k]], tr)
    geno_conf_ordered_by_edges[[k]] <-
      reorder_tip_and_node_to_edge(
        geno_recon_and_confidence_tip_node_confidence[[k]], tr)
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
} # end group_genotypes()
