gene_test_from_snps <- function(){
  # goal of this is to rebuild the gene test from snps and keep snp annotation
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
} # end gene_test_from_snps

build_gene_anc_recon_and_conf_from_snp <- function(tr, geno, g_reconstruction_and_confidence, gene_to_snp_lookup_table){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  # VALIDATE INPUTS -----------------------------------------------------------#
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  check_dimensions(geno, Ntip(tr), 2, NULL, 2)
  if (length(g_reconstruction_and_confidence) != ncol(geno)){
    stop("wrong input")
  }
  check_dimensions(gene_to_snp_lookup_table, NULL, 1, 2, 2)
  check_if_binary_vector_numeric(g_reconstruction_and_confidence[[1]]$tip_and_node_recon)
  check_if_binary_vector_numeric(g_reconstruction_and_confidence[[1]]$tip_and_node_rec_conf)

  # FUNCTION ------------------------------------------------------------------#
  tip_nodes_by_snp_mat_recon <- tip_nodes_by_snp_mat_confi <- matrix(0, nrow = (Nnode(tr) + Ntip(tr)), ncol = ncol(geno))
  if (nrow(tip_nodes_by_snp_mat_recon) != length(g_reconstruction_and_confidence[[1]]$tip_and_node_recon)){
    stop("mismatch in size")
  }
  for (k in 1:ncol(geno)){
    tip_nodes_by_snp_mat_recon[ , k] <- g_reconstruction_and_confidence[[k]]$tip_and_node_recon
    tip_nodes_by_snp_mat_confi[ , k] <- g_reconstruction_and_confidence[[k]]$tip_and_node_rec_conf
  }

  row.names(tip_nodes_by_snp_mat_recon) <- row.names(tip_nodes_by_snp_mat_confi) <- c(1:nrow(tip_nodes_by_snp_mat_recon))
  colnames(tip_nodes_by_snp_mat_recon) <- colnames(tip_nodes_by_snp_mat_confi) <- colnames(geno)

  if (nrow(gene_to_snp_lookup_table) != ncol(tip_nodes_by_snp_mat_recon)){
    stop("mismatch")
  }

  recon_times_confi <- tip_nodes_by_snp_mat_recon * tip_nodes_by_snp_mat_confi
  tip_nodes_by_snp_mat_recon_with_gene_id <- rbind(tip_nodes_by_snp_mat_recon, unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))
  tip_nodes_by_snp_mat_confi_with_gene_id <- rbind(tip_nodes_by_snp_mat_confi, unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))
  recon_times_confi_with_gene_id          <- rbind(recon_times_confi,          unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))

  if (nrow(tip_nodes_by_snp_mat_recon_with_gene_id) != (nrow(tip_nodes_by_snp_mat_recon) + 1)){
    stop("rbind didn't work")
  }
  unique_genes <- unique(gene_to_snp_lookup_table[ , 2])

  gene_mat_built_from_snps <- gene_presence_confidence <- gene_absence_confidence <- matrix(0, nrow = nrow(tip_nodes_by_snp_mat_recon), ncol = length(unique_genes))
  for (j in 1:length(unique_genes)){

    # Matrix of just the SNPs found in gene "j"
    temp_mat               <- tip_nodes_by_snp_mat_recon_with_gene_id[1:(nrow(tip_nodes_by_snp_mat_recon_with_gene_id) - 1), tip_nodes_by_snp_mat_recon_with_gene_id[nrow(tip_nodes_by_snp_mat_recon_with_gene_id), ] == unique_genes[j], drop = FALSE]
    temp_recon_times_confi <-          recon_times_confi_with_gene_id[1:(nrow(recon_times_confi_with_gene_id)          - 1),           recon_times_confi_with_gene_id[nrow(recon_times_confi_with_gene_id)        , ] == unique_genes[j], drop = FALSE]
    temp_conf              <- tip_nodes_by_snp_mat_confi_with_gene_id[1:(nrow(tip_nodes_by_snp_mat_confi_with_gene_id) - 1), tip_nodes_by_snp_mat_confi_with_gene_id[nrow(tip_nodes_by_snp_mat_confi_with_gene_id), ] == unique_genes[j], drop = FALSE]
    class(temp_mat) <- class(temp_recon_times_confi) <- class(temp_conf) <- "numeric"
    for (r in 1:nrow(gene_presence_confidence)){
      if (rowSums(temp_mat)[r] == 0 & rowSums(temp_conf)[r] > 0){
        gene_absence_confidence[r, j] <- 1
      }
    }
    gene_mat_built_from_snps[ , j] <- rowSums(temp_mat)
    gene_presence_confidence[ , j] <- rowSums(temp_recon_times_confi)
  }

  gene_mat_built_from_snps <- gene_mat_built_from_snps > 0
  gene_presence_confidence <- gene_presence_confidence > 0
  class(gene_mat_built_from_snps) <- class(gene_presence_confidence) <- "numeric"

  gene_all_confidence <- gene_presence_confidence + gene_absence_confidence

  colnames(gene_mat_built_from_snps)  <- colnames(gene_all_confidence)  <- unique_genes
  row.names(gene_mat_built_from_snps) <- row.names(gene_all_confidence) <- c(1:nrow(gene_mat_built_from_snps))

  gene_list_built_from_snps <- gene_conf_list_built_from_snps<- rep(list(0), length(unique_genes))
  for (m in 1:length(unique_genes)){
    gene_list_built_from_snps[[m]]      <- gene_mat_built_from_snps[ , m, drop = TRUE]
    gene_conf_list_built_from_snps[[m]] <- gene_all_confidence[      , m, drop = TRUE]
  }
  names(gene_list_built_from_snps) <- names(gene_conf_list_built_from_snps) <- unique_genes
  return(list("tip_node_recon" = gene_list_built_from_snps, "tip_node_conf" = gene_conf_list_built_from_snps ))
} # end build_gene_anc_recon_from_snp()


build_node_anc_recon_from_gene_list <- function(gene_list, tr){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  # TODO test this function works.
  gene_list_built_from_snps_just_node_anc_rec <- rep(list(0), length(gene_list$tip_and_node_recon))
  for (m in 1:length(gene_list$tip_and_node_recon)){
    gene_list_built_from_snps_just_node_anc_rec[[m]] <- gene_list[[m]]$tip_and_node_recon[(Ntip(tr) + 1):(Ntip(tr) + Nedge(tr))]
  }
  return(gene_list_built_from_snps_just_node_anc_rec)
} # end build_node_anc_recon_from_gene_list()


build_gene_trans_from_snp_trans <- function(tr, geno, geno_transition, gene_to_snp_lookup_table){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------

  # VALIDATE INPUTS -----------------------------------------------------------#
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  check_dimensions(geno, Ntip(tr), 2, NULL, 2)
  if (length(geno_transition) != ncol(geno)){
    stop("wrong input")
  }
  check_dimensions(gene_to_snp_lookup_table, NULL, 1, 2, 2)
  check_if_binary_vector_numeric(geno_transition[[1]]$transition)

  # FUNCTION ------------------------------------------------------------------#

  edges_by_snp_mat <- matrix(0, nrow = Nedge(tr), ncol = ncol(geno))
  if (nrow(edges_by_snp_mat) != length(geno_transition[[1]]$transition)){
    stop("mismatch in size")
  }
  for (k in 1:ncol(geno)){
    edges_by_snp_mat[ , k] <- geno_transition[[k]]$transition
  }
  row.names(edges_by_snp_mat) <- c(1:nrow(edges_by_snp_mat))
  colnames(edges_by_snp_mat) <- colnames(geno)

  if (nrow(gene_to_snp_lookup_table) != ncol(edges_by_snp_mat)){
    stop("mismatch")
  }

  edges_by_snp_mat_with_gene_id <- rbind(edges_by_snp_mat, unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))
  if (nrow(edges_by_snp_mat_with_gene_id) != (nrow(edges_by_snp_mat) + 1)){
    stop("rbind didn't work")
  }

  unique_genes <- unique(gene_to_snp_lookup_table[ , 2])
  gene_mat_built_from_snps <- matrix(0, nrow = nrow(edges_by_snp_mat), ncol = length(unique_genes))
  for (j in 1:length(unique_genes)){
    temp_mat <- edges_by_snp_mat_with_gene_id[1:(nrow(edges_by_snp_mat_with_gene_id) - 1) , edges_by_snp_mat_with_gene_id[nrow(edges_by_snp_mat_with_gene_id), ] == unique_genes[j], drop = FALSE]
    class(temp_mat) <- "numeric"
    temp_column <- rowSums(temp_mat)
    gene_mat_built_from_snps[ , j] <- temp_column
  }

  gene_mat_built_from_snps <- gene_mat_built_from_snps > 0
  class(gene_mat_built_from_snps) <- "numeric"

  colnames(gene_mat_built_from_snps) <- unique_genes
  row.names(gene_mat_built_from_snps) <- c(1:nrow(gene_mat_built_from_snps))

  gene_list_built_from_snps <- rep(list(0), length(unique_genes))
  for (m in 1:length(unique_genes)){
    gene_list_built_from_snps[[m]] <- gene_mat_built_from_snps[ , m, drop = TRUE]
  }
  names(gene_list_built_from_snps) <- unique_genes

  return(gene_list_built_from_snps)
} # end build_gene_trans_from_snp_trans()


build_gene_genotype_from_snps <- function(geno, gene_to_snp_lookup_table){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  unique_genes <- unique(gene_to_snp_lookup_table[ , 2])
  samples_by_genes <- matrix(0, nrow = nrow(geno), ncol = length(unique_genes))
  colnames(samples_by_genes) <- unique_genes
  row.names(samples_by_genes) <- row.names(geno)


  snp_geno_with_gene_id <- rbind(geno, unlist(gene_to_snp_lookup_table[ , 2, drop = TRUE]))
  if (nrow(snp_geno_with_gene_id) != (nrow(geno) + 1)){
    stop("rbind didn't work")
  }

  for (j in 1:length(unique_genes)){
    temp_mat <- snp_geno_with_gene_id[1:(nrow(snp_geno_with_gene_id) - 1) , snp_geno_with_gene_id[nrow(snp_geno_with_gene_id), ] == unique_genes[j], drop = FALSE]
    class(temp_mat) <- "numeric"
    temp_column <- rowSums(temp_mat)
    samples_by_genes[ , j] <- temp_column
  }

  samples_by_genes <- samples_by_genes > 0
  class(samples_by_genes) <- "numeric"
  return(samples_by_genes)
}
