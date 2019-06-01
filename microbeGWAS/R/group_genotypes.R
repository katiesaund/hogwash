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
} # end build_gene_genotype_from_snps()


prepare_grouped_genotype <- function(geno, lookup, group_logical){
    genotypes_to_drop_because_not_present <- c(colnames(geno)[colSums(geno) == 0], colnames(geno)[colSums(geno) == nrow(geno)])
    genotype <- geno[ , colSums(geno) > 0] # we don't want to remove snps that are too rare or too common until snps are grouped into genes, then run on the grouped genes. But we need to remove SNPs that don't occur for ace to work.
    genotype <- genotype[ , colSums(genotype) < nrow(genotype)] # we don't want to remove snps that are too rare or too common until snps are grouped into genes, then run on the grouped genes. But we need to remove SNPs that occur in all isolates for ace to work.
    # TODO replace the magic numbers in the next four lines. Group into a function?
    gene_snp_lookup <- lookup[!(lookup[ , 1] %in% genotypes_to_drop_because_not_present), , drop = FALSE]
    gene_snp_lookup <- gene_snp_lookup[gene_snp_lookup[ , 1] %in% colnames(genotype), , drop = FALSE]
    unique_genes <- unique(gene_snp_lookup[ , 2])
    snps_per_gene <- table(gene_snp_lookup[ , 2])

    # Check and return output ----------------------------------------------------
    results <- list("snps_per_gene" = snps_per_gene,
                    "unique_genes" = unique_genes,
                    "gene_snp_lookup" = gene_snp_lookup,
                    "genotype" = genotype)
    return(results)
} # end prepare_grouped_genotype()

prepare_ungrouped_genotype <- function(geno, tr, group_logical){
  # Function description -------------------------------------------------------
  # Remove genotypes that are too common or rare for ancestral reconstruction to work.
  # Given that this genotype is not grouped return NULL for the variable snps_per_gene.
  # Keep track of which genotypes got removed (convergence_not_possible_genotypes).
  #
  # Inputs:
  # tr. Phylo. Ntip = nrow(pheno).
  # pheno. Matrix. Nrow = Ntip(tr). Ncol = 1.
  # geno. Matrix. Binary.
  # disc_cont. String. Either "discrete" or "continuous"
  #
  #
  # Outputs:
  # pheno_vector. Vector. Length = Ntip(tr).
  #
  # Check input --------------------------------------------------------------
  #
  # Function
  simplified_genotype <- reduce_redundancy(geno, tr) # Remove genotypes that are too rare or too commmon for (1) convergence to be possible and (2) for ancestral reconstruction to work
  snps_per_gene <- NULL

  # Check and return output --------------------------------------------------
  results <- list("snps_per_gene" = snps_per_gene,
                  "genotype" = simplified_genotype$mat,
                  "convergence_not_possible_genotypes" = simplified_genotype$dropped_genotype_names)
  return(results)
} # end prepare_ungrouped_genotype()

prepare_genotype <- function(group_logical, geno, tr, lookup){
  if (group_logical){
    prepped_geno <- prepare_grouped_genotype(geno, lookup)
  } else {
    prepped_geno <- prepare_ungrouped_genotype(geno, tr)
  }
  return(prepped_geno)
} # end prepare_genotype()


format_and_name_grouped_transitions <- function(genotype_transition){
  dummy <- genotype_transition
  genotype_transition <- rep(list(NULL), length(dummy))
  for (i in 1:length(dummy)){
    genotype_transition[[i]]$transition <- as.numeric(as.character(unlist((dummy[i]))))
  }
  return(genotype_transition)
} # end format_and_name_grouped_transitions()



group_genotypes <- function(tr, geno, genotype_reconstruction_and_confidence, genotype_transition_con, genotype_transition_orig, lookup, uni_genes){
  # Function description -------------------------------------------------------
  # Group genotypes into larger groups.
  # Examples: SNPs into genes or genes into pathways.
  #
  # Inputs:
  # tr. Phylo.
  # geno. Binary matrix. Nrow = Ntip(tr). Ncol = Number of genotypes.
  # genotype_reconstruction_and_confidence. List of Lists of four objects.
  #     Length of list == number of genotypes.
  #     1. $node_anc_rec. Reconstructed values correspond to each ancestral node.
  #     2. $tip_and_node_recon. Reconstructed values corresponding to tips and nodes.
  #     3. $tip_and_node_rec_conf. Ancestral reconstruction confidence values corresponding to tips and nodes.
  #     4. $recon_edge_mat. Ancestral reconstruction structured as a matrix with each row correspoding to a tree edge. 1st column is parent node. 2nd column is child node.
  # genotype_transition_con. List of lists.
  #     Length of list = number of genotypes.
  #     1. $transition. Length == Nedge(tr).
  #     2. $trans_dir. Length == Nedge(tr).
  # genotype_transition_orig. List of lists.
  #     Length of list = number of genotypes.
  #     1. $transition. Length == Nedge(tr).
  #     2. $trans_dir. Length == Nedge(tr).
  # lookup. Matrix. Characters. Nrow = number of genotypes.
  # uni_genes. Character vector. Equivalent to unique(lookup[ , 2]).
  #
  # Outputs:
  # "geno_recon_ordered_by_edges" = geno_recon_ordered_by_edges.
  # "geno_conf_ordered_by_edges" = geno_conf_ordered_by_edges.
  # "geno_trans_concomitant" = genotype_transition_con.
  # "geno_trans_original" = genotype_transition_orig.
  # "convergence_not_possible_genotypes" = simplified_genotype$dropped_genotype_names.
  # "genotype" = geno.
  #
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_dimensions(geno, exact_rows = Ntip(tr), min_rows = Ntip(tr), exact_cols = NULL, min_cols = 1)
  if (length(genotype_reconstruction_and_confidence) != ncol(geno)){
    stop("Need a reconstruction/confidence list for each genotype in the matrix.")
  }
  if (length(genotype_reconstruction_and_confidence[[1]]$node_anc_rec) != Nnode(tr)){
    stop("Node genotype reconstruction must be length of Nnode(tree)")
  }
  if (length(genotype_reconstruction_and_confidence[[1]]$tip_and_node_recon) != c(Nnode(tr) + Ntip(tr))){
    stop("Tip & Node genotype reconstruction must be length of Nnode(tree) + Ntip(tree)")
  }
  if (length(genotype_reconstruction_and_confidence[[1]]$tip_and_node_rec_conf) != c(Nnode(tr) + Ntip(tr))){
    stop("Tip & Node genotype confidence must be length of Nnode(tree) + Ntip(tree)")
  }
  if (nrow(genotype_reconstruction_and_confidence[[1]]$recon_edge_mat) != Nedge(tr)){
    stop("Genotype reconstruction matrix must have a row for each tree edge.")
  }
  if (length(genotype_transition_con) != ncol(geno)){
    stop("Genotype transitions for concomitant test should correspond to each genotype")
  }
  if (length(genotype_transition_con[[1]]$transition) != Nedge(tr)){
    stop("Genotype transitions should correspond to each edge on the tree.")
  }
  if (length(genotype_transition_orig) != ncol(geno)){
    stop("Genotype transitions for original test should correspond to each genotype")
  }
  if (length(genotype_transition_orig[[1]]$transition) != Nedge(tr)){
    stop("Genotype transitions should correspond to each edge on the tree.")
  }
  check_dimensions(lookup, exact_rows = ncol(geno), min_rows = ncol(geno), exact_cols = 2, min_cols = 2)

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------


  geno_recon_and_confidence_tip_node <- build_gene_anc_recon_and_conf_from_snp(tr, geno, genotype_reconstruction_and_confidence, lookup)
  genotype_transition_con <- build_gene_trans_from_snp_trans(tr, geno, genotype_transition_con, lookup)
  genotype_transition_orig    <- build_gene_trans_from_snp_trans(tr, geno, genotype_transition_orig, lookup)

  # make new geno (just at the tips, from the snps)
  geno <- build_gene_genotype_from_snps(geno, lookup)
  simplified_genotype <- reduce_redundancy(geno, tr) # Remove genotypes that are too rare or too commmon for (1) convergence to be possible and (2) for ancestral reconstruction to work
  geno <- simplified_genotype$mat
  genes_to_keep_in_consideration <- !(uni_genes %in% simplified_genotype$dropped_genotype_names)

  # remove redundancy from geno trans, geno_recon_and_confdience_tip_node_recon, and node_confidence
  genotype_transition_con <- genotype_transition_con[genes_to_keep_in_consideration]
  genotype_transition_orig  <- genotype_transition_orig[genes_to_keep_in_consideration]

  genotype_transition_con <- format_and_name_grouped_transitions(genotype_transition_con)
  genotype_transition_orig  <- format_and_name_grouped_transitions(genotype_transition_orig)

  geno_recon_and_confidence_tip_node_recon      <- geno_recon_and_confidence_tip_node$tip_node_recon[genes_to_keep_in_consideration]
  geno_recon_and_confidence_tip_node_confidence <- geno_recon_and_confidence_tip_node$tip_node_conf[ genes_to_keep_in_consideration]
  geno_conf_ordered_by_edges <- geno_recon_ordered_by_edges <- rep(list(0), ncol(geno))
  for (k in 1:ncol(geno)){
    geno_recon_ordered_by_edges[[k]] <- reorder_tips_and_nodes_to_edges(geno_recon_and_confidence_tip_node_recon[[k]],      tr)
    geno_conf_ordered_by_edges[[k]]  <- reorder_tips_and_nodes_to_edges(geno_recon_and_confidence_tip_node_confidence[[k]], tr)
  }
  results <- list("geno_recon_ordered_by_edges" = geno_recon_ordered_by_edges,
                  "geno_conf_ordered_by_edges" = geno_conf_ordered_by_edges,
                  "geno_trans_concomitant" = genotype_transition_con,
                  "geno_trans_original" = genotype_transition_orig,
                  "convergence_not_possible_genotypes" = simplified_genotype$dropped_genotype_names,
                  "genotype" = geno)
  return(results)
} # end group_genotypes()
