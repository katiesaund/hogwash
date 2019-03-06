run_phyc <- function(args){
  # FORMAT INPUTS ---------------------------------------------------------------#
  results_object <- NULL
  results_object$log <- capture.output(sessionInfo()) # log session info

  if (!args$built_from_snps) {
    simplified_genotype <- reduce_redundancy(args$genotype, args$tree, args$output_dir, args$output_name) # Remove genotypes that are too rare or too commmon for (1) convergence to be possible and (2) for ancestral reconstruction to work
    genotype <- simplified_genotype$mat
    results_object$convergence_not_possible_genotypes <- simplified_genotype$dropped_genotype_names
  } else {
    genotypes_to_drop_because_not_present <- colnames(args$genotype)[colSums(args$genotype) == 0]
    genotype <- args$genotype[ , colSums(args$genotype) > 0] # we don't want to remove snps that are too rare or too common until snps are grouped into genes, then run on the grouped genes. But we need to remove SNPs that don't occur for ace to work.
    gene_snp_lookup <- args$gene_snp_lookup[!(args$gene_snp_lookup[ , 1] %in% genotypes_to_drop_because_not_present), , drop = FALSE]
  }

  phenotype_vector <- convert_matrix_to_vector(args$phenotype) # TODO add check that it's possible to have phenotype convergence

  # ---------------------------------------------------------------------------#
  # PHYC
  # ---------------------------------------------------------------------------#

  # ANCESTRAL RECONSTRUCTION OF PHENOTYPE -------------------------------------#
  set.seed(18591124)
  pheno_recon_and_conf  <- ancestral_reconstruction_by_ML(args$tree, args$phenotype, 1, args$discrete_or_continuous, args$bootstrap_cutoff)
  tree_conf             <- get_bootstrap_confidence(args$tree, args$bootstrap_cutoff)
  pheno_trans           <- identify_transition_edges(args$tree, args$phenotype, 1, pheno_recon_and_conf$node_anc_rec, args$discrete_or_continuous)
  pheno_recon_edge_mat  <- pheno_recon_and_conf$recon_edge_mat
  short_edges           <- identify_short_edges(args$tree)
  pheno_conf_ordered_by_edges  <- reorder_tips_and_nodes_to_edges(pheno_recon_and_conf$tip_and_node_rec_conf, args$tree)
  pheno_recon_ordered_by_edges <- reorder_tips_and_nodes_to_edges(pheno_recon_and_conf$tip_and_node_recon, args$tree)
  tree_conf_ordered_by_edges   <- reorder_tips_and_nodes_to_edges(tree_conf, args$tree)
  # pheno_trans objects and short_edges are already ordered by edges.

  # ANCESTRAL RECONSTRUCTION OF GENOTYPES ---------------------------------------#
  # INIATLIAZE DATA STRUCTS
  geno_recon_and_conf <- geno_trans <- geno_conf_ordered_by_edges <- geno_recon_ordered_by_edges <- rep(list(0), ncol(genotype))

  # PERFORM ANCESTRAL RECONSTRUCTION
  for(k in 1:ncol(genotype)){
    geno_recon_and_conf[[k]] <- ancestral_reconstruction_by_ML(args$tree, genotype, k, "discrete", args$bootstrap_cutoff)
  }

  if (args$built_from_snps){
    # TODO change this if statement into a function
    # CONVERT SNPS INTO GENES HERE

    # tip_and_node_ancestral_reconstruction
    geno_recon_and_conf$tip_and_node_recon <- build_gene_anc_recon_from_snp(args$tree, genotype, geno_recon_and_conf, gene_snp_lookup)
    # node_ancestral_reconstruction
    print("A")
    geno_recon_and_conf$node_anc_rec <- build_node_anc_recon_from_gene_list(geno_recon_and_conf, args$tree) # must be after build_gene_anc_recon_from_snp so tip_and_node_recon is updated to be gene not snps
    # tip and node confidence in ancestral reconstruction
    print("B")
    geno_recon_and_conf$tip_and_node_rec_conf <- build_gene_confidence_from_snp(geno_recon_and_conf, args$tree, gene_snp_lookup, genotype)
    print("C")
    # tip_nodes_by_snp_mat <- matrix(0, nrow = (Nnode(args$tree) + Ntip(args$tree)), ncol = ncol(genotype))
    #
    # if (nrow(tip_nodes_by_snp_mat) != length(geno_recon_and_conf[[k]]$tip_and_node_recon)){
    #   stop("mismatch in size")
    # }
    #
    # for (k in 1:ncol(genotype)){
    #   tip_nodes_by_snp_mat[ , k] <- geno_recon_and_conf[[k]]$tip_and_node_recon
    # }
    # row.names(tip_nodes_by_snp_mat) <- c(1:nrow(tip_nodes_by_snp_mat))
    # colnames(tip_nodes_by_snp_mat) <- colnames(genotype)
    #
    # if (gene_snp_lookup[ , 1, drop = TRUE] != colnames(tip_nodes_by_snp_mat)){
    #   stop("gene lookup size mismatch")
    # }
    #
    # tip_nodes_by_snp_mat_with_gene_id <- rbind(tip_nodes_by_snp_mat, unlist(gene_snp_lookup[ , 2, drop = TRUE]))
    # if (nrow(tip_nodes_by_snp_mat_with_gene_id) != (nrow(tip_nodes_by_snp_mat) + 1)){
    #   stop("rbind didn't work")
    # }
    #
    # unique_genes <- unique(gene_snp_lookup[ , 2])
    # gene_mat_built_from_snps <- matrix(0, nrow = nrow(tip_nodes_by_snp_mat), ncol = length(unique_genes))
    # for (j in 1:length(unique_genes)){
    #   temp_mat <- tip_nodes_by_snp_mat_with_gene_id[1:(nrow(tip_nodes_by_snp_mat_with_gene_id) -1) , tip_nodes_by_snp_mat_with_gene_id[nrow(tip_nodes_by_snp_mat_with_gene_id), ] == unique_genes[j], drop = FALSE]
    #   class(temp_mat) <- "numeric"
    #   temp_column <- rowSums(temp_mat)
    #   gene_mat_built_from_snps[ , j] <- temp_column
    # }
    #
    # gene_mat_built_from_snps <- gene_mat_built_from_snps > 0
    # class(gene_mat_built_from_snps) <- "numeric"
    #
    # colnames(gene_mat_built_from_snps) <- unique_genes
    # row.names(gene_mat_built_from_snps) <- c(1:nrow(gene_mat_built_from_snps))
    #
    # gene_list_built_from_snps <- rep(list(0), length(unique_genes))
    # for (m in 1:length(unique_genes)){
    #   gene_list_built_from_snps[[m]] <- gene_mat_built_from_snps[ , m, drop = TRUE]
    # }
    # names(gene_list_built_from_snps) <- unique_genes

    # end ancestral reconstruction of nodes plus tips

    # TODO create a test to check/viz that I did the above assignments correctly and started from the correct piece of data.
    # TODO how do I reconcile geno_recon_and_conf[[]] with node_anc_rec and tip_and_node_rec_conf given the new tip_and_node_recon?
      # todo repeat a similar process with confidence and just the node anc reconstruction

    #gene_list_built_from_snps_just_node_anc_rec <- rep(list(0), length(unique_genes))
    #for (m in 1:length(unique_genes)){
    #  gene_list_built_from_snps_just_node_anc_rec[[m]] <- gene_list_built_from_snps[[m]][(Ntip(args$tree) + 1):(Ntip(args$tree) + Nedge(args$tree))]
    #}
  }


  # IDENTIFY TRANSITION EDGES AND REFORMAT
  for (k in 1:ncol(genotype)){
    geno_trans[[k]] <- identify_transition_edges(args$tree, genotype, k, geno_recon_and_conf[[k]]$node_anc_rec, "discrete")
    geno_conf_ordered_by_edges[[k]]  <- reorder_tips_and_nodes_to_edges(geno_recon_and_conf[[k]]$tip_and_node_rec_conf, args$tree)
    geno_recon_ordered_by_edges[[k]] <- reorder_tips_and_nodes_to_edges(geno_recon_and_conf[[k]]$tip_and_node_recon,    args$tree)
  }

  # IDENTIFY HIGH CONFIDENCE EDGES (BOOTSTRAP, PHENOTYPE RECON, LENGTH, GENOTYPE RECONSTRUCTION)
  # TREE BOOTSTRAP, PHENOTYPUE RECONSTRUCTION CONFIDENCE, AND EDGE LENGTHS
  high_confidence_edges <- pheno_conf_ordered_by_edges + tree_conf_ordered_by_edges + short_edges == 3
  all_high_confidence_edges <- rep(list(0), ncol(genotype))

  # ADD IN GENOTYPE RECONSTRUCTION CONFIDENCE
  for (k in 1:ncol(genotype)){
    all_high_confidence_edges[[k]] <- as.numeric(geno_conf_ordered_by_edges[[k]] + high_confidence_edges == 2)
  }

  # SAVE FILE WITH NUMBER OF HIGH CONFIDENCE TRANSITION EDGES PER GENOTYPE-----#
  results_object$high_confidence_trasition_edges <- high_confidence_edges
  num_high_confidence_trasition_edges <- report_num_high_confidence_trans_edge(geno_trans, all_high_confidence_edges, colnames(genotype), args$output_dir, args$output_name)
  results_object$num_high_confidence_trasition_edges <- num_high_confidence_trasition_edges

  # KEEP ONLY GENOTYPES WITH AT LEAST TWO HIGH CONFIDENCE TRANSITION EDGES ----#
  geno_to_keep                  <- keep_at_least_two_high_conf_trans_edges(geno_trans, all_high_confidence_edges)
  geno_recon_ordered_by_edges   <- geno_recon_ordered_by_edges[geno_to_keep]
  high_conf_ordered_by_edges    <- all_high_confidence_edges[geno_to_keep]
  geno_trans                    <- geno_trans[geno_to_keep]

  dropped_genotypes <- get_dropped_genotypes(genotype, geno_to_keep)
  results_object$dropped_genotypes <- dropped_genotypes

  genotype                      <- genotype[ , geno_to_keep, drop = FALSE]

  # break following if else into two seperate functions
  if (args$discrete_or_continuous == "continuous"){
    # RUN PERMUTATION TEST ------------------------------------------------------#
    results_all_transitions <- calculate_genotype_significance(genotype, args$perm, geno_trans, args$tree, pheno_trans, high_conf_ordered_by_edges, geno_recon_ordered_by_edges)

    # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
    corrected_pvals_all_transitions <- get_sig_hits_while_correcting_for_multiple_testing(results_all_transitions$pvals, args$alpha)

    # SUBSET SIGNIFICANT HITS SO MEDIAN(DELTA PHENOTYPE) ON TRANSITION EDGES > MEDIAN(DELTA PHENOTYPE) ALL EDGES
    all_transitions_sig_hits <- keep_hits_with_more_change_on_trans_edges(results_all_transitions, corrected_pvals_all_transitions, args$alpha)

    # SAVE AND PLOT RESULTS -----------------------------------------------------#
    trans_mat_results <- plot_significant_hits("continuous", args$tree, args$alpha, args$output_dir, args$output_name, corrected_pvals_all_transitions, phenotype_vector, args$annotation, args$perm, results_all_transitions, pheno_recon_and_conf$node_anc_rec, geno_recon_ordered_by_edges, high_conf_ordered_by_edges, geno_trans, genotype, pheno_recon_edge_mat, high_confidence_edges, all_transitions_sig_hits)

    # move next five lines into a function
    results_object$transition_edge_matrix <- trans_mat_results$trans_edge_mat
    results_object$phenotype_transition_edge_matrix <- trans_mat_results$p_trans_mat
    results_object$hit_pvals <- corrected_pvals_all_transitions$hit_pvals
    results_object$sig_hits <- all_transitions_sig_hits

    save_results_as_r_object(args$output_dir, args$output_name, results_object)

  } else { # discrete phenotype

    genotype_transition_edges <- rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)){
      genotype_transition_edges[[k]] <- geno_trans[[k]]$transition
    }
    branch_overlap_trans <- count_hits_on_edges(genotype_transition_edges,   pheno_trans$transition,       high_conf_ordered_by_edges, pheno_conf_ordered_by_edges)
    branch_overlap_recon <- count_hits_on_edges(geno_recon_ordered_by_edges, pheno_recon_ordered_by_edges, high_conf_ordered_by_edges, pheno_conf_ordered_by_edges)

    disc_trans_results <- calculate_hit_pvals_corrected(branch_overlap_trans, pheno_trans$transition,       args$tree, genotype, args$perm, args$alpha, high_conf_ordered_by_edges)
    disc_recon_results <- calculate_hit_pvals_corrected(branch_overlap_recon, pheno_recon_ordered_by_edges, args$tree, genotype, args$perm, args$alpha, high_conf_ordered_by_edges)

    hit_pvals_trans <- disc_trans_results$hit_pvals
    hit_pvals_recon <- disc_recon_results$hit_pvals

    corrected_pvals_trans <- get_sig_hits_while_correcting_for_multiple_testing(hit_pvals_trans, args$alpha)
    corrected_pvals_recon <- get_sig_hits_while_correcting_for_multiple_testing(hit_pvals_recon, args$alpha)

    results_object$hit_pvals_transition     <- corrected_pvals_trans$hit_pvals
    results_object$hit_pvals_reconstruction <- corrected_pvals_recon$hit_pvals
    results_object$sig_pvals_transition     <- corrected_pvals_trans$sig_pvals
    results_object$sig_pvals_reconstruction <- corrected_pvals_recon$sig_pvals

    discrete_plots(tr = args$tree, # add a test to check that p_recon_edges and g_recon_edges have Nedge(tree)
                   dir = args$output_dir,
                   name = args$output_name,
                   a = args$alpha,
                   annot = args$annot,
                   num_perm = args$perm,
                   recon_hit_vals = corrected_pvals_recon$hit_pvals,
                   trans_hit_vals = corrected_pvals_trans$hit_pvals,
                   p_recon_edges = pheno_recon_ordered_by_edges,
                   g_recon_edges = geno_recon_ordered_by_edges,
                   recon_perm_obs_results = disc_recon_results,
                   trans_perm_obs_results = disc_trans_results,
                   tr_and_pheno_hi_conf = high_confidence_edges,
                   geno_conf = high_conf_ordered_by_edges,
                   g_trans_edges = genotype_transition_edges,
                   p_trans_edges = pheno_trans$transition)

    save_results_as_r_object(args$output_dir, args$output_name, results_object)
  }
}
