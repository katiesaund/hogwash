run_phyc <- function(args){
  # READ IN ARGUMENTS -----------------------------------------------------------#


  # FORMAT INPUTS ---------------------------------------------------------------#
  simplified_genotype <- reduce_redundancy(args$genotype, args$tree, args$output_dir, args$output_name) # Remove genotypes that are too rare or too commmon for (1) convergence to be possible and (2) for ancestral reconstruction to work
  genotype <- simplified_genotype$mat
  convergence_not_possible <- simplified_genotype$dropped_genotype_names
  results_object <- NULL
  results_object$convergence_not_possible_genotypes <- convergence_not_possible
  phenotype_vector <- convert_matrix_to_vector(args$phenotype)

  session_log <- capture.output(sessionInfo()) # log session info
  results_object$log <- session_log

  # -----------------------------------------------------------------------------#
  # PHYC
  # -----------------------------------------------------------------------------#

  # ANCESTRAL RECONSTRUCTION OF PHENOTYPE ---------------------------------------#
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

  # IDENTIFY TRANSITION EDGES AND REFORMAT
  for (k in 1:ncol(genotype)){
    geno_trans[[k]] <- identify_transition_edges(args$tree, genotype, k, geno_recon_and_conf[[k]]$node_anc_rec, "discrete")
    geno_conf_ordered_by_edges[[k]]  <- reorder_tips_and_nodes_to_edges(geno_recon_and_conf[[k]]$tip_and_node_rec_conf, args$tree)
    geno_recon_ordered_by_edges[[k]] <- reorder_tips_and_nodes_to_edges(geno_recon_and_conf[[k]]$tip_and_node_recon, args$tree)
  }

  # IDENTIFY HIGH CONFIDENCE EDGES (BOOTSTRAP, PHENOTYPE RECON, LENGTH, GENOTYPE RECONSTRUCTION)
  # TREE BOOTSTRAP, PHENOTYPUE RECONSTRUCTION CONFIDENCE, AND EDGE LENGTHS
  high_confidence_edges <- pheno_conf_ordered_by_edges + tree_conf_ordered_by_edges + short_edges == 3
  all_high_confidence_edges <- rep(list(0), ncol(genotype))

  # ADD IN GENOTYPE RECONSTRUCTION CONFIDENCE
  for (k in 1:ncol(genotype)){
    all_high_confidence_edges[[k]] <- as.numeric(geno_conf_ordered_by_edges[[k]] + high_confidence_edges == 2)
  }

  # SAVE FILE WITH NUMBER OF HIGH CONFIDENCE TRANSITION EDGES PER GENOTYPE-------#

  results_object$high_confidence_trasition_edges <- high_confidence_edges
  #write.table(high_confidence_edges, commented out on 2019-02-25
  #            file = paste0(args$output_dir, "/phyc_", args$output_name, "_tree_pheno_length_hi_conf_edges.tsv"),
  #            sep = "\t",
  #            quote = FALSE,
  #            row.names = TRUE,
  #            col.names = TRUE)

  num_high_confidence_trasition_edges <- report_num_high_confidence_trans_edge(geno_trans, all_high_confidence_edges, colnames(genotype), args$output_dir, args$output_name)
  results_object$num_high_confidence_trasition_edges <- num_high_confidence_trasition_edges

  # KEEP ONLY GENOTYPES WITH AT LEAST TWO HIGH CONFIDENCE TRANSITION EDGES ------#
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
    # 2019-02-27 changed phenotype input from pheno_recon_edge_mat to pheno_trans, which it should be.
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

    #pheno_plus_recon_by_edges <- c(phenotype_vector, pheno_recon_ordered_by_edges)
    # ^ this was wrong because it's adding the phenotype vector twice to the ancestral reconstruction, which aparently already had it.

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
                   # pheno_anc_recon  = pheno_recon_and_conf$node_anc_rec,
                   recon_perm_obs_results = disc_recon_results,
                   trans_perm_obs_results = disc_trans_results,
                   tr_and_pheno_hi_conf = high_confidence_edges,
                   geno_conf = high_conf_ordered_by_edges,
                   g_trans_edges = genotype_transition_edges,
                   p_trans_edges = pheno_trans$transition)

    # htmp_tr <- create_heatmap_compatible_tree(args$tree)
    # if(!is.null(args$annotation)){
    #   simple_annotation <- args$annotation[ , 1, drop = FALSE]
    # }else{
    #   simple_annotation <- NULL
    # }
    #
    # if (nrow(corrected_pvals_recon$sig_pvals) > 0){
    #   recon_heatmap_filename_start <- paste(args$output_dir, "/phyc_", args$output_name, "_reconstruction_", sep = "")
    #   recon_heatmap_title    <- paste(args$output_name, " reconstruction sig hits", sep = "")
    #   plot_sig_hits_summary(htmp_tr, args$tree, genotype, args$phenotype, simple_annotation, corrected_pvals_recon$sig_pvals, recon_heatmap_title, recon_heatmap_filename_start, high_confidence_edges, "recon", short_edges, pheno_conf_ordered_by_edges, tree_conf_ordered_by_edges, pheno_recon_ordered_by_edges,  high_conf_ordered_by_edges, geno_recon_ordered_by_edges, args$output_name)
    # }
    #
    # if (nrow(corrected_pvals_trans$sig_pvals) > 0){
    #   trans_heatmap_filename_start <- paste(args$output_dir, "/phyc_", args$output_name, "_transition_", sep = "")
    #   trans_heatmap_title    <- paste(args$output_name, " transition sig hits", sep = "")
    #   plot_sig_hits_summary(htmp_tr, args$tree, genotype, args$phenotype, simple_annotation, corrected_pvals_trans$sig_pvals, trans_heatmap_title, trans_heatmap_filename_start, high_confidence_edges, "trans", short_edges, pheno_conf_ordered_by_edges, tree_conf_ordered_by_edges, pheno_trans$transition,  high_conf_ordered_by_edges, genotype_transition_edges, args$output_name)
    # }

    save_results_as_r_object(args$output_dir, args$output_name, results_object)
  }
}
