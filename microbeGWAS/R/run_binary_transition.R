# TODO: go through script line by and line and write unit tests for any untested functions.
# TODO clean up script flow to make process more clear. Replace if / else with functions? How else to improve it?
run_phyc <- function(args){
  # FORMAT INPUTS -------------------------------------------------------------#
  results_object <- NULL
  results_object$log <- capture.output(sessionInfo()) # log session info

  geno <- prepare_genotype(args$group_genotype, args$genotype, args$tree, args$gene_snp_lookup)
  genotype <- geno$genotype
  results_object$convergence_not_possible_genotypes <- geno$convergence_not_possible_genotypes

  AR <- prepare_ancestral_reconstructions(args$tree, args$phenotype, genotype, args$discrete_or_continuous)

  geno_trans_concomitant <- AR$geno_trans # Include all transition edges (WT -> mutant and mutant -> WT). For discrete concomitant and continuous tests.
  geno_trans_original    <- prepare_genotype_transitions_for_original_discrete_test(args$discrete_or_continuous, genotype, AR$geno_trans) # Keep only WT -> mutant transitions.

  if (args$group_genotype){
    grouped_geno                <- group_genotypes(args$tree, genotype, AR$geno_recon_and_conf, geno_trans_concomitant, geno_trans_original, geno$gene_snp_lookup, geno$unique_genes)

    genotype                    <- grouped_geno$genotype
    geno_recon_ordered_by_edges <- grouped_geno$geno_recon_ordered_by_edges
    geno_conf_ordered_by_edges  <- grouped_geno$geno_conf_ordered_by_edges
    geno_trans_concomitant      <- grouped_geno$geno_trans_concomitant
    geno_trans_original         <- grouped_geno$geno_trans_original
    results_object$convergence_not_possible_genotypes <- grouped_geno$convergence_not_possible_genotypes
  } else {
    geno_conf_ordered_by_edges <- geno_recon_ordered_by_edges <- rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)){
      geno_conf_ordered_by_edges[[k]]  <- reorder_tips_and_nodes_to_edges(AR$geno_recon_and_conf[[k]]$tip_and_node_rec_conf, args$tree)
      geno_recon_ordered_by_edges[[k]] <- reorder_tips_and_nodes_to_edges(AR$geno_recon_and_conf[[k]]$tip_and_node_recon,    args$tree)
    }
  }

  hi_conf_concomitant <- prepare_high_confidence_objects(geno_trans_concomitant, args$tree, AR$pheno_recon_and_conf$tip_and_node_rec_conf, args$bootstrap_cutoff, genotype, geno_conf_ordered_by_edges, geno_recon_ordered_by_edges, geno$snps_per_gene)
  results_object$concomitant_high_confidence_trasition_edges     <- hi_conf_concomitant$only_high_conf_geno_trans
  results_object$concomitant_num_high_confidence_trasition_edges <- hi_conf_concomitant$num_high_confidence_trasition_edges
  results_object$concomitant_dropped_genotypes <- hi_conf_concomitant$dropped_genotypes
  hi_conf_original <- NULL

  if (args$discrete_or_continuous == "discrete"){
    hi_conf_original <- prepare_high_confidence_objects(geno_trans_original,     args$tree, AR$pheno_recon_and_conf$tip_and_node_rec_conf, args$bootstrap_cutoff, genotype, geno_conf_ordered_by_edges, geno_recon_ordered_by_edges, geno$snps_per_gene)
    results_object$original_high_confidence_trasition_edges     <- hi_conf_original$only_high_conf_geno_trans
    results_object$original_num_high_confidence_trasition_edges <- hi_conf_original$num_high_confidence_trasition_edges
    results_object$original_dropped_genotypes                   <- hi_conf_original$dropped_genotypes
  }


  # # IDENTIFY HIGH CONFIDENCE EDGES (BOOTSTRAP, PHENOTYPE RECON, LENGTH, GENOTYPE RECONSTRUCTION)
  # # TREE BOOTSTRAP, PHENOTYPUE RECONSTRUCTION CONFIDENCE, AND EDGE LENGTHS
  # pheno_conf_ordered_by_edges <- reorder_tips_and_nodes_to_edges(AR$pheno_recon_and_conf$tip_and_node_rec_conf, args$tree)
  # tree_conf                   <- get_bootstrap_confidence(args$tree, args$bootstrap_cutoff)
  # tree_conf_ordered_by_edges  <- reorder_tips_and_nodes_to_edges(tree_conf, args$tree)
  # short_edges                 <- identify_short_edges(args$tree)
  #
  #
  # high_confidence_edges <- pheno_conf_ordered_by_edges + tree_conf_ordered_by_edges + short_edges == 3
  # all_high_confidence_edges <- rep(list(0), ncol(genotype))
  #
  # # ADD IN GENOTYPE RECONSTRUCTION CONFIDENCE
  # for (k in 1:ncol(genotype)){
  #   all_high_confidence_edges[[k]] <- as.numeric(geno_conf_ordered_by_edges[[k]] + high_confidence_edges == 2)
  # }
  #
  # # as of 2019-05-15 assign_high_confidence_to_transition_edges is so stringent that no genotype is getting included after this!
  # # TODO check that assign_high_confidence_to_transition_edges is now the not stringent version (ignores parent edges) - if so, remove this and the comment above.
  # # TODO to clean up code can I move all of the results_object$dummy_name <- dummy assignments to a function at the end so as to improve readability of code / remove extra lines
  # only_high_conf_geno_trans <- assign_high_confidence_to_transition_edges(args$tree, all_high_confidence_edges, AR$geno_trans, genotype)
  # results_object$high_confidence_trasition_edges <- only_high_conf_geno_trans
  # for (i in 1:ncol(genotype)){
  #   AR$geno_trans[[i]]$transition <- only_high_conf_geno_trans[[i]]
  # }
  # # how to plot:
  # # plot_tree_with_colored_edges(args$tree, geno_trans, all_high_confidence_edges, "grey", "red", "only new transitions", args$annot, "trans", 2)
  #
  # # SAVE FILE WITH NUMBER OF HIGH CONFIDENCE TRANSITION EDGES PER GENOTYPE-----#
  # # results_object$high_confidence_trasition_edges <- high_confidence_edges 2019-03-18 this is too simplistic-- updating using assign_high_confidence_to_transition_edges()
  # # TODO follow through on replacing high_confdience_edges as necessary
  # num_high_confidence_trasition_edges <- report_num_high_confidence_trans_edge(AR$geno_trans, all_high_confidence_edges, colnames(genotype))
  # results_object$num_high_confidence_trasition_edges <- num_high_confidence_trasition_edges
  #
  # # KEEP ONLY GENOTYPES WITH AT LEAST TWO HIGH CONFIDENCE TRANSITION EDGES ----#
  # geno_to_keep                  <- keep_at_least_two_high_conf_trans_edges(AR$geno_trans, all_high_confidence_edges)
  # geno_recon_ordered_by_edges   <- geno_recon_ordered_by_edges[geno_to_keep]
  # high_conf_ordered_by_edges    <- all_high_confidence_edges[geno_to_keep]
  # AR$geno_trans                 <- AR$geno_trans[geno_to_keep]
  #
  # dropped_genotypes <- get_dropped_genotypes(genotype, geno_to_keep)
  # results_object$dropped_genotypes <- dropped_genotypes
  #
  # genotype                      <- genotype[ , geno_to_keep, drop = FALSE]
  # geno$snps_per_gene <- geno$snps_per_gene[names(geno$snps_per_gene) %in% colnames(genotype)]

  # TODO break following if else into two seperate functions
  #convergence_algorithm_and_plotting <- function(args$tree, args$perm, args$fdr, args$phenotype, args$phenotype, args$output_dir, args$output_name, args$annotation, AR, hi_conf_concomitant, hi_conf_original, results_object){

  #  }


  # Do discrete first -- if it exists it's stuff can be added to the save results object step


  if (args$discrete_or_continuous == "continuous"){
    pheno_recon_edge_mat  <- AR$pheno_recon_and_conf$recon_edge_mat
    # RUN PERMUTATION TEST ------------------------------------------------------#
    results_all_transitions <- calculate_genotype_significance(genotype, args$perm, AR$geno_trans, args$tree, pheno_recon_edge_mat, high_conf_ordered_by_edges, geno_recon_ordered_by_edges)

    # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
    corrected_pvals_all_transitions <- get_sig_hits_while_correcting_for_multiple_testing(results_all_transitions$pvals, args$fdr)
    # SUBSET SIGNIFICANT HITS SO MEDIAN(DELTA PHENOTYPE) ON TRANSITION EDGES > MEDIAN(DELTA PHENOTYPE) ALL EDGES
    all_transitions_sig_hits <- keep_hits_with_more_change_on_trans_edges(results_all_transitions, corrected_pvals_all_transitions, args$fdr)
    # SAVE AND PLOT RESULTS -----------------------------------------------------#
    phenotype_vector <- prepare_phenotype(args$phenotype, args$phenotype, args$tree)
    trans_mat_results <- plot_significant_hits("continuous", args$tree, args$fdr, args$output_dir, args$output_name, corrected_pvals_all_transitions, phenotype_vector, args$annotation, args$perm, results_all_transitions, AR$pheno_recon_and_conf$node_anc_rec, geno_recon_ordered_by_edges, high_conf_ordered_by_edges, AR$geno_trans, genotype, pheno_recon_edge_mat, high_confidence_edges, all_transitions_sig_hits)
    # move next five lines into a function
    results_object$genotype_transition_edge_matrix <- trans_mat_results$trans_dir_edge_mat
    results_object$phenotype_transition_edge_matrix <- trans_mat_results$p_trans_mat
    results_object$delta_pheno_table <- trans_mat_results$delta_pheno_table
    results_object$delta_pheno_list <- trans_mat_results$delta_pheno_list
    results_object$hit_pvals <- corrected_pvals_all_transitions$hit_pvals
    results_object$sig_hits <- all_transitions_sig_hits
    save_results_as_r_object(args$output_dir, args$output_name, results_object)
  } else { # discrete phenotype
    genotype_transition_edges <- rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)){
      genotype_transition_edges[[k]] <- AR$geno_trans[[k]]$transition
    }

    pheno_trans           <- identify_transition_edges(args$tree, args$phenotype, 1, AR$pheno_recon_and_conf$node_anc_rec, args$discrete_or_continuous)
    pheno_recon_ordered_by_edges <- reorder_tips_and_nodes_to_edges(AR$pheno_recon_and_conf$tip_and_node_recon, args$tree)

    results_object$contingency_table_trans <- create_contingency_table(genotype_transition_edges, pheno_trans$transition,       genotype)
    results_object$contingency_table_recon <- create_contingency_table(genotype_transition_edges, pheno_recon_ordered_by_edges, genotype)

    # TODO make sure that the input into discrete_calculate_pvals is appropriate for overlap vs phyc tests.
    disc_trans_results <- discrete_calculate_pvals(genotype_transition_edges, pheno_trans$transition, args$tree, genotype, args$perm, args$fdr, high_conf_ordered_by_edges)
    disc_recon_results <- discrete_calculate_pvals(genotype_transition_edges, pheno_recon_ordered_by_edges, args$tree, genotype, args$perm, args$fdr, high_conf_ordered_by_edges)

    corrected_pvals_trans <- get_sig_hits_while_correcting_for_multiple_testing(disc_trans_results$hit_pvals, args$fdr)
    corrected_pvals_recon <- get_sig_hits_while_correcting_for_multiple_testing(disc_recon_results$hit_pvals, args$fdr)

    results_object$hit_pvals_transition     <- corrected_pvals_trans$hit_pvals
    results_object$hit_pvals_reconstruction <- corrected_pvals_recon$hit_pvals
    results_object$sig_pvals_transition     <- corrected_pvals_trans$sig_pvals
    results_object$sig_pvals_reconstruction <- corrected_pvals_recon$sig_pvals

    discrete_plots(tr = args$tree, # add a test to check that p_recon_edges and g_recon_edges have Nedge(tree)
                   dir = args$output_dir,
                   name = args$output_name,
                   fdr = args$fdr,
                   annot = args$annot,
                   num_perm = args$perm,
                   recon_hit_vals = corrected_pvals_recon$hit_pvals,
                   trans_hit_vals = corrected_pvals_trans$hit_pvals,
                   p_recon_edges = pheno_recon_ordered_by_edges,
                   # g_recon_edges = geno_recon_ordered_by_edges,
                   g_recon_edges = genotype_transition_edges,
                   recon_perm_obs_results = disc_recon_results,
                   trans_perm_obs_results = disc_trans_results,
                   tr_and_pheno_hi_conf = high_confidence_edges,
                   geno_conf = high_conf_ordered_by_edges,
                   g_trans_edges = genotype_transition_edges,
                   p_trans_edges = pheno_trans$transition,
                   snp_in_gene = geno$snps_per_gene)

    save_results_as_r_object(args$output_dir, args$output_name, results_object)
  }
}
