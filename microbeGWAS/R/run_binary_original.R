# TODO: go through script line by and line and write unit tests for any untested functions.
# TODO clean up script flow to make process more clear. Replace if / else with functions? How else to improve it?
run_binary_original <- function(args){
  # FORMAT INPUTS -------------------------------------------------------------#
  results_object <- NULL
  results_object$log <- capture.output(sessionInfo()) # log session info

  geno <- prepare_genotype(args$group_genotype, args$genotype, args$tree, args$gene_snp_lookup)
  genotype <- geno$genotype
  results_object$convergence_not_possible_genotypes <- geno$convergence_not_possible_genotypes

  AR <- prepare_ancestral_reconstructions(args$tree, args$phenotype, genotype, args$discrete_or_continuous)

  geno_trans_concomitant <- AR$geno_trans # Include all transition edges (WT -> mutant and mutant -> WT). For discrete concomitant and continuous tests.
  geno_trans_original    <- prepare_genotype_transitions_for_original_discrete_test(genotype, geno_trans_concomitant) # Keep only WT -> mutant transitions.

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

  hi_conf_original <- prepare_high_confidence_objects(geno_trans_original, args$tree, AR$pheno_recon_and_conf$tip_and_node_rec_conf, args$bootstrap_cutoff, genotype, geno_conf_ordered_by_edges, geno_recon_ordered_by_edges, geno$snps_per_gene)

  genotype_transition_edges <- rep(list(0), ncol(hi_conf_original$genotype))
  for (k in 1:ncol(hi_conf_original$genotype)){
    genotype_transition_edges[[k]] <- hi_conf_original$genotype_transition[[k]]$transition
  }

  pheno_trans <- identify_transition_edges(args$tree, args$phenotype, 1, AR$pheno_recon_and_conf$node_anc_rec, args$discrete_or_continuous)
  pheno_recon_ordered_by_edges <- reorder_tips_and_nodes_to_edges(AR$pheno_recon_and_conf$tip_and_node_recon, args$tree)

  # RUN PERMUTATION TEST ------------------------------------------------------#
  disc_recon_results <- discrete_calculate_pvals(genotype_transition_edges, pheno_recon_ordered_by_edges, args$tree, hi_conf_original$genotype, args$perm, args$fdr, hi_conf_original$high_conf_ordered_by_edges)

  # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
  corrected_pvals_recon <- get_sig_hits_while_correcting_for_multiple_testing(disc_recon_results$hit_pvals, args$fdr)


  # SAVE AND PLOT RESULTS -----------------------------------------------------#
  #discrete_plots(tr = args$tree, # add a test to check that p_recon_edges and g_recon_edges have Nedge(tree)
  #               dir = args$output_dir,
  #               name = args$output_name,
  #               fdr = args$fdr,
  #               annot = args$annot,
  #               num_perm = args$perm,
  #               recon_hit_vals = corrected_pvals_recon$hit_pvals,
  #               trans_hit_vals = corrected_pvals_trans$hit_pvals,
  #               p_recon_edges = pheno_recon_ordered_by_edges,
                 # g_recon_edges = geno_recon_ordered_by_edges,
  #               g_recon_edges = genotype_transition_edges,
  #               recon_perm_obs_results = disc_recon_results,
  #               trans_perm_obs_results = disc_trans_results,
  #               tr_and_pheno_hi_conf = hi_conf_original$tr_and_pheno_hi_conf,
  #               geno_conf = hi_conf_original$high_conf_ordered_by_edges,
  #               g_trans_edges = genotype_transition_edges,
  #               p_trans_edges = pheno_trans$transition,
  #               snp_in_gene = geno$snps_per_gene)

  results_object$contingency_table_recon <- create_contingency_table(genotype_transition_edges, pheno_recon_ordered_by_edges, hi_conf_original$genotype)
  results_object$hit_pvals_reconstruction <- corrected_pvals_recon$hit_pvals
  results_object$sig_pvals_reconstruction <- corrected_pvals_recon$sig_pvals
  results_object$original_high_confidence_trasition_edges     <- hi_conf_original$high_confidence_trasition_edges
  results_object$original_num_high_confidence_trasition_edges <- hi_conf_original$num_high_confidence_trasition_edges
  results_object$original_dropped_genotypes                   <- hi_conf_original$dropped_genotypes

  save_results_as_r_object(args$output_dir, args$output_name, results_object)
} # end run_binary_original()
