run_binary_transition <- function(args){
  # FORMAT INPUTS -------------------------------------------------------------#
  results_object <- NULL
  results_object$log <- capture.output(sessionInfo()) # log session info
  args$tree <- format_tree(args$tree)


  geno <- prepare_genotype(args$group_genotype,
                           args$genotype,
                           args$tree,
                           args$gene_snp_lookup)
  genotype <- geno$genotype
  results_object$convergence_not_possible_genotypes <-
    geno$convergence_not_possible_genotypes

  AR <- prepare_ancestral_reconstructions(args$tree,
                                          args$phenotype,
                                          genotype,
                                          args$discrete_or_continuous)

  geno_trans_concomitant <- AR$geno_trans # Include all transition edges (WT -> mutant and mutant -> WT). For discrete concomitant and continuous tests.
  geno_trans_original <-
    prepare_genotype_transitions_for_original_discrete_test(genotype,
                                                            AR$geno_trans) # Keep only WT -> mutant transitions.

  if (args$group_genotype) {
    grouped_geno <- group_genotypes(args$tree,
                                    genotype,
                                    AR$geno_recon_and_conf,
                                    geno_trans_concomitant,
                                    geno_trans_original,
                                    geno$gene_snp_lookup,
                                    geno$unique_genes)
    genotype                    <- grouped_geno$genotype
    geno_recon_ordered_by_edges <- grouped_geno$geno_recon_ordered_by_edges
    geno_conf_ordered_by_edges  <- grouped_geno$geno_conf_ordered_by_edges
    geno_trans_concomitant      <- grouped_geno$geno_trans_concomitant
    geno_trans_original         <- grouped_geno$geno_trans_original
    results_object$convergence_not_possible_genotypes <-
      grouped_geno$convergence_not_possible_genotypes
  } else {
    geno_conf_ordered_by_edges <-
      geno_recon_ordered_by_edges <-
      rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)) {
      geno_conf_ordered_by_edges[[k]] <-
        reorder_tips_and_nodes_to_edges(AR$geno_recon_and_conf[[k]]$tip_and_node_rec_conf,
                                        args$tree)
      geno_recon_ordered_by_edges[[k]] <-
        reorder_tips_and_nodes_to_edges(AR$geno_recon_and_conf[[k]]$tip_and_node_recon,
                                        args$tree)
    }
  }
  hi_conf_concomitant <-
    prepare_high_confidence_objects(geno_trans_concomitant,
                                    args$tree,
                                    AR$pheno_recon_and_conf$tip_and_node_rec_conf,
                                    args$bootstrap_cutoff,
                                    genotype,
                                    geno_conf_ordered_by_edges,
                                    geno_recon_ordered_by_edges,
                                    geno$snps_per_gene)
  genotype_transition_edges <- rep(list(0), ncol(hi_conf_concomitant$genotype))
  for (k in 1:ncol(hi_conf_concomitant$genotype)) {
    genotype_transition_edges[[k]] <-
      hi_conf_concomitant$genotype_transition[[k]]$transition
  }

  pheno_trans <- identify_transition_edges(args$tree,
                                           args$phenotype,
                                           1,
                                           AR$pheno_recon_and_conf$node_anc_rec,
                                           args$discrete_or_continuous)

  # RUN PERMUTATION TEST ------------------------------------------------------#
  disc_trans_results <-
    discrete_calculate_pvals(genotype_transition_edges,
                             pheno_trans$transition,
                             args$tree,
                             hi_conf_concomitant$genotype,
                             args$perm,
                             args$fdr,
                             hi_conf_concomitant$high_conf_ordered_by_edges)

  # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
  corrected_pvals_trans <-
    get_sig_hits_while_correcting_for_multiple_testing(disc_trans_results$hit_pvals,
                                                       args$fdr)

  # SAVE AND PLOT RESULTS -----------------------------------------------------#
  discrete_plot_trans(tr = args$tree, dir = args$output_dir,
                     name = args$output_name, fdr = args$fdr,
                     annot = args$annot, num_perm = args$perm,
                     trans_hit_vals = corrected_pvals_trans$hit_pvals,
                     trans_perm_obs_results = disc_trans_results,
                     tr_and_pheno_hi_conf = hi_conf_concomitant$tr_and_pheno_hi_conf,
                     geno_confidence = hi_conf_concomitant$high_conf_ordered_by_edges,
                     g_trans_edges = genotype_transition_edges,
                     p_trans_edges = pheno_trans$transition,
                     snp_in_gene = geno$snps_per_gene)

  results_object$contingency_table_trans <-
    create_contingency_table(genotype_transition_edges,
                             pheno_trans$transition,
                             hi_conf_concomitant$genotype)
  results_object$hit_pvals_transition <- corrected_pvals_trans$hit_pvals
  results_object$sig_pvals_transition <- corrected_pvals_trans$sig_pvals
  results_object$concomitant_high_confidence_trasition_edges <-
    hi_conf_concomitant$high_confidence_trasition_edges
  results_object$concomitant_num_high_confidence_trasition_edges <-
    hi_conf_concomitant$num_high_confidence_trasition_edges
  results_object$concomitant_dropped_genotypes <-
    hi_conf_concomitant$dropped_genotypes
  save_results_as_r_object(args$output_dir, args$output_name, results_object)
} # end run_binary_transition()
