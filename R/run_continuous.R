run_continuous <- function(args){
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
                                                            geno_trans_concomitant) # Keep only WT -> mutant transitions.

  if (args$group_genotype) {
    grouped_geno <- group_genotypes(args$tree,
                                    genotype,
                                    AR$geno_recon_and_conf,
                                    geno_trans_concomitant,
                                    geno_trans_original,
                                    geno$gene_snp_lookup,
                                    geno$unique_genes)
    genotype <- grouped_geno$genotype
    geno_recon_ordered_by_edges <- grouped_geno$geno_recon_ordered_by_edges
    geno_conf_ordered_by_edges <- grouped_geno$geno_conf_ordered_by_edges
    geno_trans_concomitant <- grouped_geno$geno_trans_concomitant
    results_object$convergence_not_possible_genotypes <-
      grouped_geno$convergence_not_possible_genotypes
  } else {
    geno_conf_ordered_by_edges <- geno_recon_ordered_by_edges <-
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

  # RUN PERMUTATION TEST ------------------------------------------------------#
  results_all_transitions <-
    calculate_genotype_significance(hi_conf_concomitant$genotype,
                                    args$perm,
                                    hi_conf_concomitant$genotype_transition,
                                    args$tree,
                                    AR$pheno_recon_and_conf$recon_edge_mat,
                                    hi_conf_concomitant$high_conf_ordered_by_edges,
                                    hi_conf_concomitant$geno_recon_edge)

  # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
  corrected_pvals_all_transitions <-
    get_sig_hits_while_correcting_for_multiple_testing(results_all_transitions$pvals,
                                                       args$fdr)

  # SUBSET SIGNIFICANT HITS SO MEDIAN(DELTA PHENOTYPE) ON TRANSITION EDGES > MEDIAN(DELTA PHENOTYPE) ALL EDGES #
  all_transitions_sig_hits <-
    keep_hits_with_more_change_on_trans_edges(results_all_transitions,
                                              corrected_pvals_all_transitions,
                                              args$fdr)

  # SAVE AND PLOT RESULTS -----------------------------------------------------#
  phenotype_vector <-
    prepare_phenotype(args$phenotype, args$discrete_or_continuous, args$tree)
  trans_mat_results <-
    plot_significant_hits("continuous",
                          args$tree,
                          args$fdr,
                          args$output_dir,
                          args$output_name,
                          corrected_pvals_all_transitions,
                          phenotype_vector,
                          args$annotation,
                          args$perm,
                          results_all_transitions,
                          AR$pheno_recon_and_conf$node_anc_rec,
                          hi_conf_concomitant$geno_recon_edge,
                          hi_conf_concomitant$high_conf_ordered_by_edges,
                          hi_conf_concomitant$genotype_transition,
                          hi_conf_concomitant$genotype,
                          AR$pheno_recon_and_conf$recon_edge_mat,
                          hi_conf_concomitant$high_conf_ordered_by_edges,
                          all_transitions_sig_hits)
  results_object$concomitant_high_confidence_trasition_edges <-
    hi_conf_concomitant$high_confidence_trasition_edges
  results_object$concomitant_num_high_confidence_trasition_edges <-
    hi_conf_concomitant$num_high_confidence_trasition_edges
  results_object$concomitant_dropped_genotypes <-
    hi_conf_concomitant$dropped_genotypes
  results_object$genotype_transition_edge_matrix <-
    trans_mat_results$trans_dir_edge_mat
  results_object$phenotype_transition_edge_matrix <-
    trans_mat_results$p_trans_mat
  results_object$delta_pheno_table <-
    trans_mat_results$delta_pheno_table
  results_object$delta_pheno_list <- trans_mat_results$delta_pheno_list
  results_object$hit_pvals <- corrected_pvals_all_transitions$hit_pvals
  results_object$sig_hits <- all_transitions_sig_hits
  save_results_as_r_object(args$output_dir, args$output_name, results_object)
} # end run_continuous()
