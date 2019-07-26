#' run_continuous
#'
#' @description Run Continuous algorithm.
#'
#' @param args Object with all of the inputs necessary to run gwas.
#'
#' @return Saves two files: a pdf with plots of results and a .rda with log and
#'    all relevant results.
#' @noRd
#'
run_continuous <- function(args){
  # FORMAT INPUTS -------------------------------------------------------------#
  results_object <- NULL
  results_object$log <- utils::capture.output(utils::sessionInfo())
  args$tree <- format_tree(args$tree)

  geno <- prepare_genotype(args$group_genotype,
                           args$genotype,
                           args$tree,
                           args$gene_snp_lookup)
  genotype <- geno$genotype
  results_object$no_convergence_genotypes <- geno$no_convergence_genotypes
  AR <- prepare_ancestral_reconstructions(args$tree,
                                          args$phenotype,
                                          genotype,
                                          args$discrete_or_continuous)
  # Include all transition edges (WT -> mutant and mutant -> WT). For
  #  synchronous and continuous tests.
  geno_trans_synchronous <- AR$geno_trans

  # Keep only WT -> mutant transitions.
  geno_trans_phyc <- prep_geno_trans_for_phyc(genotype, geno_trans_synchronous)
  if (args$group_genotype) {
    grouped_geno <- group_genotypes(args$tree,
                                    genotype,
                                    AR$geno_recon_and_conf,
                                    geno_trans_synchronous,
                                    geno_trans_phyc,
                                    geno$gene_snp_lookup,
                                    geno$unique_genes)
    genotype <- grouped_geno$genotype
    geno_recon_ordered_by_edges <- grouped_geno$geno_recon_ordered_by_edges
    geno_conf_ordered_by_edges <- grouped_geno$geno_conf_ordered_by_edges
    geno_trans_synchronous <- grouped_geno$geno_trans_synchronous
    results_object$no_convergence_genotypes <-
      grouped_geno$no_convergence_genotypes
  } else {
    geno_conf_ordered_by_edges <-
      geno_recon_ordered_by_edges <-
      rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)) {
      geno_conf_ordered_by_edges[[k]] <- reorder_tip_and_node_to_edge(
        AR$geno_recon_and_conf[[k]]$tip_and_node_rec_conf, args$tree)
      geno_recon_ordered_by_edges[[k]] <- reorder_tip_and_node_to_edge(
        AR$geno_recon_and_conf[[k]]$tip_and_node_recon, args$tree)
    }
  }
  hi_conf <- prepare_high_confidence_objects(
    geno_trans_synchronous,
    args$tree,
    AR$pheno_recon_and_conf$tip_and_node_rec_conf,
    args$bootstrap_cutoff,
    genotype,
    geno_conf_ordered_by_edges,
    geno_recon_ordered_by_edges,
    geno$snps_per_gene)

  # RUN PERMUTATION TEST ------------------------------------------------------#
  results_all_transitions <-
    calc_sig(hi_conf$genotype,
             args$perm,
             hi_conf$genotype_transition,
             args$tree,
             AR$pheno_recon_and_conf$recon_edge_mat,
             hi_conf$high_conf_ordered_by_edges,
             hi_conf$geno_recon_edge)

  # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
  corrected_pvals_all_trans <-
    get_sig_hit_and_mult_test_corr(results_all_transitions$pvals, args$fdr)

  # SUBSET SIGNIFICANT HITS SO MEDIAN(DELTA PHENOTYPE) ON
  #  TRANSITION EDGES > MEDIAN(DELTA PHENOTYPE) ALL EDGES #
  all_transitions_sig_hits <-
    keep_hits_with_more_change_on_trans_edges(results_all_transitions,
                                              corrected_pvals_all_trans,
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
                          corrected_pvals_all_trans,
                          phenotype_vector,
                          args$perm,
                          results_all_transitions,
                          AR$pheno_recon_and_conf$node_anc_rec,
                          hi_conf$geno_recon_edge,
                          hi_conf$high_conf_ordered_by_edges,
                          hi_conf$genotype_transition,
                          hi_conf$genotype,
                          AR$pheno_recon_and_conf$recon_edge_mat,
                          hi_conf$high_conf_ordered_by_edges,
                          all_transitions_sig_hits,
                          args$group_genotype)
  results_object$high_confidence_trasition_edges <-
    hi_conf$high_confidence_trasition_edges
  results_object$num_high_confidence_trasition_edges <-
    hi_conf$num_high_confidence_trasition_edges
  results_object$dropped_genotypes <-
    hi_conf$dropped_genotypes
  results_object$genotype_transition_edge_matrix <-
    trans_mat_results$trans_dir_edge_mat
  results_object$phenotype_transition_edges <-
    trans_mat_results$p_trans_mat
  results_object$delta_pheno_table <-
    trans_mat_results$delta_pheno_table
  results_object$delta_pheno_list <- trans_mat_results$delta_pheno_list
  results_object$hit_pvals <- corrected_pvals_all_trans$hit_pvals
  results_object$sig_hits <- all_transitions_sig_hits
  save_results_as_r_object(args$output_dir,
                           args$output_name,
                           results_object,
                           "continuous",
                           args$group_genotype)
} # end run_continuous()
