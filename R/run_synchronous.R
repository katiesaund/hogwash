#' run_synchronous
#'
#' @description Run Synchronous algorithm.
#'
#' @param args Object with all of the inputs necessary to run gwas.
#'
#' @return Saves two files: a pdf with plots of results and a .rda with log and
#'    all relevant results.
#' @noRd
run_synchronous <- function(args){
  # FORMAT INPUTS -------------------------------------------------------------#
  results_object <- NULL
  results_object$log <- utils::capture.output(utils::sessionInfo())
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
  # Include all transition edges (WT -> mutant and mutant -> WT). For
  #  synchronous and continuous tests.
  geno_trans_synchronous <- AR$geno_trans

  # Keep only WT -> mutant transitions.
  geno_trans_phyc <- prep_geno_trans_for_phyc(genotype, AR$geno_trans)

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
    geno_conf_ordered_by_edges  <- grouped_geno$geno_conf_ordered_by_edges
    geno_trans_synchronous <- grouped_geno$geno_trans_synchronous
    geno_trans_phyc <- grouped_geno$geno_trans_phyc
    results_object$convergence_not_possible_genotypes <-
      grouped_geno$convergence_not_possible_genotypes
  } else {
    geno_conf_ordered_by_edges <-
      geno_recon_ordered_by_edges <-
      rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)) {
      geno_conf_ordered_by_edges[[k]] <-
        reorder_tips_and_nodes_to_edges(
          AR$geno_recon_and_conf[[k]]$tip_and_node_rec_conf, args$tree)
      geno_recon_ordered_by_edges[[k]] <-
        reorder_tips_and_nodes_to_edges(
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

  genotype_transition_edges <- rep(list(0), ncol(hi_conf$genotype))
  for (k in 1:ncol(hi_conf$genotype)) {
    genotype_transition_edges[[k]] <-
      hi_conf$genotype_transition[[k]]$transition
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
                             hi_conf$genotype,
                             args$perm,
                             args$fdr,
                             hi_conf$high_conf_ordered_by_edges)

  # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
  corrected_pvals_trans <-
    get_sig_hit_and_mult_test_corr(disc_trans_results$hit_pvals, args$fdr)

  # SAVE AND PLOT RESULTS -----------------------------------------------------#
  discrete_plot_trans(tr = args$tree, dir = args$output_dir,
                     name = args$output_name,
                     fdr = args$fdr,
                     num_perm = args$perm,
                     trans_hit_vals = corrected_pvals_trans$hit_pvals,
                     trans_perm_obs_results = disc_trans_results,
                     tr_and_pheno_hi_conf = hi_conf$tr_and_pheno_hi_conf,
                     geno_confidence = hi_conf$high_conf_ordered_by_edges,
                     g_trans_edges = genotype_transition_edges,
                     p_trans_edges = pheno_trans$transition,
                     snp_in_gene = geno$snps_per_gene,
                     prefix = "synchronous",
                     grouped_logical = args$group_genotype)

  results_object$contingency_table <-
    create_contingency_table(genotype_transition_edges,
                             pheno_trans$transition,
                             hi_conf$genotype,
                             "synchronous")
  results_object$hit_pvals <- corrected_pvals_trans$hit_pvals
  results_object$sig_pvals <- corrected_pvals_trans$sig_pvals
  results_object$high_confidence_transition_edges <-
    hi_conf$high_confidence_transition_edges
  results_object$num_high_confidence_transition_edges <-
    hi_conf$num_high_conf_trans_edges
  results_object$dropped_genotypes <- hi_conf$dropped_genotypes
  save_results_as_r_object(args$output_dir,
                           args$output_name,
                           results_object,
                           "synchronous",
                           args$group_genotype)
} # end run_synchronous()
