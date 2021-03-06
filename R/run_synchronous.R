#' Run Synchronous algorithm.
#'
#' @param args Object with all of the inputs necessary to run gwas.
#'
#' @return Saves two files: a pdf with plots of results and a .rda with log and
#'   all relevant results.
#' @noRd
run_synchronous <- function(args){
  # FORMAT INPUTS -------------------------------------------------------------#
  results_object <- NULL
  results_object$log <- utils::capture.output(utils::sessionInfo())

  args$tree <- format_tree(args$tree)
  args$phenotype <- match_order_to_tree_tips(args$tree, args$phenotype)
  args$genotype <- match_order_to_tree_tips(args$tree, args$genotype)

  geno <- prepare_genotype(args$group_genotype,
                           args$genotype,
                           args$tree,
                           args$gene_snp_lookup,
                           args$grouping_method)
  genotype <- geno$genotype
  results_object$no_convergence_genotypes <- geno$no_convergence_genotypes

  AR <- prepare_ancestral_reconstructions(args$tree,
                                          args$phenotype,
                                          genotype,
                                          args$discrete_or_continuous)
  print(paste0("Ancestral reconstruction complete: ", Sys.time()))

  # Include all transition edges (WT -> mutant and mutant -> WT). For
  #  synchronous and continuous tests.
  geno_trans_synchronous <- AR$geno_trans

  if (args$group_genotype == TRUE & args$grouping_method == "post-ar") {
    geno_trans_phyc <- prep_geno_trans_for_phyc(genotype, AR$geno_trans)
    grouped_geno <- group_genotypes_post_ar(args$tree,
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
    results_object$no_convergence_genotypes <-
      grouped_geno$no_convergence_genotypes
  } else {
    geno_conf_ordered_by_edges <-
      geno_recon_ordered_by_edges <-
      rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)) {
      geno_conf_ordered_by_edges[[k]] <-
        reorder_tip_and_node_to_edge(
          AR$geno_recon_and_conf[[k]]$tip_and_node_rec_conf, args$tree)
      geno_recon_ordered_by_edges[[k]] <-
        reorder_tip_and_node_to_edge(
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

  # CALCULATE CONVERGENCE -----------------------------------------------------#
  convergence <- calculate_synchronous_convergence(genotype_transition_edges,
                                                   pheno_trans,
                                                   hi_conf)

  # RUN PERMUTATION TEST ------------------------------------------------------#
  disc_trans_results <-
    discrete_calculate_pvals(genotype_transition_edges,
                             pheno_trans$transition,
                             args$tree,
                             hi_conf$genotype,
                             args$perm,
                             args$fdr,
                             hi_conf$high_conf_ordered_by_edges,
                             convergence)
  print(paste0("Permutation test complete: ", Sys.time()))

  # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
  corrected_pvals_trans <-
    get_sig_hit_and_mult_test_corr(disc_trans_results$hit_pvals, args$fdr)
  corrected_pvals_trans$hit_pvals <- -log(corrected_pvals_trans$hit_pvals)
  corrected_pvals_trans$sig_pvals <- -log(corrected_pvals_trans$sig_pvals)

  # SAVE AND PLOT RESULTS -----------------------------------------------------#
  plot_synchronous_results(
    tr = args$tree,
    dir = args$output_dir,
    name = args$output_name,
    fdr = -log(args$fdr),
    num_perm = args$perm,
    trans_hit_vals = corrected_pvals_trans$hit_pvals,
    trans_perm_obs_results = disc_trans_results,
    tr_and_pheno_hi_conf = hi_conf$tr_and_pheno_hi_conf,
    geno_confidence = hi_conf$high_conf_ordered_by_edges,
    g_trans_edges = genotype_transition_edges,
    p_trans_edges = pheno_trans$transition,
    snp_in_gene = geno$snps_per_gene,
    prefix = "synchronous",
    grouped_logical = args$group_genotype,
    tr_type = args$tree_type,
    strain_key = args$strain_key)

  results_object$contingency_table <-
    create_contingency_table(genotype_transition_edges,
                             pheno_trans$transition,
                             hi_conf$genotype,
                             "synchronous")
  results_object$hit_pvals <- corrected_pvals_trans$hit_pvals
  results_object$sig_pvals <- corrected_pvals_trans$sig_pvals
  results_object$hi_confidence_transition_edge <-
    hi_conf$hi_confidence_transition_edge
  results_object$num_hi_conf_transition_edge <-
    hi_conf$num_high_conf_trans_edges
  results_object$dropped_genotypes <- hi_conf$dropped_genotypes
  results_object$convergence <- convergence
  names(results_object$convergence)[1] <- "N"
  results_object$phylogenetic_signal <- unname(calculate_d(args$phenotype,
                                                           args$tree))

  results_object$raw_pvals <-
    as.data.frame(as.matrix(disc_trans_results$hit_pvals),
                  stringsAsFactors = FALSE)
  colnames(results_object$raw_pvals) <- "neg_log_unadjusted_pvals"
  results_object$raw_pvals$neg_log_unadjusted_pvals <-
    as.numeric(results_object$raw_pvals$neg_log_unadjusted_pvals)
  results_object$raw_pvals$neg_log_unadjusted_pvals <-
    -log(results_object$raw_pvals$neg_log_unadjusted_pvals)


  save_results_as_r_object(args$output_dir,
                           args$output_name,
                           results_object,
                           "synchronous",
                           args$group_genotype)
}
