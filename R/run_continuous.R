#' Run Continuous algorithm.
#'
#' @param args Object with all of the inputs necessary to run gwas.
#'
#' @return Saves two files: a pdf with plots of results and a .rda with log and
#'    all relevant results.
#' @noRd
run_continuous <- function(args){
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
  print("Finish ancestral reconstruction")
  print(Sys.time())

  # Include all transition edges (WT -> mutant and mutant -> WT) for both
  #  synchronous and continuous tests.
  geno_trans_synchronous <- AR$geno_trans
  if (args$group_genotype == TRUE & args$grouping_method == "post-ar") {
    geno_trans_phyc <-
      prep_geno_trans_for_phyc(genotype, geno_trans_synchronous)
    grouped_geno <- group_genotypes_post_ar(args$tree,
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
    geno_conf_ordered_by_edges <- geno_recon_ordered_by_edges <-
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

  # CALCULATE CONVERGENCE -----------------------------------------------------#
  convergence <-
    calculate_continuous_convergence(AR$pheno_recon_and_conf$recon_edge_mat,
                                      hi_conf)

  # RUN PERMUTATION TEST ------------------------------------------------------#
  results_all_transitions <-
    calc_sig(hi_conf,
             args$perm,
             args$tree,
             AR$pheno_recon_and_conf$recon_edge_mat)
  print("Finish permutation test:")
  print(Sys.time())

  # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
  corrected_pvals_all_trans <-
    get_sig_hit_and_mult_test_corr(results_all_transitions$pvals, args$fdr)

  corrected_pvals_all_trans$hit_pvals <- -log(corrected_pvals_all_trans$hit_pvals)
  corrected_pvals_all_trans$sig_pvals <- -log(corrected_pvals_all_trans$sig_pvals)

  # SAVE AND PLOT RESULTS -----------------------------------------------------#
  phenotype_vector <-
    prepare_phenotype(args$phenotype, args$discrete_or_continuous, args$tree)
  trans_mat_results <-
    plot_continuous_results(disc_cont = "continuous",
                            tr = args$tree,
                            fdr = -log(args$fdr),
                            dir = args$output_dir,
                            name = args$output_name,
                            pval_all_transition = corrected_pvals_all_trans,
                            pheno_vector = phenotype_vector,
                            perm = args$perm,
                            results_all_trans = results_all_transitions,
                            pheno_anc_rec = AR$pheno_recon_and_conf$node_anc_rec,
                            geno_reconstruction = hi_conf$geno_recon_edge,
                            geno_confidence = hi_conf$high_conf_ordered_by_edges,
                            geno_transition = hi_conf$genotype_transition,
                            geno_mat = hi_conf$genotype,
                            pheno_recon_ordered_by_edges = AR$pheno_recon_and_conf$recon_edge_mat,
                            tr_and_pheno_hi_conf = hi_conf$high_conf_ordered_by_edges,
                            all_trans_sig_hits = corrected_pvals_all_trans$hit_pvals,
                            group_logical = args$group_genotype,
                            snp_in_gene = geno$snps_per_gene)
  results_object$hi_confidence_transition_edge <-
    hi_conf$hi_confidence_transition_edge
  results_object$num_hi_conf_transition_edge <-
    hi_conf$num_hi_conf_transition_edge
  results_object$dropped_genotypes <-
    hi_conf$dropped_genotypes
  results_object$genotype_transition_edge <-
    trans_mat_results$trans_dir_edge_mat
  results_object$phenotype_transition_edges <-
    trans_mat_results$p_trans_mat
  results_object$delta_pheno_table <-
    trans_mat_results$delta_pheno_table
  results_object$delta_pheno_list <- trans_mat_results$delta_pheno_list
  results_object$hit_pvals <- corrected_pvals_all_trans$hit_pvals
  results_object$sig_hits <- corrected_pvals_all_trans$sig_pvals
  results_object$convergence <- convergence
  names(results_object$convergence)[1] <- "N"

  results_object$phylogenetic_signal <- unname(calculate_lambda(args$phenotype,
                                                                args$tree))

  results_object$raw_pvals <-
    as.data.frame(as.matrix(results_all_transitions$pvals),
                  stringsAsFactors = FALSE)
  colnames(results_object$raw_pvals) <- "neg_log_unadjusted_pvals"
  results_object$raw_pvals$neg_log_unadjusted_pvals <-
    as.numeric(results_object$raw_pvals$neg_log_unadjusted_pvals)
  results_object$raw_pvals$neg_log_unadjusted_pvals <-
    -log(results_object$raw_pvals$neg_log_unadjusted_pvals)

  save_results_as_r_object(args$output_dir,
                           args$output_name,
                           results_object,
                           "continuous",
                           args$group_genotype)
}
