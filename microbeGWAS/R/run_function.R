run_phyc <- function(args){
  # FORMAT INPUTS ---------------------------------------------------------------#
  results_object <- NULL
  results_object$log <- capture.output(sessionInfo()) # log session info

  if (!args$built_from_snps) {
    simplified_genotype <- reduce_redundancy(args$genotype, args$tree) # Remove genotypes that are too rare or too commmon for (1) convergence to be possible and (2) for ancestral reconstruction to work
    genotype <- simplified_genotype$mat
    results_object$convergence_not_possible_genotypes <- simplified_genotype$dropped_genotype_names
    snps_per_gene <- NULL
  } else {
    genotypes_to_drop_because_not_present <- colnames(args$genotype)[colSums(args$genotype) == 0]
    genotype <- args$genotype[ , colSums(args$genotype) > 0] # we don't want to remove snps that are too rare or too common until snps are grouped into genes, then run on the grouped genes. But we need to remove SNPs that don't occur for ace to work.
    gene_snp_lookup <- args$gene_snp_lookup[!(args$gene_snp_lookup[ , 1] %in% genotypes_to_drop_because_not_present), , drop = FALSE]
    gene_snp_lookup <- gene_snp_lookup[gene_snp_lookup[ , 1] %in% colnames(genotype), , drop = FALSE]
    unique_genes <- unique(gene_snp_lookup[ , 2])
    snps_per_gene <- table(gene_snp_lookup[ , 2])
  }
  check_if_phenotype_normal(args$phenotype, args$discrete_or_continuous)
  check_if_convergence_occurs(args$phenotype, args$tree, args$discrete_or_continuous)
  phenotype_vector <- convert_matrix_to_vector(args$phenotype) # TODO add check that it's possible to have phenotype convergence
  check_convergence_possible(phenotype_vector, args$discrete_or_continuous)


  # ---------------------------------------------------------------------------#
  # PHYC
  # ---------------------------------------------------------------------------#

  # ANCESTRAL RECONSTRUCTION OF PHENOTYPE -------------------------------------#
  set.seed(1)
  pheno_recon_and_conf  <- ancestral_reconstruction_by_ML(args$tree, args$phenotype, 1, args$discrete_or_continuous)
  tree_conf             <- get_bootstrap_confidence(args$tree, args$bootstrap_cutoff)
  pheno_trans           <- identify_transition_edges(args$tree, args$phenotype, 1, pheno_recon_and_conf$node_anc_rec, args$discrete_or_continuous)
  # TODO identify high confdience pheno trans edges! 2019-03-28
  pheno_recon_edge_mat  <- pheno_recon_and_conf$recon_edge_mat
  short_edges           <- identify_short_edges(args$tree)
  pheno_conf_ordered_by_edges  <- reorder_tips_and_nodes_to_edges(pheno_recon_and_conf$tip_and_node_rec_conf, args$tree)
  pheno_recon_ordered_by_edges <- reorder_tips_and_nodes_to_edges(pheno_recon_and_conf$tip_and_node_recon, args$tree)
  tree_conf_ordered_by_edges   <- reorder_tips_and_nodes_to_edges(tree_conf, args$tree)

  # ANCESTRAL RECONSTRUCTION OF GENOTYPES ---------------------------------------#
  # INIATLIAZE DATA STRUCTS
  geno_recon_and_conf <- geno_trans <- rep(list(0), ncol(genotype))

  # PERFORM ANCESTRAL RECONSTRUCTION
  for(k in 1:ncol(genotype)){
    geno_recon_and_conf[[k]] <- ancestral_reconstruction_by_ML(args$tree, genotype, k, "discrete")
  }

  for (k in 1:ncol(genotype)){
    geno_trans[[k]] <- identify_transition_edges(args$tree, genotype, k, geno_recon_and_conf[[k]]$node_anc_rec, "discrete")
  }


  if (args$discrete_or_continuous == "discrete"){ #when we're doing original phyc
    # 2019-03-28 change geno_trans to only have WT -> mutant included for $transition to better reflect original phyC
    for (k in 1:ncol(genotype)){
      # update definition of $transition to be only WT -> mutant
      parent_WT_child_mutant <- 1 # 1 implies parent < child, -1 implies parent > child, 0 implies parent == child
      geno_trans[[k]]$transition <- as.numeric(geno_trans[[k]]$trans_dir == parent_WT_child_mutant)
    }
    # TODO what does this update break?
    # it breaks the discrete transition test, but should work well for the discrete original test.
  }

  if (args$built_from_snps){
    # TODO change this if statement into a function
    print("1")
    # CONVERT SNPS INTO GENES HERE
    # tip_and_node_ancestral_reconstruction
    temp_results <- build_gene_anc_recon_and_conf_from_snp(args$tree, genotype, geno_recon_and_conf, gene_snp_lookup)
    geno_recon_and_confidence_tip_node_recon <- temp_results$tip_node_recon
    geno_recon_and_confidence_tip_node_confidence <- temp_results$tip_node_conf

    # edge based transition
    geno_trans <- build_gene_trans_from_snp_trans(args$tree, genotype, geno_trans, gene_snp_lookup)

    # make new genotype (just at the tips, from the snps)
    genotype <- build_gene_genotype_from_snps(genotype, gene_snp_lookup)

    simplified_genotype <- reduce_redundancy(genotype, args$tree) # Remove genotypes that are too rare or too commmon for (1) convergence to be possible and (2) for ancestral reconstruction to work
    genotype <- simplified_genotype$mat
    results_object$convergence_not_possible_genotypes <- simplified_genotype$dropped_genotype_names
    genes_to_keep_in_consideration <- !(unique_genes %in% simplified_genotype$dropped_genotype_names)

    # remove redundancy from geno trans, geno_recon_and_confdience_tip_node_recon, and node_confidence
    geno_trans <- geno_trans[genes_to_keep_in_consideration]

    dummy <- geno_trans
    geno_trans <- rep(list(NULL), length(dummy))
    for (i in 1:length(dummy)){
      geno_trans[[i]]$transition <- as.numeric(as.character(unlist((dummy[i]))))
    }

    geno_recon_and_confidence_tip_node_recon <- geno_recon_and_confidence_tip_node_recon[genes_to_keep_in_consideration]
    geno_recon_and_confidence_tip_node_confidence <- geno_recon_and_confidence_tip_node_confidence[genes_to_keep_in_consideration]

    geno_conf_ordered_by_edges <- geno_recon_ordered_by_edges <- rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)){
      geno_recon_ordered_by_edges[[k]] <- reorder_tips_and_nodes_to_edges(geno_recon_and_confidence_tip_node_recon[[k]],      args$tree)
      geno_conf_ordered_by_edges[[k]]  <- reorder_tips_and_nodes_to_edges(geno_recon_and_confidence_tip_node_confidence[[k]], args$tree)
    }

    # TODO create a test to check/viz that I did the above assignments correctly and started from the correct piece of data.
  } else {
    geno_conf_ordered_by_edges <- geno_recon_ordered_by_edges <- rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)){
      geno_conf_ordered_by_edges[[k]]  <- reorder_tips_and_nodes_to_edges(geno_recon_and_conf[[k]]$tip_and_node_rec_conf, args$tree)
      geno_recon_ordered_by_edges[[k]] <- reorder_tips_and_nodes_to_edges(geno_recon_and_conf[[k]]$tip_and_node_recon,    args$tree)
    }
  }

  print("A")

  # IDENTIFY HIGH CONFIDENCE EDGES (BOOTSTRAP, PHENOTYPE RECON, LENGTH, GENOTYPE RECONSTRUCTION)
  # TREE BOOTSTRAP, PHENOTYPUE RECONSTRUCTION CONFIDENCE, AND EDGE LENGTHS
  high_confidence_edges <- pheno_conf_ordered_by_edges + tree_conf_ordered_by_edges + short_edges == 3
  all_high_confidence_edges <- rep(list(0), ncol(genotype))

  # ADD IN GENOTYPE RECONSTRUCTION CONFIDENCE
  for (k in 1:ncol(genotype)){
    all_high_confidence_edges[[k]] <- as.numeric(geno_conf_ordered_by_edges[[k]] + high_confidence_edges == 2)
  }

  print("B")
  # as of 2019-05-15 assign_high_confidence_to_transition_edges is so stringent that no genotype is getting included after this!
  only_high_conf_geno_trans <- assign_high_confidence_to_transition_edges(args$tree, all_high_confidence_edges, geno_trans, genotype)
  results_object$high_confidence_trasition_edges <- only_high_conf_geno_trans
  for (i in 1:ncol(genotype)){
    geno_trans[[i]]$transition <- only_high_conf_geno_trans[[i]]
  }
  # how to plot:
  # plot_tree_with_colored_edges(args$tree, geno_trans, all_high_confidence_edges, "grey", "red", "only new transitions", args$annot, "trans", 2)

  # SAVE FILE WITH NUMBER OF HIGH CONFIDENCE TRANSITION EDGES PER GENOTYPE-----#
  # results_object$high_confidence_trasition_edges <- high_confidence_edges 2019-03-18 this is too simplistic-- updating using assign_high_confidence_to_transition_edges()
  # TODO follow through on replacing high_confdience_edges as necessary
  num_high_confidence_trasition_edges <- report_num_high_confidence_trans_edge(geno_trans, all_high_confidence_edges, colnames(genotype))
  results_object$num_high_confidence_trasition_edges <- num_high_confidence_trasition_edges

  # KEEP ONLY GENOTYPES WITH AT LEAST TWO HIGH CONFIDENCE TRANSITION EDGES ----#
  geno_to_keep                  <- keep_at_least_two_high_conf_trans_edges(geno_trans, all_high_confidence_edges)
  geno_recon_ordered_by_edges   <- geno_recon_ordered_by_edges[geno_to_keep]
  high_conf_ordered_by_edges    <- all_high_confidence_edges[geno_to_keep]
  geno_trans                    <- geno_trans[geno_to_keep]

  dropped_genotypes <- get_dropped_genotypes(genotype, geno_to_keep)
  results_object$dropped_genotypes <- dropped_genotypes

  genotype                      <- genotype[ , geno_to_keep, drop = FALSE]
  snps_per_gene <- snps_per_gene[names(snps_per_gene) %in% colnames(genotype)]

  print("C")
  # break following if else into two seperate functions
  if (args$discrete_or_continuous == "continuous"){
    # RUN PERMUTATION TEST ------------------------------------------------------#
    results_all_transitions <- calculate_genotype_significance(genotype, args$perm, geno_trans, args$tree, pheno_recon_edge_mat, high_conf_ordered_by_edges, geno_recon_ordered_by_edges)

    print("D")
    # IDENTIFY SIGNIFICANT HITS USING FDR CORRECTION ----------------------------#
    corrected_pvals_all_transitions <- get_sig_hits_while_correcting_for_multiple_testing(results_all_transitions$pvals, args$fdr)
    print("E")
    # SUBSET SIGNIFICANT HITS SO MEDIAN(DELTA PHENOTYPE) ON TRANSITION EDGES > MEDIAN(DELTA PHENOTYPE) ALL EDGES
    all_transitions_sig_hits <- keep_hits_with_more_change_on_trans_edges(results_all_transitions, corrected_pvals_all_transitions, args$fdr)
    print("F")
    # SAVE AND PLOT RESULTS -----------------------------------------------------#
    trans_mat_results <- plot_significant_hits("continuous", args$tree, args$fdr, args$output_dir, args$output_name, corrected_pvals_all_transitions, phenotype_vector, args$annotation, args$perm, results_all_transitions, pheno_recon_and_conf$node_anc_rec, geno_recon_ordered_by_edges, high_conf_ordered_by_edges, geno_trans, genotype, pheno_recon_edge_mat, high_confidence_edges, all_transitions_sig_hits)
    print("G")
    # move next five lines into a function
    results_object$genotype_transition_edge_matrix <- trans_mat_results$trans_dir_edge_mat
    results_object$phenotype_transition_edge_matrix <- trans_mat_results$p_trans_mat
    results_object$delta_pheno_table <- trans_mat_results$delta_pheno_table
    results_object$delta_pheno_list <- trans_mat_results$delta_pheno_list
    results_object$hit_pvals <- corrected_pvals_all_transitions$hit_pvals
    results_object$sig_hits <- all_transitions_sig_hits
    print("H")
    save_results_as_r_object(args$output_dir, args$output_name, results_object)
    print("I")
  } else { # discrete phenotype
    print("D")
    genotype_transition_edges <- rep(list(0), ncol(genotype))
    for (k in 1:ncol(genotype)){
      genotype_transition_edges[[k]] <- geno_trans[[k]]$transition
    }

    results_object$contingency_table_trans <- create_contingency_table(genotype_transition_edges, pheno_trans$transition,       genotype)
    results_object$contingency_table_recon <- create_contingency_table(genotype_transition_edges, pheno_recon_ordered_by_edges, genotype)

    # TODO make sure that the input into calculate_hit_pvals_corrected is appropriate for overlap vs phyc tests.
    disc_trans_results <- calculate_hit_pvals_corrected(genotype_transition_edges, pheno_trans$transition, args$tree, genotype, args$perm, args$fdr, high_conf_ordered_by_edges)
    disc_recon_results <- calculate_hit_pvals_corrected(genotype_transition_edges, pheno_recon_ordered_by_edges, args$tree, genotype, args$perm, args$fdr, high_conf_ordered_by_edges)

    corrected_pvals_trans <- get_sig_hits_while_correcting_for_multiple_testing(disc_trans_results$hit_pvals, args$fdr)
    corrected_pvals_recon <- get_sig_hits_while_correcting_for_multiple_testing(disc_recon_results$hit_pvals, args$fdr)

    results_object$hit_pvals_transition     <- corrected_pvals_trans$hit_pvals
    results_object$hit_pvals_reconstruction <- corrected_pvals_recon$hit_pvals
    results_object$sig_pvals_transition     <- corrected_pvals_trans$sig_pvals
    results_object$sig_pvals_reconstruction <- corrected_pvals_recon$sig_pvals

    print("F")
    discrete_plots(tr = args$tree, # add a test to check that p_recon_edges and g_recon_edges have Nedge(tree)
                   dir = args$output_dir,
                   name = args$output_name,
                   a = args$fdr,
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
                   snp_in_gene = snps_per_gene)

    save_results_as_r_object(args$output_dir, args$output_name, results_object)
  }
  print("end")
}


