#' Plot the continuous phenotype on tree with the phenotype coloring the branches.
#'
#' @param tr
#' @param pheno_vector
#' @param pheno_anc_rec
#'
#' @return
#'
#' @examples
plot_continuous_phenotype <- function(tr, pheno_vector, pheno_anc_rec){
  # Function description -------------------------------------------------------
  # TODO
  # Plot the continuous phenotype as a tree with phenotype as colors in branches.
  #
  # Inputs:
  # tr. Phylo.
  # pheno_vector. Vector.
  # pheno_anc_rec. Vector?
  #
  # Output:
  # None
  #
  # Check inputs ---------------------------------------------------------------
  check_tree_is_valid(tr)

  # Function -------------------------------------------------------------------
  plot_p_recon <- phytools::contMap(tr, pheno_vector, method = "user", anc.states = pheno_anc_rec, plot = FALSE)
  plot(plot_p_recon,
       add = TRUE,
       ylim = c(-1 / 25 * ape::Ntip(tr), ape::Ntip(tr)),
       colors = plot_p_recon$cols,
       lwd = 4,
       ftype = "off",
       offset = 1.7)

  # No output ------------------------------------------------------------------
} # end plot_continuous_phenotype()


histogram_raw_high_confidence_delta_pheno_highlight_transition_edges <- function(geno_transition, geno_confidence, pheno_recon_ordered_by_edges, tr, index, non_trans_color, trans_color){
  trans_index     <- c(1:ape::Nedge(tr))[as.logical(geno_transition[[index]]$transition)]
  non_trans_index <- c(1:ape::Nedge(tr))[!geno_transition[[index]]$transition]
  hi_conf_trans_index     <- trans_index[    as.logical(geno_confidence[[index]][trans_index    ])]
  hi_conf_non_trans_index <- non_trans_index[as.logical(geno_confidence[[index]][non_trans_index])]
  #hi_conf_pheno_trans_delta     <- calculate_phenotype_change_on_edge(hi_conf_trans_index,     pheno_recon_ordered_by_edges)
  #hi_conf_pheno_non_trans_delta <- calculate_phenotype_change_on_edge(hi_conf_non_trans_index, pheno_recon_ordered_by_edges)
  raw_trans_delta     <- calc_raw_diff(hi_conf_trans_index,     pheno_recon_ordered_by_edges)
  raw_non_trans_delta <- calc_raw_diff(hi_conf_non_trans_index, pheno_recon_ordered_by_edges)

  hist(raw_non_trans_delta,
       main = paste("Raw delta phenotype on only high confidence edges\n # trans edge= ",
                    length(raw_trans_delta), "\n# non trans edge =",
                    length(raw_non_trans_delta), sep = ""),
       breaks = ape::Nedge(tr)/4,
       col = non_trans_color,
       border = FALSE,
       xlim = c(min(raw_trans_delta, raw_non_trans_delta), max(raw_trans_delta, raw_non_trans_delta)),
       ylab = "Count",
       cex = .8,
       xlab = "Raw delta phenotype")

  hist(raw_trans_delta,
       breaks = ape::Nedge(tr)/4,
       col = trans_color,
       border = trans_color,
       add = TRUE,
       xlim = c(min(raw_trans_delta, raw_non_trans_delta), max(raw_trans_delta, raw_non_trans_delta)))
}

histogram_abs_high_confidence_delta_pheno_highlight_transition_edges <- function(results_all_trans, tr, index, non_trans_color, trans_color){
  par(mar = c(4, 4, 4, 4))
  hist(results_all_trans$observed_pheno_non_trans_delta[[index]],
       breaks = ape::Nedge(tr)/4,
       col = non_trans_color,
       border = FALSE,
       ylab = "Count",
       xlab = "|Delta phenotype|",
       xlim = c(0, max(results_all_trans$observed_pheno_trans_delta[[index]], results_all_trans$observed_pheno_non_trans_delta[[index]])),
       cex = .8,
       main = paste("|delta phenotype| on only high confidence edges \n # non-trans edges = ",
                    length(results_all_trans$observed_pheno_non_trans_delta[[index]]),
                    "\n # trans edges = ", length(results_all_trans$observed_pheno_trans_delta[[index]]), sep = ""))


  hist(results_all_trans$observed_pheno_trans_delta[[index]],
       breaks = ape::Nedge(tr)/4,
       col = trans_color,
       border = trans_color,
       add = TRUE,
       xlim = c(0, max(results_all_trans$observed_pheno_trans_delta[[index]], results_all_trans$observed_pheno_non_trans_delta[[index]])))
} # end histogram_abs_high_confidence_delta_pheno_highlight_transition_edges()

histogram_all_delta_pheno_overlaid_with_high_conf_delta_pheno <- function(p_trans_mat, geno_confidence, tr, index){
  edge_num <- length(unlist(p_trans_mat))
  hi_edge_num <- length(unlist(p_trans_mat)[as.logical(geno_confidence[[index]])])
  title <- paste("|delta phenotype| on all edges \n Light Green: all edges = ", edge_num, "\n Grey: high confidence edges = ", hi_edge_num, sep = "")
  delta_phenotype_on_all_edges <- as.numeric(unlist(p_trans_mat))
  hist(delta_phenotype_on_all_edges,
       breaks = ape::Nedge(tr)/4,
       col = rgb(0, 0.5, 0, 0.25),
       border = FALSE,
       ylab = "Count",
       xlab = "Delta phenotype",
       main = title)

  delta_phenotype_on_high_confidence_edges <- as.numeric(unlist(p_trans_mat))[as.logical(geno_confidence[[index]])]
  hist(delta_phenotype_on_high_confidence_edges, # plot phenotype transition only high confidence edges for this genotype
       breaks = ape::Nedge(tr)/4,
       col = rgb(0, 0, 0, 0.25),
       border = FALSE,
       ylab = "Count",
       add = TRUE)
} # end histogram_all_delta_pheno_overlaid_with_high_conf_delta_pheno()




#' Plot continuous phenotype results.
#'
#' @param disc_cont
#' @param tr
#' @param fdr
#' @param dir
#' @param name
#' @param pval_all_transition
#' @param pheno_vector
#' @param annot
#' @param perm
#' @param results_all_trans
#' @param pheno_anc_rec
#' @param geno_reconstruction
#' @param geno_confidence
#' @param geno_transition
#' @param geno
#' @param pheno_recon_ordered_by_edges
#' @param tr_and_pheno_hi_conf
#' @param all_trans_sig_hits
#'
#' @return
#'
#' @examples
plot_significant_hits <- function(disc_cont, tr, fdr, dir, name, pval_all_transition, pheno_vector, annot, perm, results_all_trans, pheno_anc_rec, geno_reconstruction, geno_confidence, geno_transition, geno, pheno_recon_ordered_by_edges, tr_and_pheno_hi_conf, all_trans_sig_hits){
  # Function description -------------------------------------------------------
  # Plot continuous phenotype results.
  # Plot all genotype's -log(p-value) as a manhattan plot.
  # If any hits are significant after FDR, then plot heatmaps & trees for each of these.
  #
  # Inputs:
  # tr.                  Phylo.
  # fdr.                 Number. False discovery rate.
  # dir.                 Character. Output path.
  # name.                Character. Output name.
  # pheno_vector.        Vector.
  # annot.               Matrix.
  # perm.                Number.

  # Output:
  # None.
  #
  # Check inputs ---------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  check_num_between_0_and_1(fdr)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_if_vector(pheno_vector)
  check_if_permutation_num_valid(perm)

  # Function -------------------------------------------------------------------
  trans_edge_mat <- NULL
  for (i in 1:length(geno_transition)) {
    trans_edge_mat <- cbind(trans_edge_mat, geno_transition[[i]]$transition)
  }
  colnames(trans_edge_mat) <- colnames(geno)

  # TODO Update trans_edge_mat to exclude low confidence  edges, currently it includes transition edges (all high and some low confidence transitions)
  for (c in 1:ncol(trans_edge_mat)) {
    trans_edge_mat[(1:ape::Nedge(tr))[geno_confidence[[c]] == 0], c] <- NA
  }
  # end update trans_edge_mat
  ph_trans <- abs(pheno_recon_ordered_by_edges[ , 1] - pheno_recon_ordered_by_edges[ , 2])

  p_trans_mat <- matrix(ph_trans, nrow = length(ph_trans), ncol = 1)
  colnames(p_trans_mat) <- "delta_pheno"
  p_trans_mat <- as.data.frame(round(p_trans_mat, 2))

  significant_loci <- data.frame("locus" = rep("not_sig", ncol(trans_edge_mat)), stringsAsFactors = FALSE)
  row.names(significant_loci) <- colnames(trans_edge_mat)
  significant_loci[row.names(significant_loci) %in% row.names(all_trans_sig_hits), ] <- "sig"

  log_p_value <- data.frame(-log(pval_all_transition$hit_pvals))
  column_annot <- cbind(significant_loci, log_p_value)

  row.names(p_trans_mat) <- row.names(trans_edge_mat) <- c(1:ape::Nedge(tr))
  # end heatmap prep
  ann_colors = list(
    locus = c(not_sig = "white", sig = "blue")
  )

  sorted_trans_edge_mat <-           trans_edge_mat[match(row.names(p_trans_mat)[order(p_trans_mat[ , 1])], row.names(trans_edge_mat)), ]
  ordered_by_p_val      <- sorted_trans_edge_mat[ , match(row.names(log_p_value)[order(log_p_value[ , 1])], colnames(sorted_trans_edge_mat))]
  column_annot_ordered_by_p_val <-     column_annot[match(row.names(log_p_value)[order(log_p_value[ , 1])], row.names(column_annot)), ]
  colnames(column_annot) <- colnames(column_annot_ordered_by_p_val) <- c("locus", "-ln(p-val)")

  fname <- create_file_name(dir, name, paste("summary_and_sig_hit_results.pdf", sep = "")) # 2019-02-25 removed counter from pdf file name because combining.
  pdf(fname, width = 16, height = 20)

  make_manhattan_plot(dir, name, pval_all_transition$hit_pvals, fdr, "continuous")
  cell_width_value <- 1.5
  if (ncol(ordered_by_p_val) < 50) {
    cell_width_value <- 10
  }
  pheatmap::pheatmap( # Plot the heatmap
    ordered_by_p_val,
    main          = paste0("Edges:\n hi conf trans vs delta pheno"),
    cluster_cols  = FALSE,
    na_col = "grey",
    cluster_rows  = FALSE,
    show_rownames = FALSE,
    color = c("white", "black"),
    annotation_col = column_annot_ordered_by_p_val,
    annotation_row = p_trans_mat,
    annotation_colors = ann_colors,
    show_colnames = TRUE,
    cellwidth = cell_width_value)

  pheatmap::pheatmap( # Plot the heatmap
    sorted_trans_edge_mat,
    main          = paste0("Edges:\n hi conf trans vs delta pheno"),
    cluster_cols  = TRUE,
    cluster_rows  = FALSE,
    show_rownames = FALSE,
    color = c("white", "black"),
    annotation_col = column_annot,
    annotation_row = p_trans_mat,
    annotation_colors = ann_colors,
    show_colnames = TRUE,
    cellwidth = cell_width_value,
    na_col = "grey")

  # ONLY MAKE THE FOLLOWING PLOTS FOR SIGNIFICANT LOCI
  counter <- 0 # TODO can I get rid of counter now? 2019-02-25
  for (j in 1:nrow(pval_all_transition$hit_pvals)) {
    if (pval_all_transition$hit_pvals[j, 1] < fdr) {
      counter <- counter + 1
      par(mfrow = c(3, 3))
      par(mgp   = c(3, 1, 0))
      par(oma   = c(0, 0, 4, 0))
      par(mar = c(4, 4, 4, 4))

      plot_continuous_phenotype(tr, pheno_vector, pheno_anc_rec)
      plot_tree_with_colored_edges(tr, geno_reconstruction, geno_confidence, "grey", "red", paste0(row.names(pval_all_transition$hit_pvals)[j], "\n Genotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", j)
      plot_tree_with_colored_edges(tr, geno_transition,     geno_confidence, "grey", "red", "Genotype transition edge:\n Red = transition; Black = No transition", annot, "trans", j)
      histogram_all_delta_pheno_overlaid_with_high_conf_delta_pheno(p_trans_mat, geno_confidence, tr, j)
      histogram_abs_high_confidence_delta_pheno_highlight_transition_edges(results_all_trans, tr, j, "grey", "red") # 6. Delta phenotype histogram. Note that results_all_trans$observed_pheno_non_trans_delta[[j]] is already subset down to only high confidence edges
      histogram_raw_high_confidence_delta_pheno_highlight_transition_edges(geno_transition, geno_confidence, pheno_recon_ordered_by_edges, tr, j, "grey", "red")

      hist(log(results_all_trans$ks_statistics[[j]]),
           breaks = perm/10, col = "grey", border = FALSE,
           main = paste("Null distribution of KS statistic for all transitions.\n Red = Observed KS statistic.\n p-value = ", round(pval_all_transition$hit_pvals[j, 1], 10), "\np-value rank = ", rank(pval_all_transition$hit_pvals)[j], sep = ""),
           ylab = "Count",
           xlab = "ln(KS statistic)",
           xlim = c(min(log(as.numeric(results_all_trans$observed_ks_stat[j])), log(results_all_trans$ks_statistics[[j]])), 0))
      abline(v = log(as.numeric(results_all_trans$observed_ks_stat[j])), col = "red")

      pheatmap::pheatmap(
        sorted_trans_edge_mat[ , j, drop = FALSE],
        main          = paste0(row.names(pval_all_transition$hit_pvals)[j], "\n Tree edges: hi conf trans vs delta pheno"),
        cluster_cols  = FALSE,
        cluster_rows  = FALSE,
        na_col = "grey",
        show_rownames = FALSE,
        color = c("white", "black"),
        annotation_row = p_trans_mat,
        show_colnames = TRUE,
        cellwidth = 20)
    }
  }

  dev.off()

  trans_dir_edge_mat <- NULL
  for (i in 1:length(geno_transition)) {
    trans_dir_edge_mat <- cbind(trans_dir_edge_mat, geno_transition[[i]]$trans_dir)
  }

  colnames(trans_dir_edge_mat) <- colnames(geno)

  # TODO Update trans_edge_mat to exclude low confidence  edges, currently it includes transition edges (all high and some low confidence transitions)
  for (c in 1:ncol(trans_dir_edge_mat)) {
    trans_dir_edge_mat[(1:ape::Nedge(tr))[geno_confidence[[c]] == 0], c] <- NA
  }

  all_tables <- all_lists <- rep(list(NULL), ncol(trans_dir_edge_mat))
  delta_pheno_table <- matrix(0, nrow = 3, ncol = 1)
  row.names(delta_pheno_table) <- c("geno_parent_0_child_1", "geno_parent_1_child_0", "geno_no_change")
  colnames(delta_pheno_table) <- c("sum(|delta_phenotype|)")
  for (i in 1:ncol(trans_dir_edge_mat)) {
    temp_table <- delta_pheno_table
    temp_table[1, 1] <- sum(p_trans_mat[which(trans_dir_edge_mat[ , i] == 1),  1], na.rm = TRUE)
    temp_table[2, 1] <- sum(p_trans_mat[which(trans_dir_edge_mat[ , i] == -1), 1], na.rm = TRUE)
    temp_table[3, 1] <- sum(p_trans_mat[which(trans_dir_edge_mat[ , i] == 0),  1], na.rm = TRUE)
    all_tables[[i]] <- temp_table
  }
  names(all_tables) <- colnames(geno)
  delta_pheno_table <- all_tables
  delta_pheno_list <- rep(list(0), 3)
  names(delta_pheno_list) <- c("geno_parent_0_child_1", "geno_parent_1_child_0", "geno_no_change")
  for (i in 1:ncol(trans_dir_edge_mat)) {
    temp_list <- delta_pheno_list
    temp_list[[1]] <- p_trans_mat[which(trans_dir_edge_mat[ , i] == 1),  1]
    temp_list[[2]] <- p_trans_mat[which(trans_dir_edge_mat[ , i] == -1), 1]
    temp_list[[3]] <- p_trans_mat[which(trans_dir_edge_mat[ , i] == 0),  1]
    all_lists[[i]] <- temp_list
    names(all_lists)[i] <- colnames(trans_dir_edge_mat)[i]
  }
  delta_pheno_list <- all_lists

  # Return output --------------------------------------------------------------
  results <- list("trans_dir_edge_mat" = trans_dir_edge_mat,
                  "p_trans_mat" = p_trans_mat,
                  "delta_pheno_table" = delta_pheno_table,
                  "delta_pheno_list" = delta_pheno_list )
  return(results)
} # end plot_significant_hits()

#' Create a manhattan plot of GWAS hit p-values.
#'
#' @param outdir
#' @param geno_pheno_name
#' @param pval_hits
#' @param alpha
#' @param trans_or_recon
#'
#' @return
#' @export
#'
#' @examples
make_manhattan_plot <- function(outdir, geno_pheno_name, pval_hits, alpha, trans_or_recon){
  # Function description -------------------------------------------------------
  # Create a manhattan plot of GWAS hit p-values.
  # TODO
  # Inputs:
  #
  # Output:
  # None
  #
  # Check inputs ---------------------------------------------------------------

  # Function -------------------------------------------------------------------
  # Create negative log p-values with arbitrary locus numbers
  neg_log_p_value <- data.frame(-log(pval_hits))
  neg_log_p_with_num <- cbind(1:nrow(neg_log_p_value), neg_log_p_value)
  colnames(neg_log_p_with_num)[1] <- "locus"
  sig_temp <- subset(neg_log_p_with_num, neg_log_p_with_num[ , 2] > -log(alpha))
  ymax <- max(-log(0.01), neg_log_p_with_num[ , 2, drop = TRUE])
  with(neg_log_p_with_num,
       plot(x = neg_log_p_with_num[ , 1],
            y = neg_log_p_with_num[ , 2, drop = TRUE],
            type = "p",
            main = paste(trans_or_recon, geno_pheno_name, sep = " "),
            col = rgb(0, 0, 0, 0.3),
            pch = 19,
            xlab = "Genetic locus",
            ylim = c(0, ymax),
            ylab = "-ln(p-val)" ))

  abline(h = -log(alpha), col = "red")
  if (nrow(sig_temp) > 0) {
    text(x = sig_temp[ , 1], y = sig_temp[ , 2], labels = row.names(sig_temp), pos = 1, cex = 0.7)
  }
  #dev.off()
} #end make_manhattan_plot()

#' Plot a phylogenetic tree with certain edges highlighted.
#'
#' @param tr
#' @param edges_to_highlight
#' @param geno_confidence
#' @param edge_color_na
#' @param edge_color_bright
#' @param title
#' @param annot
#' @param trans_or_recon
#' @param index
#'
#' @return
#'
#' @examples
plot_tree_with_colored_edges <- function(tr, edges_to_highlight, geno_confidence, edge_color_na, edge_color_bright, title, annot, trans_or_recon, index){
  # Function description -------------------------------------------------------
  # Plot a phylogenetic tree with certain edges highlighted.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  edge_color <- rep("black", ape::Nedge(tr))
  if (trans_or_recon == "recon") {
    edge_color[edges_to_highlight[[index]] == 1] <- edge_color_bright
  } else if (trans_or_recon == "trans") {
    edge_color[edges_to_highlight[[index]]$transition == 1] <- edge_color_bright
  }
  edge_color[geno_confidence[[index]] == 0] <- edge_color_na # grey out long edges and low ML bootstrap support
  par(mar = c(4, 4, 4, 4))
  plot(tr,
       font = 1,
       edge.color = edge_color,
       main = title,
       use.edge.length = FALSE,
       label.offset = 3,
       adj = 0)
  ape::tiplabels(pch = 21,
            col = annot[ , 2],
            adj = 2,
            bg = annot[ , 2],
            cex = 0.75)
  if (!is.null(annot)) {
    legend("bottomleft",
           legend = unique(annot[ , 1]),
           col = unique(annot[ , 2]),
           lty = 1,
           ncol = length(unique(annot[ , 1])),
           lwd = 5,
           cex = 0.6)
  }
} # end plot_tree_with_colored_edges()

#' Create object to annotate columns in the significant hit results.
#'
#' @param geno_matrix
#'
#' @return
#' @export
#'
#' @examples
make_ann_colors <- function(geno_matrix){
  # Function description -------------------------------------------------------
  # Create object to annotate columns in the significant hit results.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  ones <- sum(geno_matrix == 1, na.rm = TRUE) > 1
  zeros <- sum(geno_matrix == 0, na.rm = TRUE) > 1
  nas <- sum(is.na(geno_matrix)) > 1

  if (ones + zeros + nas == 3) {
    ann_colors = list(pheno_presence = c(na = "grey", absent = "white", present = "red"))
  } else if (ones == 1 && zeros == 1 && nas == 0) {
    ann_colors = list(pheno_presence = c(absent = "white", present = "red"))
  } else if (ones == 1 && zeros == 0 && nas == 1) {
    ann_colors = list(pheno_presence = c(na = "grey", present = "red"))
  } else if (ones == 0 && zeros == 1 && nas == 1) {
    ann_colors = list(pheno_presence = c(na = "grey", absent = "white"))
  } else if (ones == 0 && zeros == 0 && nas == 1) {
    ann_colors = list(pheno_presence = c(na = "grey"))
  } else if (ones == 0 && zeros == 1 && nas == 0) {
    ann_colors = list(pheno_presence = c(absent = "white"))
  } else if (ones == 1 && zeros == 0 && nas == 0) {
    ann_colors = list(pheno_presence = c(present = "red"))
  } else {
    stop("no ones, zeroes, or NAs present in g_mat")
  }
  return(ann_colors)
}

#' Plot convergence test results.
#'
#' @param tr
#' @param dir
#' @param name
#' @param fdr
#' @param annot
#' @param num_perm
#' @param recon_hit_vals
#' @param p_recon_edges
#' @param recon_perm_obs_results
#' @param tr_and_pheno_hi_conf
#' @param geno_confidence
#' @param g_trans_edges
#' @param p_trans_edges
#' @param snp_in_gene
#'
#' @return
#' @export
#'
#' @examples
discrete_plot_orig <- function(tr, dir, name, fdr, annot, num_perm,
                                recon_hit_vals, p_recon_edges,
                                recon_perm_obs_results, tr_and_pheno_hi_conf,
                                geno_confidence, g_trans_edges, p_trans_edges,
                                snp_in_gene){
  # Function description -------------------------------------------------------
  # Plot the discrete test results.
  #
  # Inputs:
  # tr. Phylo.
  # dir. Directory where to save plots.
  # name. Prefix in plot file name.
  # fdr. Numeric. False discovery rate. Between 0 and 1.
  # TODO add description of annotation and make sure that including an annotation doesn't break anything OR remove annotation.
  # annot. ??
  # num_perm. Numeric. Number of permutations.
  # recon_hit_vals. Dataframe. Nrows = number of genotypes. Ncol = 1. Corrected p-values for each genotype tested.
  # p_recon_edges. Vector. Length = Nedge(tree). Reconstruction of phenotype.
  # recon_perm_obs_results. List of many results.
  #     $hit_pvals. Character. P-val for each genotype. Length = number of tested genotypes.
  #     $permuted_count. List of vectors. 1 vector for each tested genotype. Length of each vector = number of permuations.
  #     $observed_overlap. Integer vector. Length = number of tested genotypes.
  # tr_and_pheno_hi_conf. Vector of logicals. TRUE = high confidence. FALSE = low confidence. Length = Nedge(tree).
  # geno_confidence. List of vectors. Length of individual vector = Nedge(tree). Genotype high confidence edges. Either 1 (high confidence) or 0 (low confidence).
  # g_trans_edges. List of vectors. Length of individual vector = Nedge(tree). Genotype transition edges. Either 1 (transition) or 0 (no transition).
  # p_trans_edges. Vector. Length = Nedge(tree). Transitions marked as 1, not transition marked as 0.
  # snp_in_gene. Either NULL or Table of integers where each entry corresponds to one genotype.
  #
  # Outputs:
  # Plots printed into one pdf.
  #
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_num_between_0_and_1(fdr)
  check_if_permutation_num_valid(num_perm)
  if (ncol(recon_hit_vals) != 1 | nrow(recon_hit_vals) != length(geno_confidence)) {
    stop("Dimensions of hit p-values dataframe are incorrect.")
  }
  check_equal(length(p_recon_edges), ape::Nedge(tr))
  if (!is.null(snp_in_gene)) {
    if (class(snp_in_gene) != "table" | typeof(snp_in_gene) != "integer") {
      stop("snp_in_gene should be a table of integers")
    }
  }
  check_equal(length(p_trans_edges), ape::Nedge(tr))
  check_equal(length(geno_confidence[[1]]), ape::Nedge(tr))
  check_equal(length(g_trans_edges[[1]]), ape::Nedge(tr))
  check_if_binary_vector(geno_confidence[[1]])
  check_if_binary_vector(p_trans_edges)
  check_if_binary_vector(g_trans_edges[[1]])
  check_equal(length(tr_and_pheno_hi_conf), ape::Nedge(tr))
  check_equal(length(recon_perm_obs_results$permuted_count[[1]]), num_perm)
  check_class(recon_perm_obs_results$hit_pvals, "character")
  check_class(recon_perm_obs_results$observed_overlap, "integer")

  # Function -------------------------------------------------------------------
  image_width <- 250
  pdf(paste0(dir, "/phyc_",  name, ".pdf"))

  par(mfrow = c(1,1))
  make_manhattan_plot(dir, name, recon_hit_vals, fdr, "phyc")

  g_trans_mat <- matrix(0, nrow = ape::Nedge(tr), ncol = length(g_trans_edges))

  for (i in 1:length(g_trans_edges)) {
    g_trans_mat[ , i] <- g_trans_edges[[i]]
    g_trans_mat[geno_confidence[[i]] == 0, i] <- NA
  }

  p_recon_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
  p_mat <- matrix(p_recon_edges, nrow = length(p_recon_edges), ncol = 1)
  colnames(p_mat) <- "pheno_presence"
  phenotype_annotation <- as.data.frame(p_mat)
  row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)

  temp_g_trans_mat <- cbind(phenotype_annotation, g_trans_mat)
  g_trans_mat <- temp_g_trans_mat[order(temp_g_trans_mat[,1], na.last = FALSE, decreasing = FALSE ), 2:ncol(temp_g_trans_mat), drop = FALSE]
  colnames(g_trans_mat) <- row.names(recon_hit_vals)

  cell_width_value <- image_width / ncol(g_trans_mat)

  significant_loci <- data.frame("locus" = rep("not_sig", ncol(g_trans_mat)), stringsAsFactors = FALSE)
  row.names(significant_loci) <- row.names(recon_hit_vals)
  log_p_value <- data.frame(-log(recon_hit_vals))
  significant_loci[log_p_value > -log(fdr)] <- "sig"

  if (!is.null(snp_in_gene)) {
    snp_in_gene <- as.data.frame(snp_in_gene, row.names = 1)
    colnames(snp_in_gene) <- "SNPs in gene"
    snp_in_gene <- snp_in_gene[row.names(snp_in_gene) %in% row.names(log_p_value), , drop = FALSE]
    column_annot <- cbind(significant_loci, log_p_value, snp_in_gene) # TODO add a test to make sure order doesn't matter here
  } else {
    column_annot <- cbind(significant_loci, log_p_value)
  }

  ordered_by_p_val      <-           g_trans_mat[ , match(row.names(log_p_value)[order(log_p_value[ , 1])], colnames(g_trans_mat)), drop = FALSE]
  column_annot_ordered_by_p_val <-     column_annot[match(row.names(log_p_value)[order(log_p_value[ , 1])], row.names(column_annot)), ]

  # TODO: make sure this kind of annotation steps gets generalized to more locatoins -- I think I wrote a function like this somewhere else. Double check!
  if (length(unique(phenotype_annotation[ ,1])) == 3) {
    pheno_presence_col = c( na = "grey", absence = "white", presence = "red")
  } else if (length(unique(phenotype_annotation[ ,1])) == 2) {
    if (sum(unique(phenotype_annotation[ ,1]) %in% c(-1, 0)) == 2) {
      pheno_presence_col = c(na = "grey", absence = "white")
    } else if (sum(unique(phenotype_annotation[ ,1]) %in% c(-1, 1)) == 2) {
      pheno_presence_col = c( na = "grey", presence = "red")
    } else if (sum(unique(phenotype_annotation[ ,1]) %in% c(1, 0)) == 2) {
      pheno_presence_col = c(absence = "white", presence = "red")
    }
  } else if (length(unique(phenotype_annotation[ ,1])) == 1) {
    pheno_presence_col = c(presence = "red")
    if (unique(phenotype_annotation[ ,1]) == -1) {
      pheno_presence_col = c(na = "grey")
    } else if (unique(phenotype_annotation[ ,1]) == 0) {
      pheno_presence_col = c(absence = "white")
    }
  }

  if (length(unique(column_annot_ordered_by_p_val[ ,1])) == 2) {
    locus_col = c(not_sig = "white", sig = "blue")
  } else if (length(unique(column_annot_ordered_by_p_val[ ,1])) == 1) {
    locus_col = c(sig = "blue")
    if (unique(column_annot_ordered_by_p_val[ ,1]) == "not_sig") {
      locus_col = c(not_sig = "white")
    }
  }
  ann_colors = list(locus = locus_col, pheno_presence = pheno_presence_col)

  plotting_logical <- check_if_g_mat_can_be_plotted(ordered_by_p_val)
  if (plotting_logical) {
    # phyc loci summary heat maps
    pheatmap::pheatmap( # Plot the heatmap
      ordered_by_p_val,
      main          = paste0("Edges:\n Genotype transition with phenotype presence/absence"),
      cluster_cols  = FALSE,
      na_col = "grey",
      cluster_rows  = FALSE,
      show_rownames = FALSE,
      color = c("white", "black"),
      annotation_col = column_annot_ordered_by_p_val,
      annotation_row = phenotype_annotation,
      annotation_colors = ann_colors,
      show_colnames = TRUE,
      cellwidth = cell_width_value)
  }

  # loop through reconstruction sig hits:
  pheno_as_list <- list(p_recon_edges)
  pheno_conf_as_list <- list(tr_and_pheno_hi_conf)
  # TODO break these plots into more functions b/c lots of redundant code between recon and transition plots
  for (j in 1:nrow(recon_hit_vals)) {
    if (recon_hit_vals[j, 1] < fdr) {
      par(mfrow = c(3, 2), mgp = c(3, 1, 0), oma = c(0, 0, 4, 0), mar = c(4, 4, 4, 4))
      # pheno
      plot_tree_with_colored_edges(tr, pheno_as_list, pheno_conf_as_list, "grey", "red", paste0("\n Phenotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", 1)
      # geno
      plot_tree_with_colored_edges(tr, g_trans_edges, geno_confidence, "grey", "red", paste0(row.names(recon_hit_vals)[j], "\n Genotype transition:\n Red = transition; Black = no transition"), annot, "recon", j)
      # Permutation test
      max_x <- max(recon_perm_obs_results$permuted_count[[j]], recon_perm_obs_results$observed_overlap[j]) # change to loop through sig hits
      hist(recon_perm_obs_results$permuted_count[[j]],
           breaks = num_perm/10,
           xlim = c(0, max_x),
           col = "grey",
           border = FALSE,
           ylab = "Count",
           xlab = "# edges with genotype transition & phenotype presence",
           main = paste0("Phyc: Overlap of genotype transition edge\n& phenotype presence \npval=", round(recon_hit_vals[j, 1], 4), "\nRed=observed,Grey=permutations", sep = ""))
      abline(v = recon_perm_obs_results$observed_overlap[j], col = "red")

      p_recon_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
      p_mat <- matrix(p_recon_edges, nrow = length(p_recon_edges), ncol = 1)
      colnames(p_mat) <- "pheno_presence"
      phenotype_annotation <- as.data.frame(p_mat)
      row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)
      temp_g_trans_edges <- g_trans_edges[[j]]
      temp_g_trans_edges[geno_confidence[[j]] == 0] <- NA
      g_mat <- as.matrix(temp_g_trans_edges)
      row.names(g_mat) <- c(1:nrow(g_mat))
      colnames(g_mat) <- "genotype_transition"
      temp_g_mat <- cbind(g_mat, phenotype_annotation)
      g_mat <- temp_g_mat[order(temp_g_mat[,2], temp_g_mat[,1], na.last = FALSE, decreasing = FALSE ), 1, drop = FALSE]
      plotting_logical <- check_if_g_mat_can_be_plotted(g_mat)

      if (plotting_logical) {
        ann_colors <- make_ann_colors(g_mat)

        if (!is.null(snp_in_gene)) {
          num_snps <- snp_in_gene[row.names(recon_hit_vals)[j], , drop = FALSE]
          row.names(num_snps) <- "genotype_transition"
          colnames(num_snps) <- "# grouped genotypes"
          ann_colors <- c(ann_colors, list(SNPs_in_gene = c(num_snps_in_gene = "blue")))
        } else {
          num_snps <- NULL
        }

        cell_width_value <- image_width / ncol(g_mat)

        pheatmap::pheatmap(
          mat               = g_mat,
          main              = paste0(row.names(recon_hit_vals)[j], "\n Tree edges clustered by edge type\n Genotype transition edge\n & phenotype present edge"),
          cluster_cols      = FALSE,
          cluster_rows      = FALSE,
          na_col            = "grey",
          show_rownames     = FALSE,
          color             = c("white", "red"),
          annotation_row    = phenotype_annotation,
          annotation_col    = num_snps,
          annotation_legend = TRUE,
          annotation_colors = ann_colors,
          show_colnames     = TRUE,
          legend            = TRUE,
          cellwidth         = cell_width_value)
      }
    }
  }

  dev.off()
} # end discrete_plot_orig()

#' Plot the synchronous test results.
#'
#' @param tr
#' @param dir
#' @param name
#' @param fdr
#' @param annot
#' @param num_perm
#' @param trans_hit_vals
#' @param trans_perm_obs_results
#' @param tr_and_pheno_hi_conf
#' @param geno_confidence
#' @param g_trans_edges
#' @param p_trans_edges
#' @param snp_in_gene
#'
#' @return
#'
#' @examples
discrete_plot_trans  <- function(tr, dir, name, fdr, annot, num_perm,
                                 trans_hit_vals, trans_perm_obs_results,
                                 tr_and_pheno_hi_conf, geno_confidence,
                                 g_trans_edges, p_trans_edges, snp_in_gene){
  # Function description -------------------------------------------------------
  # Plot the discrete test results.
  #
  # Inputs:
  # tr. Phylo.
  # dir. Directory where to save plots.
  # name. Prefix in plot file name.
  # fdr. Numeric. False discovery rate. Between 0 and 1.
  # TODO add description of annotation and make sure that including an annotation doesn't break anything OR remove annotation.
  # annot. ??
  # num_perm. Numeric. Number of permutations.
  # trans_hit_vals. Dataframe. Nrows = number of genotypes. Ncol = 1. Corrected p-values for each genotype tested.
  # trans_perm_obs_results. List of many results.
  #     $hit_pvals. Character. P-val for each genotype. Length = number of tested genotypes.
  #     $permuted_count. List of vectors. 1 vector for each tested genotype. Length of each vector = number of permuations.
  #     $observed_overlap. Integer vector. Length = number of tested genotypes.
  # tr_and_pheno_hi_conf. Vector of logicals. TRUE = high confidence. FALSE = low confidence. Length = Nedge(tree).
  # geno_confidence. List of vectors. Length of individual vector = Nedge(tree). Genotype high confidence edges. Either 1 (high confidence) or 0 (low confidence).
  # g_trans_edges. List of vectors. Length of individual vector = Nedge(tree). Genotype transition edges. Either 1 (transition) or 0 (no transition).
  # p_trans_edges. Vector. Length = Nedge(tree). Transitions marked as 1, not transition marked as 0.
  # snp_in_gene. Either NULL or Table of integers where each entry corresponds to one genotype.
  #
  # Outputs:
  # Plots printed into one pdf.
  #
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_num_between_0_and_1(fdr)
  check_if_permutation_num_valid(num_perm)
  if (ncol(trans_hit_vals) != 1 | nrow(trans_hit_vals) != length(geno_confidence)) {
    stop("Dimensions of hit p-values dataframe are incorrect.")
  }
  if (!is.null(snp_in_gene)) {
    if (class(snp_in_gene) != "table" | typeof(snp_in_gene) != "integer") {
      stop("snp_in_gene should be a table of integers")
    }
  }
  check_equal(length(p_trans_edges), ape::Nedge(tr))
  check_equal(length(geno_confidence[[1]]), ape::Nedge(tr))
  check_equal(length(g_trans_edges[[1]]), ape::Nedge(tr))
  check_if_binary_vector(geno_confidence[[1]])
  check_if_binary_vector(p_trans_edges)
  check_if_binary_vector(g_trans_edges[[1]])
  check_equal(length(tr_and_pheno_hi_conf), ape::Nedge(tr))
  check_equal(length(trans_perm_obs_results$permuted_count[[1]]), num_perm)
  check_class(trans_perm_obs_results$hit_pvals, "character")
  check_class(trans_perm_obs_results$observed_overlap,"integer")

  # Function -------------------------------------------------------------------
  image_width <- 250
  pdf(paste0(dir, "/synchronous_",  name, ".pdf"))

  par(mfrow = c(1,1))
  make_manhattan_plot(dir, name, trans_hit_vals, fdr, "synchronous")

  # heatmaps
  g_trans_mat <- matrix(0, nrow = ape::Nedge(tr), ncol = length(g_trans_edges))

  for (i in 1:length(g_trans_edges)) {
    g_trans_mat[ , i] <- g_trans_edges[[i]]
    g_trans_mat[geno_confidence[[i]] == 0, i] <- NA
  }

  p_trans_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
  p_mat <- matrix(p_trans_edges, nrow = length(p_trans_edges), ncol = 1)
  colnames(p_mat) <- "pheno_transitions"
  phenotype_annotation <- as.data.frame(p_mat)
  row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)

  temp_g_trans_mat <- cbind(phenotype_annotation, g_trans_mat)
  g_trans_mat <- temp_g_trans_mat[order(temp_g_trans_mat[,1], na.last = FALSE, decreasing = FALSE ), 2:ncol(temp_g_trans_mat), drop = FALSE]
  colnames(g_trans_mat) <- row.names(trans_hit_vals)

  significant_loci <- data.frame("locus" = rep("not_sig", ncol(g_trans_mat)), stringsAsFactors = FALSE)
  row.names(significant_loci) <- row.names(trans_hit_vals)
  log_p_value <- data.frame(-log(trans_hit_vals))
  significant_loci[log_p_value > -log(fdr)] <- "sig"

  if (!is.null(snp_in_gene)) {
    snp_in_gene <- as.data.frame(snp_in_gene, row.names = 1)
    colnames(snp_in_gene) <- "# grouped genotypes"
    snp_in_gene <- snp_in_gene[row.names(snp_in_gene) %in% row.names(log_p_value), , drop = FALSE]
    column_annot <- cbind(significant_loci, log_p_value, snp_in_gene) # TODO add a test to make sure order doesn't matter here
  } else {
    column_annot <- cbind(significant_loci, log_p_value)
  }

  ordered_by_p_val      <-           g_trans_mat[ , match(row.names(log_p_value)[order(log_p_value[ , 1])], colnames(g_trans_mat)), drop = FALSE]
  column_annot_ordered_by_p_val <-     column_annot[match(row.names(log_p_value)[order(log_p_value[ , 1])], row.names(column_annot)), ]

  if (length(unique(phenotype_annotation[ ,1])) == 3) {
    pheno_presence_col = c("-1" = "grey", "0" = "white", "1" = "red")
  } else if (length(unique(phenotype_annotation[ ,1])) == 2) {
    if (sum(unique(phenotype_annotation[ ,1]) %in% c(-1, 0)) == 2) {
      pheno_presence_col = c("-1" = "grey", "0" = "white")
    } else if (sum(unique(phenotype_annotation[ ,1]) %in% c(-1, 1)) == 2) {
      pheno_presence_col = c("-1" = "grey", "1" = "red")
    } else if (sum(unique(phenotype_annotation[ ,1]) %in% c(1, 0)) == 2) {
      pheno_presence_col = c("0" = "white", "1" = "red")
    }
  } else if (length(unique(phenotype_annotation[ ,1])) == 1) {
    pheno_presence_col = c("1" = "red")
    if (unique(phenotype_annotation[ ,1]) == -1) {
      pheno_presence_col = c("-1" = "grey")
    } else if (unique(phenotype_annotation[ ,1]) == 0) {
      pheno_presence_col = c("0" = "white")
    }
  }

  if (length(unique(column_annot_ordered_by_p_val[ ,1])) == 2) {
    locus_col = c(not_sig = "white", sig = "blue")
  } else if (length(unique(column_annot_ordered_by_p_val[ ,1])) == 1) {
    locus_col = c(sig = "blue")
    if (unique(column_annot_ordered_by_p_val[ ,1]) == "not_sig") {
      locus_col = c(not_sig = "white")
    }
  }

  ann_colors = list(locus = locus_col, pheno_presence = pheno_presence_col)
  can_be_plotted <- check_if_g_mat_can_be_plotted(ordered_by_p_val)
  if (can_be_plotted) {
    cell_width_value <- image_width / ncol(ordered_by_p_val)

    # Transition loci summary heat maps
    pheatmap::pheatmap( # Plot the heatmap
      ordered_by_p_val,
      main          = paste0("Edges:\n Genotype transitions with phenotype transitions"),
      cluster_cols  = FALSE,
      na_col = "grey",
      cluster_rows  = FALSE,
      show_rownames = FALSE,
      color = c("white", "black"),
      annotation_col = column_annot_ordered_by_p_val,
      annotation_row = phenotype_annotation,
      annotation_colors = ann_colors,
      show_colnames = TRUE,
      cellwidth = cell_width_value)
  }
  pheno_conf_as_list <- list(tr_and_pheno_hi_conf)
  # loop through transition sig hits:
  for (j in 1:nrow(trans_hit_vals)) {
    if (trans_hit_vals[j, 1] < fdr) {
      par(mfrow = c(3, 2), mgp = c(3, 1, 0), oma = c(0, 0, 4, 0), mar = c(4, 4, 4, 4))
      # Plot pheno
      p_trans_edges_as_list <- list(p_trans_edges)
      #plot_tree_with_colored_edges(tr, pheno_as_list,         pheno_conf_as_list, "grey", "red", paste0("\n Phenotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", 1)
      plot_tree_with_colored_edges(tr, p_trans_edges_as_list, pheno_conf_as_list, "grey", "red", paste0("\n Phenotype transitions:\n Red = transition; Black = no change"), annot, "recon", 1)
      # Plot geno
      #plot_tree_with_colored_edges(tr, g_trans_edges, geno_confidence, "grey", "red", paste0(row.names(recon_hit_vals)[j], "\n Genotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", j)
      plot_tree_with_colored_edges(tr, g_trans_edges, geno_confidence, "grey", "red", paste0(row.names(trans_hit_vals)[j], "\n Genotype transitions:\n Red = transition; Black = no change"), annot, "recon", j)
      # Permutation test
      max_x <- max(trans_perm_obs_results$permuted_count[[j]], trans_perm_obs_results$observed_overlap[j]) # change to loop through sig hits
      hist(trans_perm_obs_results$permuted_count[[j]],
           breaks = num_perm/10,
           xlim = c(0, max_x),
           col = "grey",
           border = FALSE,
           ylab = "Count",
           xlab = "# edges where genotype-phenotype transitions co-occur",
           main = paste0("Geno & pheno transition overlap\npval=", round(trans_hit_vals[j, 1], 4), "\nRed=observed,Grey=permutations", sep = "")) # TODO add rank pvalue
      abline(v = trans_perm_obs_results$observed_overlap[j], col = "red")

      # edge heatmap - heatmap is tree edges, annotation is phenotype edges
      par(mfrow = c(1,1))
      p_trans_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
      p_mat <- matrix(p_trans_edges, nrow = length(p_trans_edges), ncol = 1)
      colnames(p_mat) <- "pheno_transition"
      phenotype_annotation <- as.data.frame(p_mat)
      row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)

      temp_g_trans_edges <- g_trans_edges[[j]]
      temp_g_trans_edges[geno_confidence[[j]] == 0] <- NA
      g_mat <- as.matrix(temp_g_trans_edges)
      row.names(g_mat) <- c(1:nrow(g_mat))
      colnames(g_mat) <- "genotype_transition"
      temp_g_mat <- cbind(g_mat, phenotype_annotation)
      g_mat <- temp_g_mat[order(temp_g_mat[ ,2], temp_g_mat[,1], na.last = FALSE, decreasing = FALSE ), 1, drop = FALSE]

      ann_colors = list(pheno_transition = c( na = "grey", no_transition = "white", transition = "red"))
      plotting_logical <- check_if_g_mat_can_be_plotted(g_mat)

      if (plotting_logical) {
        ann_colors <- make_ann_colors(g_mat)
        if (!is.null(snp_in_gene)) {
          num_snps <- snp_in_gene[row.names(trans_hit_vals)[j], , drop = FALSE]
          row.names(num_snps) <- "genotype_transition"
          colnames(num_snps) <- "# grouped genotypes"
          ann_colors <- c(ann_colors, list(SNPs_in_gene = c(num_snps_in_gene = "blue")))
        }
        cell_width_value <- image_width / ncol(g_mat)

        pheatmap::pheatmap(
          g_mat,
          main              = paste0(row.names(trans_hit_vals)[j], "\n Tree edges: genotype & phenotype transitions"),
          cluster_cols      = FALSE,
          cluster_rows      = FALSE,
          na_col            = "grey",
          show_rownames     = FALSE,
          color             = c("white", "red"),
          annotation_row    = phenotype_annotation,
          annotation_col    = num_snps,
          annotation_colors = ann_colors,
          annotation_legend = TRUE,
          show_colnames     = TRUE,
          legend            = FALSE,
          cellwidth         = cell_width_value)
      }
    }
  }
  dev.off()
} # end discrete_plot_trans()
# End of script ----------------------------------------------------------------
