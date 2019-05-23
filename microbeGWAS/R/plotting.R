# CONTIUOUS PLOTS
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
  plot_p_recon <- contMap(tr, pheno_vector, method = "user", anc.states = pheno_anc_rec, plot = FALSE)
  plot(plot_p_recon,
       add = TRUE,
       ylim = c(-1 / 25 * Ntip(tr), Ntip(tr)),
       colors = plot_p_recon$cols,
       lwd = 4,
       ftype = "off",
       offset = 1.7)

  # No output ------------------------------------------------------------------
} # end plot_continuous_phenotype()

plot_significant_hits <- function(disc_cont, tr, fdr, dir, name, pval_all_transition, pheno_vector, annot, perm, results_all_trans, pheno_anc_rec, geno_reconstruction, geno_confidence, geno_transition, geno, pheno_recon_ordered_by_edges, tr_and_pheno_hi_conf, all_trans_sig_hits){
  # Function description -------------------------------------------------------
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
  for (i in 1:length(geno_transition)){
    trans_edge_mat <- cbind(trans_edge_mat, geno_transition[[i]]$transition)
  }
  colnames(trans_edge_mat) <- colnames(geno)

  # TODO Update trans_edge_mat to exclude low confidence  edges, currently it includes transition edges (all high and some low confidence transitions)
  for (c in 1:ncol(trans_edge_mat)){
    trans_edge_mat[(1:Nedge(tr))[geno_confidence[[c]] == 0], c] <- NA
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

  row.names(p_trans_mat) <- row.names(trans_edge_mat) <- c(1:Nedge(tr))
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

  make_manhattan_plot(dir, name, pval_all_transition$hit_pvals, fdr, "transition")
  cell_width_value <- 1.5
  if (ncol(ordered_by_p_val) < 50){
    cell_width_value <- 10
  }
  pheatmap( # Plot the heatmap
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

  pheatmap( # Plot the heatmap
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
  for (j in 1:nrow(pval_all_transition$hit_pvals)){
    if (pval_all_transition$hit_pvals[j, 1] < fdr){
      print("making significant plots")
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

      pheatmap(
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
  for (i in 1:length(geno_transition)){
    trans_dir_edge_mat <- cbind(trans_dir_edge_mat, geno_transition[[i]]$trans_dir)
  }

  colnames(trans_dir_edge_mat) <- colnames(geno)

  # TODO Update trans_edge_mat to exclude low confidence  edges, currently it includes transition edges (all high and some low confidence transitions)
  for (c in 1:ncol(trans_dir_edge_mat)){
    trans_dir_edge_mat[(1:Nedge(tr))[geno_confidence[[c]] == 0], c] <- NA
  }

  all_tables <- all_lists <- rep(list(NULL), ncol(trans_dir_edge_mat))
  delta_pheno_table <- matrix(0, nrow = 3, ncol = 1)
  row.names(delta_pheno_table) <- c("geno_parent_0_child_1", "geno_parent_1_child_0", "geno_no_change")
  colnames(delta_pheno_table) <- c("sum(|delta_phenotype|)")
  for (i in 1:ncol(trans_dir_edge_mat)){
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
  for (i in 1:ncol(trans_dir_edge_mat)){
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

histogram_raw_high_confidence_delta_pheno_highlight_transition_edges <- function(geno_transition, geno_confidence, pheno_recon_ordered_by_edges, tr, index, non_trans_color, trans_color){
  # Function description ------------------------------------------------------#
  # Plot a historgram of delta phenotype, distinguishing transition and non-transition edges.
  #
  # Inputs:
  # TODO
  # geno_transition
  # geno_confidence
  # pheno_recon_ordered_by_edges
  # tr
  # index
  # non_trans_color
  # trans_color
  #
  # Output:
  # None.
  #
  # Check inputs ---------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_is_number(index)
  check_is_string(non_trans_color)
  check_is_string(trans_color)
  # TODO Add more checks

  # Function -------------------------------------------------------------------
  trans_index     <- c(1:Nedge(tr))[as.logical(geno_transition[[index]]$transition)]
  non_trans_index <- c(1:Nedge(tr))[!geno_transition[[index]]$transition]
  hi_conf_trans_index     <- trans_index[    as.logical(geno_confidence[[index]][trans_index    ])]
  hi_conf_non_trans_index <- non_trans_index[as.logical(geno_confidence[[index]][non_trans_index])]
  raw_trans_delta     <- calc_raw_diff(hi_conf_trans_index,     pheno_recon_ordered_by_edges)
  raw_non_trans_delta <- calc_raw_diff(hi_conf_non_trans_index, pheno_recon_ordered_by_edges)

  hist(raw_non_trans_delta,
       main = paste("Raw delta phenotype on only high confidence edges\n # trans edge= ",
                    length(raw_trans_delta), "\n# non trans edge =",
                    length(raw_non_trans_delta), sep = ""),
       breaks = Nedge(tr)/4,
       col = non_trans_color,
       border = FALSE,
       xlim = c(min(raw_trans_delta, raw_non_trans_delta), max(raw_trans_delta, raw_non_trans_delta)),
       ylab = "Count",
       cex = .8,
       xlab = "Raw delta phenotype")

  hist(raw_trans_delta,
       breaks = Nedge(tr)/4,
       col = trans_color,
       border = trans_color,
       add = TRUE,
       xlim = c(min(raw_trans_delta, raw_non_trans_delta), max(raw_trans_delta, raw_non_trans_delta)))

  # No outputs -----------------------------------------------------------------
}

histogram_abs_high_confidence_delta_pheno_highlight_transition_edges <- function(results_all_trans, tr, index, non_trans_color, trans_color){
  # Function description ------------------------------------------------------#
  # Plot a histogram of absolute value of delta phenotype, distinguishing transition and non-transition edges.
  #
  # Inputs:
  # TODO
  # results_all_trans
  # tr
  # index
  # non_trans_color
  # trans_color
  #
  # Output:
  # None.
  #
  # Check inputs ---------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_is_number(index)
  check_is_string(non_trans_color)
  check_is_string(trans_color)
  # TODO Add more checks

  # Function -------------------------------------------------------------------
  par(mar = c(4, 4, 4, 4))
  hist(results_all_trans$observed_pheno_non_trans_delta[[index]],
       breaks = Nedge(tr)/4,
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
       breaks = Nedge(tr)/4,
       col = trans_color,
       border = trans_color,
       add = TRUE,
       xlim = c(0, max(results_all_trans$observed_pheno_trans_delta[[index]], results_all_trans$observed_pheno_non_trans_delta[[index]])))

  # No output ------------------------------------------------------------------
} # end histogram_abs_high_confidence_delta_pheno_highlight_transition_edges()

histogram_all_delta_pheno_overlaid_with_high_conf_delta_pheno <- function(p_trans_mat, geno_confidence, tr, index){
  # Function description ------------------------------------------------------#
  # Plot a histogram of absolute value of delta phenotype, distinguishing transition and non-transition edges.
  #
  # Inputs:
  # TODO
  # p_trans_mat
  # geno_confidence
  # tr
  # index
  #
  # Output:
  # None.
  #
  # Check inputs ---------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_is_number(index)
  # TODO Add more checks

  # Function -------------------------------------------------------------------
  edge_num <- length(unlist(p_trans_mat))
  hi_edge_num <- length(unlist(p_trans_mat)[as.logical(geno_confidence[[index]])])
  title <- paste("|delta phenotype| on all edges \n Light Green: all edges = ", edge_num, "\n Grey: high confidence edges = ", hi_edge_num, sep = "")
  delta_phenotype_on_all_edges <- as.numeric(unlist(p_trans_mat))
  hist(delta_phenotype_on_all_edges,
       breaks = Nedge(tr)/4,
       col = rgb(0, 0.5, 0, 0.25),
       border = FALSE,
       ylab = "Count",
       xlab = "Delta phenotype",
       main = title)

  delta_phenotype_on_high_confidence_edges <- as.numeric(unlist(p_trans_mat))[as.logical(geno_confidence[[index]])]
  hist(delta_phenotype_on_high_confidence_edges, # plot phenotype transition only high confidence edges for this genotype
       breaks = Nedge(tr)/4,
       col = rgb(0, 0, 0, 0.25),
       border = FALSE,
       ylab = "Count",
       add = TRUE)

  # No output ------------------------------------------------------------------
} # end histogram_all_delta_pheno_overlaid_with_high_conf_delta_pheno()

make_manhattan_plot <- function(outdir, geno_pheno_name, pval_hits, alpha, trans_or_recon){
  # Function description -------------------------------------------------------
  # Create a file name and save results to that file name.
  #
  # Inputs:
  # hits.        Vector. Pvals.
  # output_dir.  Character.
  # output_name. Character.
  # pval_name.   Character.
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
  sig_temp <-subset(neg_log_p_with_num, neg_log_p_with_num[ , 2] > -log(alpha))
  ymax <- max(-log(0.01), neg_log_p_with_num[ , 2, drop = TRUE])
  #pdf(paste0(outdir, "/phyc_", trans_or_recon, "_", geno_pheno_name, "_manhattan_plot.pdf"))
  with(neg_log_p_with_num,
       plot(x = neg_log_p_with_num[ , 1],
            y = neg_log_p_with_num[ , 2, drop = TRUE],
            type = "p",
            main = paste(trans_or_recon, "phyC", geno_pheno_name, sep = " "),
            col = rgb(0, 0, 0, 0.3),
            pch = 19,
            xlab = "Genetic locus",
            ylim = c(0, ymax),
            ylab = "-ln(p-val)" ))

  abline(h = -log(alpha), col = "red")
  if (nrow(sig_temp) > 0){
    text(x = sig_temp[ , 1], y = sig_temp[ , 2], labels = row.names(sig_temp), pos = 1, cex = 0.7)
  }
  #dev.off()
} #end make_manhattan_plot()

plot_tree_with_colored_edges <- function(tr, edges_to_highlight, geno_confidence, edge_color_na, edge_color_bright, title, annot, trans_or_recon, index){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
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
  edge_color <- rep("black", Nedge(tr))
  if (trans_or_recon == "recon"){
    edge_color[edges_to_highlight[[index]] == 1] <- edge_color_bright
  } else if (trans_or_recon == "trans"){
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
  tiplabels(pch = 21,
            col = annot[ , 2],
            adj = 2,
            bg = annot[ , 2],
            cex = 0.75)
  if (!is.null(annot)){
    legend("bottomleft",
           legend = unique(annot[ , 1]),
           col = unique(annot[ , 2]),
           lty = 1,
           ncol = length(unique(annot[ , 1])),
           lwd = 5,
           cex = 0.6)
  }
} # end plot_tree_with_colored_edges()


# DISCRETE:
plot_sig_hits_summary <- function(heat_tr, tr, g_mat, p_mat, annot, sig_hits, heatmap_title, filename_start, high_conf_edges, recon_or_trans, short, pheno_conf, bootstrap, pheno_recon_or_trans_by_edge, all_high_conf_edges, geno_trans_or_recon, output_name){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
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
  hits <- g_mat[ , colnames(g_mat) %in%  row.names(sig_hits), drop = FALSE]

  if (recon_or_trans == "trans"){
    filename <- paste(filename_start, "tree.pdf", sep = "")
  } else {
    filename <- paste(filename_start, "tree.pdf", sep = "")
  }
  pdf(filename, width = 16, height = 20)
  par(mfrow = c(1, 2))
  par(mgp = c(3, 1, 0))
  par(oma = c(0, 0, 4, 0))
  # 1. High confidence phenotype (reconstruction or transition)
  edge_color <- pheno_recon_or_trans_by_edge
  edge_color[edge_color == 1] <- "green"
  edge_color[high_conf_edges == 0] <- "grey"
  edge_color[edge_color == 0] <- "black"
  if(recon_or_trans == "trans"){
    title <- paste("Phenotype transitions\n Black=no trans Green=trans Grey=low conf", sep = "")
  } else if (recon_or_trans == "recon"){
    title <- paste("Phenotype reconstruction\n Black=absent Blue=present Grey=low conf", sep = "")
  }
  plot(tr,
       main = title,
       type = "phylogram",
       use.edge.length = TRUE,
       edge.width = 2.0,
       cex = 0.4,
       cex.main = 1.0,
       no.margin = FALSE,
       label.offset = 0.0001,
       edge.color = edge_color)

  # 2. Low confidence edges shown (Low confidence due to either: ML phenotype, bootstrap values, and/or long tree edges)
  edge_color <- rep("black", Nedge(tr))
  edge_color[pheno_conf == 0] <- "blue"
  edge_color[bootstrap == 0] <- "purple"
  edge_color[short == 0] <- "red"
  title <- paste("Low confidence edges\nRed=long Purple=bootstrap Blue=pheno\nBlack=high conf all metrics", sep = "")
  plot(tr,
       main = title,
       type = "phylogram",
       use.edge.length = TRUE,
       edge.width = 2.0,
       cex = 0.4,
       cex.main = 1.0,
       no.margin = FALSE,
       label.offset = 0.0001,
       edge.color = edge_color)


  #mtext(output_name,
  #      outer = TRUE,
  #      side = 3,
  #      cex = 1.2,
  #      line = 1)
  #dev.off()

  # Now 2 PDF per significant hit: 1 is the heatmap and 1 is the geno & pheno trans/recon

  counter <- 0
  for (i in 1:ncol(g_mat)){
    # Subset to only significant hits:
    if (colnames(g_mat)[i] %in% row.names(sig_hits)){
      counter <- counter + 1
      hit_name <- colnames(g_mat)[i]
      # First heatmap
      single_hit <- g_mat[ , i, drop = FALSE]
      # Second/Third: Geno recon or trans & Pheno recon or trans
      #if (recon_or_trans == "trans"){
      #  filename <- paste(filename_start, "tree_sig_hit_", counter , ".pdf", sep = "")
      #} else {
      #  filename <- paste(filename_start, "tree_sig_hit_", counter , ".pdf", sep = "")
      #}
      #pdf(filename, width = 16, height = 20)
      par(mfrow = c(1, 2))
      par(mgp = c(3, 1, 0))
      par(oma = c(0, 0, 4, 0))

      # Genotype
      edge_color <- geno_trans_or_recon[[i]]
      edge_color[edge_color == 1] <- "orange"
      edge_color[edge_color == 0] <- "black"
      edge_color[all_high_conf_edges[[i]] == 0] <- "grey"

      if (recon_or_trans == "trans"){
        title <- paste("Genotype transition edges\nBlack=no trans Orange=trans\nGrey=low confidence", sep = "")
      } else {
        title <- paste("Genotype reconstruction\nBlack=absent Orange=present\nGrey=low confidence", sep = "")
      }
      plot(tr,
           main = title,
           type = "phylogram",
           use.edge.length = TRUE,
           edge.width = 2.0,
           cex = 0.5,
           cex.main = 1.0,
           no.margin = FALSE,
           label.offset = 0.0001,
           edge.color = edge_color)

      # Phenotype
      edge_color <- pheno_recon_or_trans_by_edge
      edge_color[edge_color == 1] <- "green"
      edge_color[edge_color == 0] <- "black"
      edge_color[all_high_conf_edges[[i]] == 0] <- "grey"
      if(recon_or_trans == "trans"){
        title <- paste("Phenotype transition edges\n Black=no trans Green=trans\nGrey=low confidence", sep = "")
      } else {
        title <- paste("Phenotype reconstruction\n Black=absent Green=present\nGrey=low confidence", sep = "")
      }
      plot(tr,
           main = title,
           type = "phylogram",
           use.edge.length = TRUE,
           edge.width = 2.0,
           cex = 0.5,
           cex.main = 1.0,
           no.margin = FALSE,
           label.offset = 0.0001,
           edge.color = edge_color)

      # mtext(paste("Sig hit: ", hit_name, sep = ""),
      #        outer = TRUE,
      #        side = 3,
      #        cex = 1.2,
      #        line = 1)
      #  dev.off()
    } # end subset on significant
  } # end for loop
} # end plot_sig_hits_summary()

create_heatmap_compatible_tree <- function(tree){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
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
  heatmap_tree <- tree
  heatmap_tree$edge.length[which(heatmap_tree$edge.length == 0)] <- 0.00001
  heatmap_tree <- chronopl(heatmap_tree,
                           lambda = 0.1,
                           tol = 0)
  heatmap_tree <- as.dendrogram(as.hclust.phylo(heatmap_tree))
  return(heatmap_tree)
} # end create_heatmap_compatible_tree()


discrete_plots <- function(tr, dir, name, fdr,
                           annot, num_perm, recon_hit_vals,
                           trans_hit_vals, p_recon_edges,
                           g_recon_edges, pheno_anc_rec,
                           recon_perm_obs_results, trans_perm_obs_results,
                           tr_and_pheno_hi_conf, geno_confidence,
                           g_trans_edges, p_trans_edges, snp_in_gene){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
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

  pdf(paste0(dir, "/phyc_",  name, ".pdf"))
  # reconstruction first
  par(mfrow = c(1,1))
  make_manhattan_plot(dir, name, recon_hit_vals, fdr, "reconstruction")

  # TODO 2019-03-18 fix all references to reconstruction a genotype transition-- because that's what it should be
  g_recon_mat <- matrix(0, nrow = Nedge(tr), ncol = length(g_recon_edges))

  for (i in 1:length(g_recon_edges)){
    g_recon_mat[ , i] <- g_recon_edges[[i]]
    g_recon_mat[geno_confidence[[i]] == 0, i] <- NA
  }

  p_recon_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
  p_mat <- matrix(p_recon_edges, nrow = length(p_recon_edges), ncol = 1)
  colnames(p_mat) <- "pheno_presence"
  phenotype_annotation <- as.data.frame(p_mat)
  row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)

  temp_g_recon_mat <- cbind(phenotype_annotation, g_recon_mat)
  g_recon_mat <- temp_g_recon_mat[order(temp_g_recon_mat[,1], na.last = FALSE, decreasing = FALSE ), 2:ncol(temp_g_recon_mat), drop = FALSE]


  cell_width_value <- 1.5
  if (ncol(g_recon_mat) < 50){
    cell_width_value <- 10
  }

  colnames(g_recon_mat) <- row.names(recon_hit_vals)

  significant_loci <- data.frame("locus" = rep("not_sig", ncol(g_recon_mat)), stringsAsFactors = FALSE)
  row.names(significant_loci) <- row.names(recon_hit_vals)
  log_p_value <- data.frame(-log(recon_hit_vals))
  significant_loci[log_p_value > -log(fdr)] <- "sig"

  if (!is.null(snp_in_gene)){
    snp_in_gene <- as.data.frame(snp_in_gene, row.names = 1)
    colnames(snp_in_gene) <- "SNPs in gene"
    column_annot <- cbind(significant_loci, log_p_value, snp_in_gene) # TODO add a test to make sure order doesn't matter here
  } else {
    column_annot <- cbind(significant_loci, log_p_value)
  }


  ann_colors = list(
    locus = c(not_sig = "white", sig = "blue"),
    pheno_presence = c( na = "grey", absence = "white", presence = "red")
  )

  ordered_by_p_val      <-           g_recon_mat[ , match(row.names(log_p_value)[order(log_p_value[ , 1])], colnames(g_recon_mat))]
  column_annot_ordered_by_p_val <-     column_annot[match(row.names(log_p_value)[order(log_p_value[ , 1])], row.names(column_annot)), ]

  # reconstruction loci summary heat maps
  pheatmap( # Plot the heatmap
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

  # loop through reconstruction sig hits:
  pheno_as_list <- list(p_recon_edges)
  pheno_conf_as_list <- list(tr_and_pheno_hi_conf)
  # TODO break these plots into more functions b/c lots of redundant code between recon and transition plots
  for (j in 1:nrow(recon_hit_vals)){
    if (recon_hit_vals[j, 1] < fdr){
      par(mfrow = c(3, 2), mgp = c(3, 1, 0), oma = c(0, 0, 4, 0), mar = c(4, 4, 4, 4))
      # pheno
      plot_tree_with_colored_edges(tr, pheno_as_list, pheno_conf_as_list, "grey", "red", paste0("\n Phenotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", 1)
      # geno
      plot_tree_with_colored_edges(tr, g_recon_edges, geno_confidence, "grey", "red", paste0(row.names(recon_hit_vals)[j], "\n Genotype transition:\n Red = transition; Black = no transition"), annot, "recon", j)
      # Permutation test
      max_x <- max(recon_perm_obs_results$permuted_count[[j]], recon_perm_obs_results$observed_overlap[j]) # change to loop through sig hits
      hist(recon_perm_obs_results$permuted_count[[j]],
           breaks = num_perm/10,
           xlim = c(0, max_x),
           col = "grey",
           border = FALSE,
           ylab = "Count",
           xlab = "# edges where genotype-phenotype co-occur",
           main = paste0("Overlap of genotype transition edge\n& phenotype presence \npval=", round(recon_hit_vals[j, 1], 4), "\nRed=observed,Grey=permutations", sep = ""))
      abline(v = recon_perm_obs_results$observed_overlap[j], col = "red")

      p_recon_edges[tr_and_pheno_hi_conf == 0] <- -1 # should be NA but it won't work correctedly TODO
      p_mat <- matrix(p_recon_edges, nrow = length(p_recon_edges), ncol = 1)
      colnames(p_mat) <- "pheno_presence"
      phenotype_annotation <- as.data.frame(p_mat)
      row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)


      temp_g_recon_edges <- g_recon_edges[[j]]
      temp_g_recon_edges[geno_confidence[[j]] == 0] <- NA
      g_mat <- as.matrix(temp_g_recon_edges)
      row.names(g_mat) <- c(1:nrow(g_mat))
      colnames(g_mat) <- "genotype_transition"
      temp_g_mat <- cbind(g_mat, phenotype_annotation)
      g_mat <- temp_g_mat[order(temp_g_mat[,2], temp_g_mat[,1], na.last = FALSE, decreasing = FALSE ), 1, drop = FALSE]
      plotting_logical <- check_if_g_mat_can_be_plotted(g_mat)

      if (plotting_logical){
        ann_colors <- make_ann_colors(g_mat)

        if (!is.null(snp_in_gene)){
          num_snps <- snp_in_gene[row.names(recon_hit_vals)[j], , drop = FALSE]
          row.names(num_snps) <- "genotype_transition"
          colnames(num_snps) <- "SNPs_in_gene"
          ann_colors <- c(ann_colors, list(SNPs_in_gene = c(num_snps_in_gene = "blue")))
        } else {
          num_snps <- NULL
        }

        pheatmap(
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
          cellwidth         = 20)
      }
    }
  }

  par(mfrow = c(1,1))
  make_manhattan_plot(dir, name, trans_hit_vals, fdr, "transition")

  # start transition heatmaps
  g_trans_mat <- matrix(0, nrow = Nedge(tr), ncol = length(g_recon_edges))

  for (i in 1:length(g_recon_edges)){
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


  if (!is.null(snp_in_gene)){
    column_annot <- cbind(significant_loci, log_p_value, snp_in_gene) # TODO add a test to make sure order doesn't matter here
  } else {
    column_annot <- cbind(significant_loci, log_p_value)
  }

  ann_colors = list(
    locus = c(not_sig = "white", sig = "blue"),
    pheno_transitions = c( na = "grey", no_transition = "white", transition = "red")
  )

  ordered_by_p_val      <-           g_trans_mat[ , match(row.names(log_p_value)[order(log_p_value[ , 1])], colnames(g_trans_mat))]
  column_annot_ordered_by_p_val <-     column_annot[match(row.names(log_p_value)[order(log_p_value[ , 1])], row.names(column_annot)), ]

  print(ordered_by_p_val)
  print(column_annot_ordered_by_p_val)
  print(phenotype_annotation)

  # Transition loci summary heat maps
  pheatmap( # Plot the heatmap
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

  # loop through transition sig hits:
  for (j in 1:nrow(trans_hit_vals)){
    if (trans_hit_vals[j, 1] < fdr){
      par(mfrow = c(3, 2), mgp = c(3, 1, 0), oma = c(0, 0, 4, 0), mar = c(4, 4, 4, 4))
      # Plot pheno
      p_trans_edges_as_list <- list(p_trans_edges)
      plot_tree_with_colored_edges(tr, pheno_as_list,         pheno_conf_as_list, "grey", "red", paste0("\n Phenotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", 1)
      plot_tree_with_colored_edges(tr, p_trans_edges_as_list, pheno_conf_as_list, "grey", "red", paste0("\n Phenotype transitions:\n Red = transition; Black = no change"), annot, "recon", 1)
      # Plot geno
      plot_tree_with_colored_edges(tr, g_recon_edges, geno_confidence, "grey", "red", paste0(row.names(recon_hit_vals)[j], "\n Genotype reconstruction:\n Red = Variant; Black = WT"), annot, "recon", j)
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
      g_mat<- temp_g_mat[order(temp_g_mat[,2], temp_g_mat[,1], na.last = FALSE, decreasing = FALSE ), 1, drop = FALSE]

      ann_colors = list(pheno_transition = c( na = "grey", no_transition = "white", transition = "red"))

      plotting_logical <- check_if_g_mat_can_be_plotted(g_mat)


      if (plotting_logical){
        ann_colors <- make_ann_colors(g_mat)
        if (!is.null(snp_in_gene)){
          num_snps <- snp_in_gene[row.names(trans_hit_vals)[j], , drop = FALSE]
          row.names(num_snps) <- "genotype_transition"
          colnames(num_snps) <- "SNPs_in_gene"
          ann_colors <- c(ann_colors, list(SNPs_in_gene = c(num_snps_in_gene = "blue")))
        }

        pheatmap(
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
          legend = FALSE,
          cellwidth         = 20)
      }
    }
  }
  dev.off()
} # end discrete_plots()


make_ann_colors <- function(geno_matrix){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
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

  if (ones + zeros + nas == 3){
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

