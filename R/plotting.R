#' Plot the continuous phenotype as a phylogenetic tree filled in with a color
#' gradient
#'
#' @description Plot the continuous phenotype on tree with the phenotype
#'  coloring the branches.
#'
#' @param tr Phylo
#' @param pheno_vector Vector. Length = Ntip(tr).
#'  Names(pheno_vector) == tr$tip.label.
#' @param pheno_anc_rec Vector. Length = Nnode(tr).
#'
#' @return A tree plot where the tree edges are filled in with colors
#'  corresponding to the phenotypic trait values.
#'
#' @noRd
plot_continuous_phenotype <- function(tr, pheno_vector, pheno_anc_rec){
  # Check inputs ---------------------------------------------------------------
  check_tree_is_valid(tr)
  check_equal(length(pheno_vector), ape::Ntip(tr))
  check_equal(length(pheno_anc_rec), ape::Nnode(tr))
  if (!identical(names(pheno_vector), tr$tip.label)) {
    stop("Phenotype vector must have same names as tree tips.")
  }

  # Function -------------------------------------------------------------------
  plot_p_recon <- phytools::contMap(tr,
                                    pheno_vector,
                                    method = "user",
                                    anc.states = pheno_anc_rec,
                                    plot = FALSE)
  graphics::plot(plot_p_recon,
       add = TRUE,
       ylim = c(-1 / 25 * ape::Ntip(tr), ape::Ntip(tr)),
       colors = plot_p_recon$cols,
       lwd = 4,
       ftype = "off",
       offset = 1.7)

} # end plot_continuous_phenotype()

#' Draw histogram of the raw phenotype change on each high confidence tree edge
#'
#' @description Plot a histogram to show the change in phenotype on each tree
#'  edge. Plot phenotype change as raw value change (not absolute value).
#'  Exclude low confidence tree edges.
#'
#' @noRd
#' @param geno_transition List of list of 2 vectors: $transition and $trans_dir.
#'   Each list corresponds to 1 genotype. In $transition vector each entry
#'   corresponds to 1 tree edge.
#' @param geno_confidence List of vectors. Each list corresponds to 1 genotype.
#'   Each vector entry corresponds to 1 tree edge. 1 == high confidence edge.
#'   0 == low confidence edge.
#' @param pheno_recon_ordered_by_edges Matrix. Dimensions: Nedge x 2 matrix.
#'   Entries are the phenotype value at the node, where row is the edge, 1st
#'   column is the parent node and 2nd column is the child node.
#' @param tr Phylo.
#' @param index Integer. Indexes which genotype is being plotted.
#' @param non_trans_color Character. Color of non-transition histogram bars.
#' @param trans_color Character. Color of transition histogram bars.
#'
#' @return Histogram showing the change in phenotype on each high confidence
#'   tree edge. Transition edges are different color than non-transition edges.
#'   Low confidence edges are excluded. Change in phenotype is raw (not absolute
#'   value).
hist_raw_hi_conf_delta_pheno <- function(geno_transition,
                                         geno_confidence,
                                         pheno_recon_ordered_by_edges,
                                         tr,
                                         index,
                                         non_trans_color,
                                         trans_color){
  # Check inputs ---------------------------------------------------------------
  check_is_number(index)
  if (index > length(geno_transition)) {
    stop("Index must correspond to one of the genotypes")
  }
  check_equal(length(geno_transition[[1]]$transition), ape::Nedge(tr))
  check_equal(length(geno_transition), length(geno_confidence))
  check_dimensions(pheno_recon_ordered_by_edges,
                   ape::Nedge(tr),
                   ape::Nedge(tr),
                   2,
                   2)
  check_for_root_and_bootstrap(tr)
  check_is_string(non_trans_color)
  check_is_string(trans_color)
  if (non_trans_color == trans_color) {
    stop("These tree edges need to be different colors")
  }

  # Function -------------------------------------------------------------------
  hist_cex_size <- 1.3

  trans_index <-
    c(1:ape::Nedge(tr))[as.logical(geno_transition[[index]]$transition)]
  non_trans_index <- c(1:ape::Nedge(tr))[!geno_transition[[index]]$transition]
  hi_conf_trans_index <-
    trans_index[as.logical(geno_confidence[[index]][trans_index])]
  hi_conf_non_trans_index <-
    non_trans_index[as.logical(geno_confidence[[index]][non_trans_index])]
  raw_trans_delta <- calc_raw_diff(hi_conf_trans_index,
                                   pheno_recon_ordered_by_edges)
  raw_non_trans_delta <- calc_raw_diff(hi_conf_non_trans_index,
                                       pheno_recon_ordered_by_edges)

  graphics::hist(raw_non_trans_delta,
       main = "Raw Delta Phenotype",
       breaks = ape::Nedge(tr) / 4,
       col = non_trans_color,
       border = FALSE,
       xlim = c(min(raw_trans_delta, raw_non_trans_delta),
                1.2 * max(raw_trans_delta, raw_non_trans_delta)),
       ylim = c(0, 0.75 * sum(length(raw_non_trans_delta),
                              length(raw_trans_delta))),
       ylab = "Count",
       xlab = "Raw Delta Phenotype",
       cex.main = hist_cex_size,
       cex.lab = hist_cex_size,
       cex.axis = hist_cex_size)

  graphics::hist(raw_trans_delta,
       breaks = ape::Nedge(tr) / 4,
       col = trans_color,
       border = trans_color,
       add = TRUE)

  legend("topright",
         title = "Hi Conf. Geno. Edge Type",
         legend = c("Transition", "Non-transition"),
         col = c(trans_color, non_trans_color),
         pch = 15,
         cex = hist_cex_size,
         bg = rgb(0, 0, 0, 0.01))
} # end hist_raw_hi_conf_delta_pheno()

#' Draw histogram of the absolute value of phenotype change on each high
#' confidence tree edge
#'
#' @description Plot a histogram to show the change in phenotype on each tree
#'  edge. Plot phenotype change as absolute value change. Exclude low confidence
#'  tree edges.
#'
#' @noRd
#' @param all_trans List of 2 lists. $observed_pheno_non_trans_delta and
#'  $observed_pheno_trans_delta. Position in each corresponds to one genotype.
#'  Contains informatino about amount of phenotype change.
#' @param tr Phylo.
#' @param index Integer. Indexes which genotype is being plotted.
#' @param non_trans_color Character. Color of non-transition histogram bars.
#' @param trans_color Character. Color of transition histogram bars.
#'
#' @return Histogram showing the change in phenotype on each high confidence
#'   tree edge. Transition edges are different color than non-transition edges.
#'   Low confidence edges are excluded. Change in phenotype is absolute value.
hist_abs_hi_conf_delta_pheno <- function(all_trans,
                                         tr,
                                         index,
                                         non_trans_color,
                                         trans_color){
  # Check inputs ---------------------------------------------------------------
  check_is_number(index)
  if (index > length(all_trans)) {
    stop("Index must correspond to one of the genotypes")
  }
  check_for_root_and_bootstrap(tr)
  check_is_string(non_trans_color)
  check_is_string(trans_color)
  if (non_trans_color == trans_color) {
    stop("These tree edges need to be different colors")
  }

  # Function -------------------------------------------------------------------
  hist_cex_size <- 1.3
  graphics::par(mar = c(4, 4, 4, 4))
  graphics::hist(all_trans$observed_pheno_non_trans_delta[[index]],
       breaks = ape::Nedge(tr) / 4,
       col = non_trans_color,
       border = FALSE,
       ylab = "Count",
       xlab = "|Delta Phenotype|",
       xlim = c(0,
                1.2 * max(all_trans$observed_pheno_trans_delta[[index]],
                          all_trans$observed_pheno_non_trans_delta[[index]])),
       main = "|Delta Phenotype|",
       ylim = c(0, 0.75 * sum(length(all_trans$observed_pheno_trans_delta[[index]]),
                              length(all_trans$observed_pheno_non_trans_delta[[index]]))),
       cex.main = hist_cex_size,
       cex.lab = hist_cex_size,
       cex.axis = hist_cex_size)

  graphics::hist(all_trans$observed_pheno_trans_delta[[index]],
       breaks = ape::Nedge(tr) / 4,
       col = trans_color,
       border = trans_color,
       add = TRUE)

  legend("topright",
         title = "Hi Conf. Geno. Edge Type",
         legend = c("Transition", "Non-transition"),
         col = c(trans_color, non_trans_color),
         pch = 15,
         cex = hist_cex_size,
         bg = rgb(0, 0, 0, 0.01))
} # end hist_abs_hi_conf_delta_pheno()

#' Draw histogram of the absolute value of the phenotype change on each edge
#'
#' @description Plot a histogram to show the change in phenotype on each tree
#'  edge. Plot phenotype change as absolute value change. First, plot all tree
#'  edge values then overlay only the high confidence edges.
#'
#' @noRd
#' @param p_trans_mat Numeric matrix. Dim: ncol == 1, nrow == Nedge(tree).
#'  Values are the absolute value in the change of the phenotype on that edge.
#'   Each row corresponds to one tree edge.
#' @param geno_confidence List of vectors. Each list corresponds to 1 genotype.
#'   Each vector entry corresponds to 1 tree edge. 1 == high confidence edge.
#'   0 == low confidence edge.
#' @param tr Phylo.
#' @param index Integer. Indexes which genotype is being plotted.
#'
#' @return Histogram showing the change in phenotype on each edge. Emphasizes
#'   edges lost due to low confidence. Change in phenotype is absolute value.
hist_abs_delta_pheno_all_edges <- function(p_trans_mat,
                                           geno_confidence,
                                           tr,
                                           index){
  # Check inputs ---------------------------------------------------------------
  check_is_number(index)
  if (index > length(geno_confidence)) {
    stop("Index must correspond to one of the genotypes")
  }
  check_for_root_and_bootstrap(tr)
  check_dimensions(p_trans_mat,
                   exact_rows = ape::Nedge(tr),
                   min_rows = ape::Nedge(tr),
                   exact_cols = 1,
                   min_cols = 1)
  check_equal(length(geno_confidence[[1]]), ape::Nedge(tr))

  # Function -------------------------------------------------------------------
  hist_cex_size <- 1.3
  transparent_purple <- rgb(1, 0, 1, 0.25)
  transparent_blue <- rgb(0, 0, 1, 0.25)

  edge_num <- length(unlist(p_trans_mat))
  hi_edge_num <-
    length(unlist(p_trans_mat)[as.logical(geno_confidence[[index]])])
  delta_phenotype_on_all_edges <- as.numeric(unlist(p_trans_mat))
  graphics::hist(delta_phenotype_on_all_edges,
       breaks = ape::Nedge(tr) / 4,
       col = transparent_blue,
       border = FALSE,
       ylab = "Count",
       xlab = "|Delta Phenotype|",
       main = "|Delta Phenotype|",
       cex.main = hist_cex_size,
       cex.lab = hist_cex_size,
       cex.axis = hist_cex_size,
       xlim = c(0, 1.2 * max(delta_phenotype_on_all_edges)),
       ylim = c(0, 0.75 * (length(delta_phenotype_on_all_edges))))

  delta_pheno_on_hi_conf_edges <-
    as.numeric(unlist(p_trans_mat))[as.logical(geno_confidence[[index]])]
  graphics::hist(delta_pheno_on_hi_conf_edges,
  # plot phenotype transition only high confidence edges for this genotype
       breaks = ape::Nedge(tr) / 4,
       col = transparent_purple,
       border = FALSE,
       ylab = "Count",
       add = TRUE)

  legend("topright",
         title = "Edge Type",
         legend = c("Hi Confidence", "Any"),
         col = c(transparent_purple, transparent_blue),
         pch = 15,
         cex = hist_cex_size,
         bg = rgb(0, 0, 0, 0.01))
} # end hist_abs_delta_pheno_all_edges()

#' Plot all results from the Continuous Test
#'
#' @details Plot all genotype's -log(p-value) as a manhattan plot. If any hits
#' are significant after FDR, then plot histrograms & trees for each of these.
#'
#' @noRd
#' @param disc_cont String, should be "continuous."
#' @param tr Phylo.
#' @param fdr Number. False discovery rate.
#' @param dir Character. Output path.
#' @param name Character. Output name.
#' @param pval_all_transition List of two data.frames.
#'  * $hit_pvals. Dataframe. 1 column. Nrow = number of genotypes. Row.names =
#'       genotypes. Column name = "fdr_corrected_pvals". Values between 1 and 0.
#'  * $sig_pvals. Dataframe. 1 column. Nrow = number of genotypes that are
#'       significant after FDR correction. Column name =
#'       "fdr_corrected_pvals[fdr_corrected_pvals < fdr]".
#'       Row.names = genotypes. Nrow = is variable-- could be between 0 and
#'       max number of genotypes. It will only have rows if the corrected
#'       p-value is less than the fdr value.
#' @param pheno_vector Vector. Length = Ntip(tr).
#' @param perm Number.
#' @param results_all_trans List of 8.
#'  * $pvals. Named numeric vector. Length == number of genotypes. Values
#'      between 1 and 0. Names are genotype names.
#'  * $ks_statistics. List of numeric vectors. Length of list == number of
#'      genotypes. Each vector has length == number of permutations. Values
#'      between 1 and 0.
#'  * $observed_pheno_trans_delta. List of numeric vectors. Length of list ==
#'      number of genotypes. Vectors are of variable length because length is
#'      the number of transition edges for that particular genotype. Vectors are
#'      numeric.
#'  * $observed_pheno_non_trans_delta. List of numeric vectors. Length of list
#'       == number of genotypes. Vectors are of variable length because length
#'       is the number of non-transition edges for that particular genotype.
#'       Vectors are numeric.
#'  * $trans_median. Numberic. Vector. Length = number of genotypes. Describes
#'      median delta phenotype on all transition edges.
#'  * $all_edges_median. Numeric vector. Length = number of genotypes. Describes
#'      median delta phenotype on all edges.
#'  * $num_genotypes. Integer. The number of genotypes.
#'  * $observed_ks_stat. Numeric Vector. Length = number of genotypes. Values
#'      between 1 and 0.
#' @param pheno_anc_rec Vector. The values of the ancestral reconstruction of
#'  the phenotype at each internal node. Length = Nnode(tr).
#' @param geno_reconstruction List of lists. Binary. Number of lists = number of
#'   genotypes. Length(each individual list) == Nedge(tree).
#' @param geno_confidence List. Each entry corresponds to one genotype.
#'  Length = number of genotypes.
#' @param geno_transition Object with two lists: $trans_dir and $transition.
#'  Each list has an entry for each genotype. Each sublist has one value for
#'  each tree edge.
#' @param geno Matrix. Columns = genotypes, rows = samples. Binary.
#' @param pheno_recon_ordered_by_edges Matrix. Dim: nrow = Nedge(tr) x ncol = 2.
#'  Parent (older) node is 1st column. Child (younger) node is the 2nd column.
#'  Ancestral reconstruction value of each node.
#' @param tr_and_pheno_hi_conf List. Length(list) = ncol(mat) == number of
#'   genotypes. Each entry is a vector with length == Nedge(tr). All entries are
#'   0 (low confidence) or 1 (high confidence).
#' @param all_trans_sig_hits Data.frame. 1 column. Colnum names =
#'  "fdr_corrected_pvals". Nrow = variable. Number of genotypes that are (1)
#'  significant after multiple test correction and (2) have higher median delta
#'  phenotype on transition edges than on all edges. Values are between 1 and 0.
#'  Rownames are genotypes.
#' @param group_logical Logical. Inidicates whether or not genotypes were
#'  grouped.
#'
#' @return Plots continuous test resutlts.
plot_continuous_results <- function(disc_cont,
                                    tr,
                                    fdr,
                                    dir,
                                    name,
                                    pval_all_transition,
                                    pheno_vector,
                                    perm,
                                    results_all_trans,
                                    pheno_anc_rec,
                                    geno_reconstruction,
                                    geno_confidence,
                                    geno_transition,
                                    geno,
                                    pheno_recon_ordered_by_edges,
                                    tr_and_pheno_hi_conf,
                                    all_trans_sig_hits,
                                    group_logical){
  # Check inputs ---------------------------------------------------------------
  check_str_is_discrete_or_continuous(disc_cont)
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  check_num_between_0_and_1(fdr)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_class(pval_all_transition$hit_pvals, "data.frame")
  check_equal(length(pheno_vector), ape::Ntip(tr))
  check_if_permutation_num_valid(perm)
  check_equal(length(results_all_trans$ks_statistics), ncol(geno))
  check_equal(length(pheno_anc_rec), ape::Nnode(tr))
  check_equal(length(geno_reconstruction), ncol(geno))
  check_equal(length(geno_reconstruction[[1]]), ape::Nedge(tr))
  check_equal(length(geno_confidence), ncol(geno))
  check_equal(length(geno_transition[[1]]$transition), ape::Nedge(tr))
  check_if_binary_matrix(geno)
  check_dimensions(pheno_recon_ordered_by_edges,
                   ape::Nedge(tr),
                   ape::Nedge(tr),
                   2,
                   2)
  check_equal(length(tr_and_pheno_hi_conf), ncol(geno))
  check_class(group_logical, "logical")

  # Function -------------------------------------------------------------------
  image_width <- 250
  hist_cex_size <- 1.3
  trans_edge_mat <- NULL
  for (i in 1:length(geno_transition)) {
    trans_edge_mat <- cbind(trans_edge_mat, geno_transition[[i]]$transition)
  }
  colnames(trans_edge_mat) <- colnames(geno)

  for (c in 1:ncol(trans_edge_mat)) {
    trans_edge_mat[(1:ape::Nedge(tr))[geno_confidence[[c]] == 0], c] <- NA
  }

  ph_trans <-
    abs(pheno_recon_ordered_by_edges[, 1] - pheno_recon_ordered_by_edges[, 2])

  p_trans_mat <- matrix(ph_trans, nrow = length(ph_trans), ncol = 1)
  colnames(p_trans_mat) <- "|Delta Phenotype|"
  p_trans_mat <- as.data.frame(round(p_trans_mat, 2))

  significant_loci <-
    data.frame("Locus Significance" = rep("Not Significant",
                                          ncol(trans_edge_mat)),
                                 stringsAsFactors = FALSE)
  row.names(significant_loci) <- colnames(trans_edge_mat)
  significant_loci[row.names(significant_loci) %in%
                     row.names(all_trans_sig_hits), ] <- "Significant"

  log_p_value <- data.frame(-log(pval_all_transition$hit_pvals))
  column_annot <- cbind(significant_loci, log_p_value)

  row.names(p_trans_mat) <- row.names(trans_edge_mat) <- c(1:ape::Nedge(tr))

  ann_colors <- list(
    `Locus Significance` = c(`Not Significant` = "white", Significant = "blue")
  )

  sorted_trans_edge_mat <-
    trans_edge_mat[match(row.names(p_trans_mat)[order(p_trans_mat[, 1])],
                         row.names(trans_edge_mat)), , drop = FALSE]
  ordered_by_p_val <-
    sorted_trans_edge_mat[, match(
      row.names(log_p_value)[order(log_p_value[, 1])],
      colnames(sorted_trans_edge_mat)), drop = FALSE]

  column_annot_ordered_by_p_val <-
    column_annot[match(row.names(log_p_value)[order(log_p_value[, 1])],
                       row.names(column_annot)), , drop = FALSE]

  if (ncol(column_annot_ordered_by_p_val) == 2) {
    colnames(column_annot) <- colnames(column_annot_ordered_by_p_val) <-
      c("Locus Significance", "-ln(FDR Corrected P-value)")
  } else {
    colnames(column_annot) <- colnames(column_annot_ordered_by_p_val) <-
      c("Locus Significance", "-ln(FDR Corrected P-value)", "Variants in Group")
  }

  if (group_logical) {
    fname <- paste0(dir, "/hogwash_continuous_grouped_", name, ".pdf")
  } else {
    fname <- paste0(dir, "/hogwash_continuous_", name, ".pdf")
  }

  grDevices::pdf(fname, width = 16, height = 20)
  graphics::par(mfrow = c(1, 1), mar = c(10, 10, 10, 10))
  make_manhattan_plot(name,
                      pval_all_transition$hit_pvals,
                      fdr,
                      "continuous")

  cell_width_value <- image_width / ncol(ordered_by_p_val)

  colnames(ordered_by_p_val) <- substr(colnames(ordered_by_p_val), 1, 20)

  pheatmap::pheatmap(
    ordered_by_p_val,
    main = paste0("Edges:\n hi conf trans vs delta pheno"),
    cluster_cols = TRUE,
    na_col = "grey",
    cluster_rows = FALSE,
    show_rownames = FALSE,
    color = c("white", "black"),
    annotation_col = column_annot_ordered_by_p_val,
    annotation_row = p_trans_mat,
    annotation_colors = ann_colors,
    show_colnames = TRUE,
    fontsize = 8,
    cellwidth = cell_width_value)

  print("continuous plotting bug fixing")
  print("ann_colors")
  print(ann_colors)

  print("column_annot_ordered_by_p_val")
  print(column_annot_ordered_by_p_val)

  transparent_grey <- rgb(0, 0, 0, 0.25)
  transparent_red <- rgb(1, 0, 0, 0.25)
  # ONLY MAKE THE FOLLOWING PLOTS FOR SIGNIFICANT LOCI
  for (j in 1:nrow(pval_all_transition$hit_pvals)) {
    if (pval_all_transition$hit_pvals[j, 1] < fdr) {
      graphics::par(mfrow = c(3, 3),
                    mgp   = c(3, 1, 0),
                    oma   = c(0, 0, 4, 0),
                    mar = c(4, 4, 4, 4),
                    xpd = FALSE)
      plot_continuous_phenotype(tr, pheno_vector, pheno_anc_rec)
      plot_tr_w_color_edges(tr,
                            geno_reconstruction,
                            geno_confidence,
                            "grey",
                            "red",
                            paste0(row.names(pval_all_transition$hit_pvals)[j],
                                   "\n Genotype reconstruction"),
                            "recon",
                            j,
                            "Wild type",
                            "Variant")
      plot_tr_w_color_edges(tr,
                            geno_transition,
                            geno_confidence,
                            "grey",
                            "red",
                            "Genotype transition",
                            "trans",
                            j,
                            "No transition",
                            "Transition")
      mat_p_trans_mat <- as.matrix(p_trans_mat)
      hist_abs_delta_pheno_all_edges(mat_p_trans_mat,
                                     geno_confidence,
                                     tr,
                                     j)
      hist_abs_hi_conf_delta_pheno(results_all_trans,
                                   tr,
                                   j,
                                   transparent_grey,
                                   transparent_red)
      hist_raw_hi_conf_delta_pheno(geno_transition,
                                   geno_confidence,
                                   pheno_recon_ordered_by_edges,
                                   tr,
                                   j,
                                   transparent_grey,
                                   transparent_red)

      graphics::hist(log(results_all_trans$ks_statistics[[j]]),
           breaks = perm / 10,
           col = "grey",
           border = FALSE,
           main =
             paste("KS Test Statistic Null Distribution\n-ln(FDR Corrected P-value) = ",
                   formatC(-log(pval_all_transition$hit_pvals[j, 1]),format = "e", digits = 1),
                   " P-value rank = ",
                   rank(pval_all_transition$hit_pvals)[j],
                   sep = ""),
           cex.main = hist_cex_size,
           cex.lab = hist_cex_size,
           cex.axis = hist_cex_size,
           ylab = "Count",
           xlab = "ln(KS Test Statistic)",
           ylim = c(0, 0.6 * length(results_all_trans$ks_statistics[[j]])),
           xlim = c(min(log(as.numeric(results_all_trans$observed_ks_stat[j])),
                        log(results_all_trans$ks_statistics[[j]])), 0))
      graphics::abline(v =
                         log(as.numeric(results_all_trans$observed_ks_stat[j])),
                       col = rgb(1, 0, 0, 0.25),
                       lwd = 4)

      legend("topleft",
             title = "KS Test Statistics",
             legend = c("Null", "Observed"),
             col = c("grey", rgb(1, 0, 0, 0.25)),
             pch = 15,
             cex = hist_cex_size,
             bg = rgb(0, 0, 0, 0.01))

    }
  }

  grDevices::dev.off()

  trans_dir_edge_mat <- NULL
  for (i in 1:length(geno_transition)) {
    trans_dir_edge_mat <- cbind(trans_dir_edge_mat,
                                geno_transition[[i]]$trans_dir)
  }

  colnames(trans_dir_edge_mat) <- colnames(geno)

  for (c in 1:ncol(trans_dir_edge_mat)) {
    trans_dir_edge_mat[(1:ape::Nedge(tr))[geno_confidence[[c]] == 0], c] <- NA
  }

  all_tables <- all_lists <- rep(list(NULL), ncol(trans_dir_edge_mat))
  delta_pheno_table <- matrix(0, nrow = 3, ncol = 1)
  row.names(delta_pheno_table) <- c("geno_parent_0_child_1",
                                    "geno_parent_1_child_0",
                                    "geno_no_change")
  colnames(delta_pheno_table) <- c("sum(|delta_phenotype|)")
  for (i in 1:ncol(trans_dir_edge_mat)) {
    temp_table <- delta_pheno_table
    temp_table[1, 1] <-
      sum(p_trans_mat[which(trans_dir_edge_mat[, i] == 1),  1], na.rm = TRUE)
    temp_table[2, 1] <-
      sum(p_trans_mat[which(trans_dir_edge_mat[, i] == -1), 1], na.rm = TRUE)
    temp_table[3, 1] <-
      sum(p_trans_mat[which(trans_dir_edge_mat[, i] == 0),  1], na.rm = TRUE)
    all_tables[[i]] <- temp_table
  }
  names(all_tables) <- colnames(geno)
  delta_pheno_table <- all_tables
  delta_pheno_list <- rep(list(0), 3)
  names(delta_pheno_list) <- c("geno_parent_0_child_1",
                               "geno_parent_1_child_0",
                               "geno_no_change")
  for (i in 1:ncol(trans_dir_edge_mat)) {
    temp_list <- delta_pheno_list
    temp_list[[1]] <- p_trans_mat[which(trans_dir_edge_mat[, i] == 1),  1]
    temp_list[[2]] <- p_trans_mat[which(trans_dir_edge_mat[, i] == -1), 1]
    temp_list[[3]] <- p_trans_mat[which(trans_dir_edge_mat[, i] == 0),  1]
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
} # end plot_continuous_results()

#' make_manhattan_plot
#'
#' @description Create a manhattan plot of GWAS hit p-values.
#'
#' @noRd
#' @param geno_pheno_name Character. Will be printed in title of manhattan plot.
#' @param pval_hits Vector. P-values for all loci.
#' @param fdr Number. False discorvery rate.
#' @param test_name Character. Either "continuous", "synchronous", or "phyc".
#'  Will be printed in title of manhattan plot.
#'
#' @return A manhattan plot of all of the hit values.
make_manhattan_plot <- function(geno_pheno_name,
                                pval_hits,
                                fdr,
                                test_name){
  # Check inputs ---------------------------------------------------------------
  check_is_string(geno_pheno_name)
  check_num_between_0_and_1(fdr)
  check_str_is_test_name(test_name)

  # Function -------------------------------------------------------------------
  # Create negative log p-values with arbitrary locus numbers
  manhattan_cex <- 2
  neg_log_p_value <- data.frame(-log(pval_hits))
  neg_log_p_with_num <- cbind(1:nrow(neg_log_p_value), neg_log_p_value)
  colnames(neg_log_p_with_num)[1] <- "Locus Significance"
  sig_temp <- subset(neg_log_p_with_num, neg_log_p_with_num[, 2] > -log(fdr))
  ymax <- max(-log(0.01), neg_log_p_with_num[, 2, drop = TRUE])
  with(neg_log_p_with_num,
       graphics::plot(x = neg_log_p_with_num[, 1],
            y = jitter(neg_log_p_with_num[, 2, drop = TRUE]),
            type = "p",
            cex.main = manhattan_cex,
            cex.lab = manhattan_cex,
            cex.axis = manhattan_cex,
            main = paste(test_name, geno_pheno_name, sep = " "),
            col = grDevices::rgb(0, 0, 0, 0.3),
            pch = 19,
            xaxt = 'n',
            xlab = "Genetic loci",
            ylim = c(0, ymax),
            ylab = "-ln(FDR Corrected P-value)"))

  graphics::abline(h = -log(fdr),
                   col = "red")
  if (nrow(sig_temp) > 0) {
    graphics::text(x = sig_temp[, 1],
                   y = sig_temp[, 2],
                   labels = row.names(sig_temp),
                   pos = 1,
                   cex = 1)
  }
  legend("topright",
         bty = "n",
         legend = "Significance Threshold",
         col = "red",
         pch = "-",
         cex = manhattan_cex)
} #end make_manhattan_plot()

#' plot_tr_w_color_edges
#'
#' Plot a phylogenetic tree with certain edges highlighted.
#'
#' @noRd
#' @param tr Phylo.
#' @param edges_to_highlight List of vectors.
#' @param geno_confidence List of vectors.
#' @param edge_color_na. Character. Color.
#' @param edge_color_bright Character. Color.
#' @param title Character. Plot title.
#' @param trans_or_recon Character. Either "recon" or "trans."
#' @param index Number.
#' @param legend_baseline Character. Legend name for black lines.
#' @param legend_highlighted Character. Legend name for highlighted lines.
#'
#' @return Plot of a tree with certain edges colored.
#'
plot_tr_w_color_edges <- function(tr,
                                  edges_to_highlight,
                                  geno_confidence,
                                  edge_color_na,
                                  edge_color_bright,
                                  title,
                                  trans_or_recon,
                                  index,
                                  legend_baseline,
                                  legend_highlighted){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_is_string(edge_color_na)
  check_is_string(edge_color_bright)
  if (edge_color_na == edge_color_bright) {
    stop("These tree edges need to be different colors")
  }
  check_is_string(legend_baseline)
  check_is_string(legend_highlighted)
  check_is_number(index)
  if (index > length(edges_to_highlight)) {
    stop("Index must select an item from edges_to_highlight")
  }
  check_equal(length(geno_confidence), length(edges_to_highlight))
  check_equal(length(geno_confidence[[index]]), ape::Nedge(tr))
  check_is_string(trans_or_recon)
  if (!trans_or_recon %in% c("recon", "trans")) {
    stop("String must be trans or recon")
  }
  check_is_string(title)

  # Function -------------------------------------------------------------------
  tree_legend_cex <- 1.3
  edge_color_baseline <- "black"
  edge_color <- rep(edge_color_baseline, ape::Nedge(tr))
  if (trans_or_recon == "recon") {
    check_equal(length(edges_to_highlight[[index]]), ape::Nedge(tr))
    edge_color[edges_to_highlight[[index]] == 1] <- edge_color_bright
  } else if (trans_or_recon == "trans") {
    check_equal(length(edges_to_highlight[[index]]$transition), ape::Nedge(tr))
    edge_color[edges_to_highlight[[index]]$transition == 1] <- edge_color_bright
  }
  edge_color[geno_confidence[[index]] == 0] <- edge_color_na # grey out long
  # edges and low ML bootstrap support
  graphics::par(mar = c(4, 4, 4, 4))
  graphics::plot(tr,
                 font = 1,
                 edge.color = edge_color,
                 main = title,
                 use.edge.length = FALSE,
                 label.offset = 0.25,
                 adj = 0,
                 cex = tree_legend_cex)
  graphics::legend("topleft",
                   bty = "n",
                   legend = c(legend_baseline,
                              legend_highlighted,
                              "Low confidence"),
                   col = c(edge_color_baseline,
                           edge_color_bright,
                           edge_color_na),
                   lty = 1,
                   ncol = 1,
                   lwd = 1,
                   cex = tree_legend_cex)
} # end plot_tr_w_color_edges()

#' make_ann_colors
#'
#' @description Create object to annotate columns in the significant hit
#'  results.
#'
#' @noRd
#' @param geno_matrix Genotype matrix to plot in heatmap.
#'
#' @return Annotation color object for heatmap.
make_ann_colors <- function(geno_matrix){
  # Check input ----------------------------------------------------------------
  check_dimensions(geno_matrix, min_rows = 1, min_cols = 1)

  # Function -------------------------------------------------------------------
  ones <- sum(geno_matrix == 1, na.rm = TRUE) > 0
  zeros <- sum(geno_matrix == 0, na.rm = TRUE) > 0
  nas <- sum(is.na(geno_matrix)) > 0

  if (ones + zeros + nas == 3) {
    ann_colors <- list(`Phenotype Presence/Absence` = c(na = "grey",
                                         absent = "white",
                                         present = "red"))
  } else if (ones == 1 && zeros == 1 && nas == 0) {
    ann_colors <- list(`Phenotype Presence/Absence` = c(absent = "white", present = "red"))
  } else if (ones == 1 && zeros == 0 && nas == 1) {
    ann_colors <- list(`Phenotype Presence/Absence` = c(na = "grey", present = "red"))
  } else if (ones == 0 && zeros == 1 && nas == 1) {
    ann_colors <- list(`Phenotype Presence/Absence` = c(na = "grey", absent = "white"))
  } else if (ones == 0 && zeros == 0 && nas == 1) {
    ann_colors <- list(`Phenotype Presence/Absence` = c(na = "grey"))
  } else if (ones == 0 && zeros == 1 && nas == 0) {
    ann_colors <- list(`Phenotype Presence/Absence` = c(absent = "white"))
  } else if (ones == 1 && zeros == 0 && nas == 0) {
    ann_colors <- list(`Phenotype Presence/Absence` = c(present = "red"))
  } else {
    stop("No ones, zeroes, or NAs present in g_mat")
  }
  # Return output --------------------------------------------------------------
  return(ann_colors)
}

#' Plot PhyC results
#'
#' @noRd
#' @param tr Phylo.
#' @param dir Directory where to save plots.
#' @param name Prefix in plot file name.
#' @param fdr Numeric. False discovery rate. Between 0 and 1.
#' @param num_perm Numeric. Number of permutations.
#' @param recon_hit_vals Dataframe. Nrows = number of genotypes. Ncol = 1.
#'   Corrected p-values for each genotype tested.
#' @param p_recon_edges Vector. Length = Nedge(tree). Reconstruction of
#'   phenotype.
#' @param recon_perm_obs_results List of many results.
#'   $hit_pvals. Character. P-val for each genotype. Length = number of tested
#'   genotypes. $permuted_count. List of vectors. 1 vector for each tested
#'   genotype. Length of each vector = number of permuations. $observed_overlap.
#'   Integer vector. Length = number of tested genotypes.
#' @param tr_and_pheno_hi_conf Vector of logicals. TRUE = high confidence.
#'   FALSE = low confidence. Length = Nedge(tree).
#' @param geno_confidence List of vectors. Length of individual vector =
#'   Nedge(tree). Genotype high confidence edges. Either 1 (high confidence) or
#'   0 (low confidence).
#' @param g_trans_edges List of vectors. Length of individual vector =
#'   Nedge(tree). Genotype transition edges. Either 1 (transition) or 0 (no
#'   transition).
#' @param p_trans_edges Vector. Length = Nedge(tree). Transitions marked as 1,
#'   not transition marked as 0.
#' @param snp_in_gene Either NULL or table of integers where each entry
#'   corresponds to one genotype.
#'
#' @return  Plots printed into one pdf.
plot_phyc_results <- function(tr,
                               dir,
                               name,
                               fdr,
                               num_perm,
                               recon_hit_vals,
                               p_recon_edges,
                               recon_perm_obs_results,
                               tr_and_pheno_hi_conf,
                               geno_confidence,
                               g_trans_edges,
                               p_trans_edges,
                               snp_in_gene,
                               prefix,
                               grouped_logical){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_num_between_0_and_1(fdr)
  check_if_permutation_num_valid(num_perm)
  if (ncol(recon_hit_vals) != 1 |
      nrow(recon_hit_vals) != length(geno_confidence)) {
    stop("Dimensions of hit p-values dataframe are incorrect.")
  } # Don't change to check_dimensions because input is dataframe, not matrix.
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
  check_class(grouped_logical, "logical")
  check_str_is_test_name(prefix)

  # Function -------------------------------------------------------------------
  image_width <- 250
  hist_cex_size <- 1.3
  if (grouped_logical) {
    fname <- paste0(dir, "/hogwash_", prefix, "_grouped_", name, ".pdf")
  } else {
    fname <- paste0(dir, "/hogwash_", prefix, "_", name, ".pdf")
  }

  grDevices::pdf(fname, width = 16, height = 20)

  graphics::par(mfrow = c(1, 1), mar = c(10, 10, 10, 10))
  make_manhattan_plot(name, recon_hit_vals, fdr, prefix)

  g_trans_mat <- matrix(0, nrow = ape::Nedge(tr), ncol = length(g_trans_edges))

  for (i in 1:length(g_trans_edges)) {
    g_trans_mat[, i] <- g_trans_edges[[i]]
    g_trans_mat[geno_confidence[[i]] == 0, i] <- NA
  }

  p_recon_edges[tr_and_pheno_hi_conf == 0] <- -1
  p_mat <- matrix(p_recon_edges, nrow = length(p_recon_edges), ncol = 1)
  colnames(p_mat) <- "Phenotype Presence/Absence"
  phenotype_annotation <- as.data.frame(p_mat)
  row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)

  temp_g_trans_mat <- cbind(phenotype_annotation, g_trans_mat)
  g_trans_mat <- temp_g_trans_mat[order(temp_g_trans_mat[, 1],
                                        na.last = FALSE,
                                        decreasing = FALSE ),
                                  2:ncol(temp_g_trans_mat),
                                  drop = FALSE]
  colnames(g_trans_mat) <- row.names(recon_hit_vals)

  cell_width_value <- image_width / ncol(g_trans_mat)

  significant_loci <-
    data.frame("Locus Significance" = rep("Not Significant", ncol(g_trans_mat)),
                                 stringsAsFactors = FALSE)
  row.names(significant_loci) <- row.names(recon_hit_vals)
  log_p_value <- data.frame(-log(recon_hit_vals))
  significant_loci[log_p_value > -log(fdr)] <- "Significant"

  if (!is.null(snp_in_gene)) {
    snp_in_gene <- as.data.frame(snp_in_gene, row.names = 1)
    colnames(snp_in_gene) <- "Variants in Group"
    snp_in_gene <-
      snp_in_gene[row.names(snp_in_gene) %in% row.names(log_p_value), ,
                  drop = FALSE]
    column_annot <- cbind(significant_loci, log_p_value, snp_in_gene)
  } else {
    column_annot <- cbind(significant_loci, log_p_value)
  }

  ordered_by_p_val <-
    g_trans_mat[,
                 match(row.names(log_p_value)[order(log_p_value[, 1])],
                       colnames(g_trans_mat)),
                 drop = FALSE]
  column_annot_ordered_by_p_val <-
    column_annot[match(row.names(log_p_value)[order(log_p_value[, 1])],
                        row.names(column_annot)), , drop = FALSE]

  if (ncol(column_annot_ordered_by_p_val) == 2) {
    colnames(column_annot_ordered_by_p_val) <- c("Locus Significance",
                                                 "-ln(FDR Corrected P-value)")
  } else {
    colnames(column_annot_ordered_by_p_val) <- c("Locus Significance",
                                                 "-ln(FDR Corrected P-value)",
                                                 "Variants in Group")
  }


  if (length(unique(phenotype_annotation[, 1])) == 3) {
    pheno_presence_col <-
      c( `Low Confidence` = "grey", Absence = "white", Presence = "red")
  } else if (length(unique(phenotype_annotation[, 1])) == 2) {
    if (sum(unique(phenotype_annotation[, 1]) %in% c(-1, 0)) == 2) {
      pheno_presence_col <- c(`Low Confidence` = "grey", Absence = "white")
    } else if (sum(unique(phenotype_annotation[, 1]) %in% c(-1, 1)) == 2) {
      pheno_presence_col <- c( `Low Confidence` = "grey", Presence = "red")
    } else if (sum(unique(phenotype_annotation[, 1]) %in% c(1, 0)) == 2) {
      pheno_presence_col <- c(Absence = "white", Presence = "red")
    }
  } else if (length(unique(phenotype_annotation[, 1])) == 1) {
    pheno_presence_col <- c(Presence = "red")
    if (unique(phenotype_annotation[, 1]) == -1) {
      pheno_presence_col <- c(`Low Confidence` = "grey")
    } else if (unique(phenotype_annotation[, 1]) == 0) {
      pheno_presence_col <- c(Absence = "white")
    }
  }

  ann_colors <-
    list(`Locus Significance` = c(`Not Significant` = "white",
                                  Significant = "blue"),
         `Phenotype Presence/Absence` = pheno_presence_col)


  print("p-values")
  print(log_p_value)

  print("PHYC plotting bug fixing")
  print("ann_colors")
  print(ann_colors)
  print("column_annot_ordered_by_p_val")
  print(column_annot_ordered_by_p_val)

  phenotype_annotation[phenotype_annotation == -1] <- "Low Confidence"
  phenotype_annotation[phenotype_annotation == 0] <- "Absence"
  phenotype_annotation[phenotype_annotation == 1] <- "Presence"

  print("phenotype_annotation")
  print(phenotype_annotation)

  print("matrix")
  print("ordered_by_p_val")
  print(ordered_by_p_val)

  plotting_logical <- check_if_g_mat_can_be_plotted(ordered_by_p_val)

  print("grouped_logical")
  print(grouped_logical)
  print("plotting_logical")
  print(plotting_logical)

  if (plotting_logical) {
    colnames(ordered_by_p_val) <- substr(colnames(ordered_by_p_val), 1, 20)

    pheatmap::pheatmap(
      ordered_by_p_val,
      main = "PhyC Convergence Summary",
      cluster_cols  = TRUE,
      na_col = "grey",
      cluster_rows  = FALSE,
      show_rownames = FALSE,
      color = c("white", "black"),
      annotation_col = column_annot_ordered_by_p_val,
      annotation_row = phenotype_annotation,
      annotation_colors = ann_colors,
      show_colnames = TRUE,
      fontsize = 8,
      cellwidth = cell_width_value)
  }

  # Loop through significant hits:
  pheno_as_list <- list(p_recon_edges)
  pheno_conf_as_list <- list(tr_and_pheno_hi_conf)
  for (j in 1:nrow(recon_hit_vals)) {
    if (recon_hit_vals[j, 1] < fdr) {
      graphics::par(mfrow = c(3, 2),
                    mgp = c(3, 1, 0),
                    oma = c(0, 0, 4, 0),
                    mar = c(4, 4, 4, 4))
      # Phenotype
      plot_tr_w_color_edges(tr,
                            pheno_as_list,
                            pheno_conf_as_list,
                            "grey",
                            "red",
                            paste0("\n Phenotype Reconstruction"),
                            "recon",
                            1,
                            "Wild type",
                            "Variant")
      blank_plot()

      # Genotype
      plot_tr_w_color_edges(tr,
                            g_trans_edges,
                            geno_confidence,
                            "grey",
                            "red",
                            paste0(row.names(recon_hit_vals)[j],
                                   "\nGenotype Transitions"),
                            "recon",
                            j,
                            "No transition",
                            "Transition")
      blank_plot()

      # Permutation test
      max_x <- 1.2 * max(recon_perm_obs_results$permuted_count[[j]],
                   recon_perm_obs_results$observed_overlap[j])
      max_y <- 0.85 * length(recon_perm_obs_results$permuted_count[[j]])
      graphics::hist(recon_perm_obs_results$permuted_count[[j]],
                     breaks = num_perm / 10,
                     xlim = c(0, max_x),
                     ylim = c(0, max_y),
                     col = "grey",
                     border = FALSE,
                     ylab = "Count",
                     xlab =
                       "Genotype Transition & Phenotype Presence Co-occurrence",
                     main = paste0(
                       "Co-occurence Null Distribution\n -ln(FDR Corrected P-value) = ",
                       formatC(-log(recon_hit_vals[j, 1]), format = "e", digits = 1),
                       " P-value rank = ",
                       rank(recon_hit_vals[j, ]),
                       sep = ""))
      graphics::abline(v = recon_perm_obs_results$observed_overlap[j],
                       col = rgb(1, 0, 0, 0.25),
                       lwd = 4)
      graphics::legend("topleft",
                       title = "Co-occurence",
                       legend = c("Null", "Observed"),
                       col = c( "grey", rgb(1, 0, 0, 0.25)),
                       pch = 15,
                       cex = hist_cex_size)
    }
  }

  grDevices::dev.off()
} # end plot_phyc_results()

#' plot_synchronous_results
#'
#' @description Plot the synchronous test results.
#'
#' @noRd
#' @param tr Phylo.
#' @param dir Directory where to save plots.
#' @param name Prefix in plot file name.
#' @param fdr Numeric. False discovery rate. Between 0 and 1.
#' @param num_perm Numeric. Number of permutations.
#' @param trans_hit_vals Dataframe. Nrows = number of genotypes. Ncol = 1.
#'   Corrected p-values for each genotype tested.
#' @param trans_perm_obs_results List of many results.  $hit_pvals. Character.
#'   P-val for each genotype. Length = number of tested genotypes.
#'   $permuted_count. List of vectors. 1 vector for each tested genotype. Length
#'   of each vector = number of permuations. $observed_overlap. Integer vector.
#'   Length = number of tested genotypes.
#' @param tr_and_pheno_hi_conf Vector of logicals. TRUE = high confidence. FALSE
#'   = low confidence. Length = Nedge(tree).
#' @param geno_confidence List of vectors. Length of individual vector =
#'   Nedge(tree). Genotype high confidence edges. Either 1 (high confidence) or
#'   0 (low confidence).
#' @param g_trans_edges List of vectors. Length of individual vector =
#'   Nedge(tree). Genotype transition edges. Either 1 (transition) or 0 (no
#'   transition).
#' @param p_trans_edges Vector. Length = Nedge(tree). Transitions marked as 1,
#'   not transition marked as 0.
#' @param snp_in_gene Either NULL or Table of integers where each entry
#'   corresponds to one genotype.
#'
#' @return  Plots printed into one pdf.
plot_synchronous_results  <- function(tr,
                                 dir,
                                 name,
                                 fdr,
                                 num_perm,
                                 trans_hit_vals,
                                 trans_perm_obs_results,
                                 tr_and_pheno_hi_conf,
                                 geno_confidence,
                                 g_trans_edges,
                                 p_trans_edges,
                                 snp_in_gene,
                                 prefix,
                                 grouped_logical){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_num_between_0_and_1(fdr)
  check_if_permutation_num_valid(num_perm)
  if (ncol(trans_hit_vals) != 1 |
      nrow(trans_hit_vals) != length(geno_confidence)) {
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
  check_class(trans_perm_obs_results$observed_overlap, "integer")
  check_str_is_test_name(prefix)
  check_class(grouped_logical, "logical")

  # Function -------------------------------------------------------------------
  image_width <- 250
  hist_cex_size <- 1.3
  if (grouped_logical) {
    fname <- paste0(dir, "/hogwash_", prefix, "_grouped_", name, ".pdf")
  } else {
    fname <- paste0(dir, "/hogwash_", prefix, "_", name, ".pdf")
  }
  grDevices::pdf(fname, width = 16, height = 20)

  graphics::par(mfrow = c(1, 1), mar = c(10, 10, 10, 10))
  make_manhattan_plot(name, trans_hit_vals, fdr, prefix)

  # heatmaps
  g_trans_mat <- matrix(0, nrow = ape::Nedge(tr), ncol = length(g_trans_edges))

  for (i in 1:length(g_trans_edges)) {
    g_trans_mat[, i] <- g_trans_edges[[i]]
    g_trans_mat[geno_confidence[[i]] == 0, i] <- NA
  }

  # TODO -1 should be NA but it won't work correctedly
  p_trans_edges[tr_and_pheno_hi_conf == 0] <- -1
  p_mat <- matrix(p_trans_edges, nrow = length(p_trans_edges), ncol = 1)
  colnames(p_mat) <- "Phenotype Transitions/Non"
  phenotype_annotation <- as.data.frame(p_mat)
  row.names(phenotype_annotation) <- 1:nrow(phenotype_annotation)

  temp_g_trans_mat <- cbind(phenotype_annotation, g_trans_mat)
  g_trans_mat <-
    temp_g_trans_mat[order(temp_g_trans_mat[, 1],
                           na.last = FALSE, decreasing = FALSE ),
                     2:ncol(temp_g_trans_mat),
                     drop = FALSE]
  colnames(g_trans_mat) <- row.names(trans_hit_vals)

  significant_loci <-
    data.frame("Locus Significance" = rep("Not Significant", ncol(g_trans_mat)),
               stringsAsFactors = FALSE)
  row.names(significant_loci) <- row.names(trans_hit_vals)
  log_p_value <- data.frame(-log(trans_hit_vals))
  significant_loci[log_p_value > -log(fdr)] <- "Significant"

  if (!is.null(snp_in_gene)) {
    snp_in_gene <- as.data.frame(snp_in_gene, row.names = 1)
    colnames(snp_in_gene) <- "# grouped genotypes"
    snp_in_gene <-
      snp_in_gene[row.names(snp_in_gene) %in% row.names(log_p_value), ,
                  drop = FALSE]
    column_annot <- cbind(significant_loci, log_p_value, snp_in_gene)
  } else {
    column_annot <- cbind(significant_loci, log_p_value)
  }

  ordered_by_p_val <-
    g_trans_mat[, match(row.names(log_p_value)[order(log_p_value[, 1])],
                         colnames(g_trans_mat)),
                 drop = FALSE]
  column_annot_ordered_by_p_val <-
    column_annot[match(row.names(log_p_value)[order(log_p_value[, 1])],
                       row.names(column_annot)), , drop = FALSE]
  if (ncol(column_annot_ordered_by_p_val) == 2) {
    colnames(column_annot_ordered_by_p_val) <- c("Locus Significance",
                                                 "-ln(FDR Corrected P-value)")
  } else {
    colnames(column_annot_ordered_by_p_val) <- c("Locus Significance",
                                                 "-ln(FDR Corrected P-value)",
                                                 "Variants in Group")
  }

  if (length(unique(phenotype_annotation[, 1])) == 3) {
    pheno_presence_col <- c(`Low Confidence` = "grey", `Non-Transition` = "white", Transition = "red")
  } else if (length(unique(phenotype_annotation[, 1])) == 2) {
    if (sum(unique(phenotype_annotation[, 1]) %in% c(-1, 0)) == 2) {
      pheno_presence_col <- c(`Low Confidence` = "grey", `Non-Transition` = "white")
    } else if (sum(unique(phenotype_annotation[, 1]) %in% c(-1, 1)) == 2) {
      pheno_presence_col <- c(`Low Confidence` = "grey", Transition = "red")
    } else if (sum(unique(phenotype_annotation[, 1]) %in% c(1, 0)) == 2) {
      pheno_presence_col <- c(`Non-Transition` = "white", Transition = "red")
    }
  } else if (length(unique(phenotype_annotation[, 1])) == 1) {
    pheno_presence_col <- c(Transition = "red")
    if (unique(phenotype_annotation[, 1]) == -1) {
      pheno_presence_col <- c(`Low Confidence` = "grey")
    } else if (unique(phenotype_annotation[, 1]) == 0) {
      pheno_presence_col <- c(`Non-Transition` = "white")
    }
  }

  if (length(unique(column_annot_ordered_by_p_val[, 1])) == 2) {
    locus_col <- c(`Not Significant` = "white", Significant = "blue")
  } else if (length(unique(column_annot_ordered_by_p_val[, 1])) == 1) {
    locus_col <- c(Significant = "blue")
    if (unique(column_annot_ordered_by_p_val[, 1]) == "Not Significant") {
      locus_col <- c(`Not Significant` = "white")
    }
  }

  phenotype_annotation[phenotype_annotation == -1] <- "Low Confidence"
  phenotype_annotation[phenotype_annotation == 0] <- "Non-Transition"
  phenotype_annotation[phenotype_annotation == 1] <- "Transition"

  ann_colors <- list(`Locus Significance` = locus_col,
                     `Phenotype Transitions/Non` = pheno_presence_col)
  can_be_plotted <- check_if_g_mat_can_be_plotted(ordered_by_p_val)
  if (can_be_plotted) {
    cell_width_value <- image_width / ncol(ordered_by_p_val)
    colnames(ordered_by_p_val) <- substr(colnames(ordered_by_p_val), 1, 20)

    # Transition loci summary heat maps
    pheatmap::pheatmap( # Plot the heatmap
      ordered_by_p_val,
      main = paste0("Edges:\nGenotype transitions with phenotype transitions"),
      cluster_cols = TRUE,
      na_col = "grey",
      cluster_rows  = FALSE,
      show_rownames = FALSE,
      color = c("white", "black"),
      annotation_col = column_annot_ordered_by_p_val,
      annotation_row = phenotype_annotation,
      annotation_colors = ann_colors,
      fontsize = 8,
      show_colnames = TRUE,
      cellwidth = cell_width_value)
  }
  pheno_conf_as_list <- list(tr_and_pheno_hi_conf)
  p_trans_edges_as_list <- list(p_trans_edges)

   graphics::par(mfrow = c(1, 1),
                mgp   = c(3, 1, 0),
                oma   = c(0, 0, 4, 0),
                mar = c(4, 4, 4, 4),
                xpd = FALSE)
  plot_tr_w_color_edges(tr,
                        p_trans_edges_as_list,
                        pheno_conf_as_list,
                        "grey",
                        "red",
                        paste0("\n Phenotype"),
                        "recon",
                        1,
                        "No transition",
                        "Transition")


  # Loop through significant hits:
  for (j in 1:nrow(trans_hit_vals)) {
    if (trans_hit_vals[j, 1] < fdr) {
      graphics::par(mfrow = c(2, 2),
                    mgp   = c(3, 1, 0),
                    oma   = c(0, 0, 4, 0),
                    mar = c(4, 4, 4, 4),
                    xpd = FALSE)

      # Plot genotype
      plot_tr_w_color_edges(tr,
                            g_trans_edges,
                            geno_confidence,
                            "grey",
                            "red",
                            paste0(row.names(trans_hit_vals)[j],
                                   "\n Genotype transitions"),
                            "recon",
                            j,
                            "No transition",
                            "Transition")

      blank_plot()

      # Plot permutation test
      max_x <- max(trans_perm_obs_results$permuted_count[[j]],
                   trans_perm_obs_results$observed_overlap[j])
      graphics::hist(trans_perm_obs_results$permuted_count[[j]],
           breaks = num_perm / 10,
           xlim = c(0, max_x),
           ylim = c(0, .85 * length(trans_perm_obs_results$permuted_count[[j]])),
           col = "grey",
           border = FALSE,
           ylab = "Count",
           xlab = "Genotype & Phenotype Transition Edge Co-occurrence",
           main = paste0("Co-occurence Null Distribution\n -ln(FDR Corrected P-value) = ",
                         formatC(-log(trans_hit_vals[j, 1]), format = "e", digits = 1),
                         " P-value rank = ",
                         rank(trans_hit_vals[j, ]),
                         sep = ""))
      graphics::abline(v = trans_perm_obs_results$observed_overlap[j],
                       col = rgb(1, 0, 0, 0.25),
                       lwd = 4)
      legend("topleft",
             title = "Co-occurence",
             legend = c("Null", "Observed"),
             col = c("grey", rgb(1, 0, 0, 0.25)),
             pch = 15,
             cex = hist_cex_size,
             bg = rgb(0, 0, 0, 0.01))
    }
  }
  grDevices::dev.off()
} # end plot_synchronous_results()


#' Draw a blank plot
blank_plot <- function(){
  plot(0,
       type = 'n',
       axes = FALSE,
       ann = FALSE)
}
