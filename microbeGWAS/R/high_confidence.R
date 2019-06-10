discretize_confidence_using_threshold <- function(confidence_vector, threshold){
  # Function description -------------------------------------------------------
  # Given a vector with values that describe confidence, binarize the vector a
  # accoriding to a cutoff value.
  #
  # Input:
  # Confidence vector. Numeric vector.
  # Threshold. Number.
  #
  # Output:
  # Confidence vector. Binary vector.
  #
  # Check inputs ---------------------------------------------------------------
  check_is_number(threshold)
  if (!is.vector(confidence_vector)){stop("input must be a numeric vector")}
  check_is_number(confidence_vector[1])

  # Function -------------------------------------------------------------------
  confidence_vector[confidence_vector  < threshold] <- 0
  confidence_vector[confidence_vector >= threshold] <- 1

  # Check and return output ----------------------------------------------------
  check_if_binary_vector(confidence_vector)
  return(confidence_vector)
} # end discretize_confidence_using_threshold()


report_num_high_confidence_trans_edge <- function(genotype_transition,
                                                  high_conf_edges, geno_names){
  # Function description -------------------------------------------------------
  # Given a genotype for which you have: a list of vectors that stores if there
  # is a genotype transition or not for each edge (genotype_transition), a list
  # of vectors that stores if that edge is high confidence or not
  # (high_conf_edges), and a character vector of the genotype names -- create an
  # object that stores the number of high confidence transition edges per
  # genotype.
  #
  # Inputs:
  # genotype_transition. List of numeric vectors. Number of lists == number of
  #                      genotypes. Length of vector == Nedge(tr)
  # high_conf_edges. Binary vector. List of numeric vectors. Number of lists ==
  #                  number of genotypes. Length of vector == Nedge(tr)
  # geno_names. Character vector. Length == ncol(genotype_matrix)
  #
  # Outputs:
  # num_high_confidence_transition_edges. Numeric vector. Count of number of
  #                                       high confidence transitions per
  #                                       genotype. Vector is named with
  #                                       genotype names.
  #
  # Check input ----------------------------------------------------------------
  if (length(genotype_transition) != length(geno_names)){
    stop("genotype_transition should have one vector for each genotype")
  }
  if (length(high_conf_edges) != length(geno_names)){
    stop("high_conf_edges should have one vector for each genotype")
  }
  if(!is.vector(genotype_transition[[1]]$transition)){stop("Input must have a vector.")}
  if(!is.vector(high_conf_edges[[1]])){stop("Input must have a vector.")}
  if(!is.vector(geno_names)){stop("Input must be a vector.")}
  if(!is.character(geno_names[1])){stop("Input must be a character vector.")}

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  num_high_confidence_transition_edges <- rep(0, length(high_conf_edges))
  for (p in 1:length(high_conf_edges)){
    num_high_confidence_transition_edges[p] <- sum(genotype_transition[[p]]$transition * high_conf_edges[[p]])
  }

  names(num_high_confidence_transition_edges) <- geno_names
  return(num_high_confidence_transition_edges)
} # end report_num_high_confidence_trans_edge

assign_high_confidence_to_transition_edges_including_parent_info  <- function(tr,
                                                                             all_confidence_by_edge,
                                                                             genotype_transition_by_edges,
                                                                             geno){
  # This function used to be called: assign_high_confidence_to_transition_edges
  # but on 2019-05-15 I figured out that this is waaay too stringent and was dropping
  # all or nearly all of my genotypes. I rewrote the function to only look at the
  # edge in question and to ignore the parent edge.
  # I don't know if this function will ever get used for anything.

  # Function description -------------------------------------------------------
  # Identify all edges for which the edge and the parent edge are both high confidence
  #
  # Inputs:
  # tr. Phylo.
  # all_confidence_by_edge. List of vectors. Each vector is binary. Length(list) == number of genomes.
  # genotype_transition_by_edges. List of vectors. Each vector is binary. Length(list) == number of genomes.
  # geno. Matrix. Binary.
  #
  # Outputs:
  # edge_and_parent_confident_and_trans_edge. List of vector. Each vector is binary. Length(list) == number of genomes.
  #
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  if (length(genotype_transition_by_edges[[1]]$transition) != Nedge(tr)){
    stop("Dimension mismatch")
  }
  check_for_root_and_bootstrap(tr)
  if (length(all_confidence_by_edge[[1]]) != Nedge(tr)){
    stop("Dimension mismatch")
  }
  # Function -------------------------------------------------------------------
  # Identify all edges for which the edge and the parent edge are both high confidence
  edge_and_parent_both_confident <- rep(list(rep(0, Nedge(tr))), ncol(geno))
  for (ge in 1:ncol(geno)){
    for (ed in 1:Nedge(tr)){
      parent_edge <- find_parent_edge(tr, ed)
      if (!is.na(parent_edge)){
        if (all_confidence_by_edge[[ge]][ed] == 1){
          if (all_confidence_by_edge[[ge]][parent_edge] == 1){
            edge_and_parent_both_confident[[ge]][ed] <- 1
          }
        }
      } else{
        # have to account for the fact that there isn't a parent edge two tree edges
        edge_and_parent_both_confident[[ge]][ed] <- all_confidence_by_edge[[ge]][ed]
      } # end if (!is.na(parent_edge))
    }
  }

  # Identify high confidence transition edges by overlapping the above and transitions
  edge_and_parent_confident_and_trans_edge <- rep(list(NULL), ncol(geno))
  for (k in 1:ncol(geno)){
    edge_and_parent_confident_and_trans_edge[[k]] <- as.numeric((edge_and_parent_both_confident[[k]] + genotype_transition_by_edges[[k]]$transition) == 2)
  }

  # Return output --------------------------------------------------------------
  return(edge_and_parent_confident_and_trans_edge)
} # end assign_high_confidence_to_transition_edges_including_parent_info()

assign_high_confidence_to_transition_edges <- function(tr,
                                                       all_confidence_by_edge,
                                                       genotype_transition_by_edges,
                                                       geno){
  # Function description -------------------------------------------------------
  # Identify all edges for which the edge is high confidence and a transition edge.
  #
  # Inputs:
  # tr. Phylo.
  # all_confidence_by_edge. List of vectors. Each vector is binary. Length(list) == number of genomes.
  # genotype_transition_by_edges. List of vectors. Each vector is binary. Length(list) == number of genomes.
  # geno. Matrix. Binary.
  #
  # Outputs:
  # edge_confident_and_trans_edge. List of vector. Each vector is binary. Length(list) == number of genomes.
  #
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(geno)
  if (length(genotype_transition_by_edges[[1]]$transition) != Nedge(tr)){
    stop("Dimension mismatch")
  }
  check_for_root_and_bootstrap(tr)
  if (length(all_confidence_by_edge[[1]]) != Nedge(tr)){
    stop("Dimension mismatch")
  }
  # Function -------------------------------------------------------------------
  edge_confident_and_trans_edge <- rep(list(NULL), ncol(geno))
  for (k in 1:ncol(geno)){
    edge_confident_and_trans_edge[[k]] <- as.numeric((all_confidence_by_edge[[k]] + genotype_transition_by_edges[[k]]$transition) == 2)
  }

  # Return output --------------------------------------------------------------
  return(edge_confident_and_trans_edge)
} # end assign_high_confidence_to_transition_edges()


prepare_high_confidence_objects <- function(genotype_transition, tr,
                                            pheno_tip_node_recon_conf,
                                            boot_threshold, geno, geno_conf_edge,
                                            geno_recon_edge, snps_in_each_gene){

  # Function description -------------------------------------------------------
  # Identify high confidence edges (considering: tree bootstrap values,
  # phenotype reconstruction, tree edge lengths, and ancestral reconstruction of
  # genotype).
  #
  # Inputs:
  # genotype_transition. List of lists.
  #                      Number of lists = number of genotypes.
  #                      Each list is made of a $transition and $trans_dir list.
  #                      Length(transition) == Length(trans_dir) == Nedge(tree)
  # tr. Phylo.
  # pheno_tip_node_recon_conf. List of confidence values. Binary.
  #                            Length(list) == Ntip() + Nedge()
  # boot_threshold. Numeric. Between 0 and 1.
  # geno. Matrix. Binary. Nrow = Ntip(tree). Ncol = Number of genotypes.
  # geno_conf_edge. List of lists.
  #                 Binary.
  #                 Number of lists = number of genotypes.
  #                 Length(each individual list) == Nedge(Tree)
  # geno_recon_edge. List of lists.
  #                  Binary.
  #                  Number of lists = number of genotypes.
  #                  Length(each individual list) == Nedge(Tree)
  # TODO Fill in descript and checks for snps_in_each_gene!
  # snps_in_each_gene. Either NULL or ?????.
  #
  # TODO Add descripts for each output.
  # Outputs:
  # list(c("dropped_genotypes" = dropped_genotypes,
  # "high_confidence_trasition_edges" = only_high_conf_geno_trans,
  # "genotype" = geno,
  # "snps_per_gene" = snps_in_each_gene,
  # "genotype_transition" = genotype_transition,
  # "geno_recon_edge" = geno_recon_edge,
  # "high_conf_ordered_by_edges" = high_conf_ordered_by_edges,
  # "num_high_confidence_trasition_edges" = num_high_confidence_trasition_edges))
  #
  # Check input ----------------------------------------------------------------
  if (length(genotype_transition) != ncol(geno)){
    stop("Need 1 genotype transition for each genotype")
  }
  if (length(genotype_transition[[1]]$transition) != Nedge(tr)){
    stop("Need genotype transition information for each tree edge.")
  }
  check_for_root_and_bootstrap(tr)
  if (length(pheno_tip_node_recon_conf) != c(Ntip(tr) + Nnode(tr))){
    stop("Phenotype confidence list should correspond to each tip and node of tree.")
  }
  check_num_between_0_and_1(boot_threshold)
  check_dimensions(geno, exact_rows = Ntip(tr), min_rows = Ntip(tr), exact_cols = NULL, min_cols = 1)
  if (length(geno_conf_edge) != ncol(geno)){
    stop("Need 1 genotype transition for each genotype")
  }
  if (length(geno_conf_edge[[1]]) != Nedge(tr)){
    stop("Need genotype transition information for each tree edge.")
  }
  if (length(geno_recon_edge) != ncol(geno)){
    stop("Need 1 genotype transition for each genotype")
  }
  if (length(geno_recon_edge[[1]]) != Nedge(tr)){
    stop("Need genotype transition information for each tree edge.")
  }
  check_if_binary_vector(geno_conf_edge[[1]])
  check_if_binary_vector(geno_recon_edge[[1]])

  # Function -------------------------------------------------------------------
  pheno_conf_ordered_by_edges <- reorder_tips_and_nodes_to_edges(pheno_tip_node_recon_conf, tr)
  tree_conf                   <- get_bootstrap_confidence(tr, boot_threshold)
  tree_conf_ordered_by_edges  <- reorder_tips_and_nodes_to_edges(tree_conf, tr)
  short_edges                 <- identify_short_edges(tr)


  high_confidence_edges <- pheno_conf_ordered_by_edges + tree_conf_ordered_by_edges + short_edges == 3
  if (length(high_confidence_edges) != Nedge(tr)){stop("Confidence should correspond to each tree edge")}
  all_high_confidence_edges <- rep(list(0), ncol(geno))

  # ADD IN GENO RECONSTRUCTION CONFIDENCE
  for (k in 1:ncol(geno)){
    all_high_confidence_edges[[k]] <- as.numeric(geno_conf_edge[[k]] + high_confidence_edges == 2)
  }
  only_high_conf_geno_trans <- assign_high_confidence_to_transition_edges(tr, all_high_confidence_edges, genotype_transition, geno) # here
  for (i in 1:ncol(geno)){
    genotype_transition[[i]]$transition <- only_high_conf_geno_trans[[i]]
    genotype_transition[[i]]$trans_dir <- only_high_conf_geno_trans[[i]] * genotype_transition[[i]]$trans_dir
  }
  num_high_confidence_trasition_edges <- report_num_high_confidence_trans_edge(genotype_transition, all_high_confidence_edges, colnames(geno))

  # KEEP ONLY genoS WITH AT LEAST TWO HIGH CONFIDENCE TRANSITION EDGES ----#
  geno_to_keep                <- keep_at_least_two_high_conf_trans_edges(genotype_transition, all_high_confidence_edges)
  genotype_transition         <- genotype_transition[geno_to_keep]
  geno_recon_edge             <- geno_recon_edge[geno_to_keep]
  high_conf_ordered_by_edges  <- all_high_confidence_edges[geno_to_keep]
  dropped_genotypes           <- get_dropped_genotypes(geno, geno_to_keep)
  geno                        <- geno[ , geno_to_keep, drop = FALSE]
  snps_in_each_gene           <- snps_in_each_gene[names(snps_in_each_gene) %in% colnames(geno)]

  # Return output --------------------------------------------------------------
  results = list("dropped_genotypes" = dropped_genotypes,
                  "high_confidence_trasition_edges" = only_high_conf_geno_trans,
                  "genotype" = geno,
                  "snps_per_gene" = snps_in_each_gene,
                  "genotype_transition" = genotype_transition,
                  "geno_recon_edge" = geno_recon_edge,
                  "high_conf_ordered_by_edges" = high_conf_ordered_by_edges,
                  "num_high_confidence_trasition_edges" = num_high_confidence_trasition_edges,
                 "tr_and_pheno_hi_conf" = high_confidence_edges)
  return(results)

} # end prepare_high_confidence_objects()
