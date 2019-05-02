read_in_tsv_matrix <- function(mat){
  # Function description -------------------------------------------------------
  # Read in the standardized GWAS matrix format: rows correspond to samples, columns correspond to genotypes/phenotypes.
  # Data are tab-separated.
  #
  # Inputs:
  # mat. Character. Path to matrix file.
  #
  # Output:
  # temp. Matrix.
  #
  # Check inputs ---------------------------------------------------------------
  check_is_string(mat)
  check_file_exists(mat)

  # Function -------------------------------------------------------------------
  temp <- read.table(mat,
                     sep = "\t",
                     row.names = 1,
                     header = TRUE,
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
  temp <- as.matrix(temp)

  # Check and return output ----------------------------------------------------
  if (class(temp) != "matrix"){stop("matrix is incorrectly formatted")}
  return(temp)
} # end read_in_tsv_matrix()


read_in_arguments <- function(args){
  # TODO: Update description & reduce if statements down to functions
  # Function description ------------------------------------------------------#
  # Read in the commandline arguments.
  #
  # Inputs:
  # args. Command line inputs.
  #
  # Output:
  # test
  # tree
  # phenotype
  # genotype
  # output_name
  # output_dir
  # alpha
  # discrete_or_continuous
  # annotation
  #
  # Check inputs, function, & check and return outputs -------------------------
  if (length(args) == 2){
    if (args[1] == "test"){
      test <- TRUE
    } else {
      stop("The first argument for a test run should be \"test.\"")
    }
    test_data              <- create_test_data()
    phenotype              <- test_data$phenotype
    tree                   <- test_data$tree
    genotype               <- test_data$genotype
    output_name            <- "test_data"
    output_dir             <- args[2]
    perm                   <- 1000
    alpha                  <- 0.05
    discrete_or_continuous <- check_input_format(phenotype, tree, genotype, output_name, output_dir, perm, alpha)
    annotation             <- NULL
    results <- list("test" = test, "tree" = tree, "phenotype" = phenotype,
                    "genotype" = genotype, "output_name" = output_name,
                    "output_dir" = output_dir, "perm" = perm, "alpha" = alpha,
                    "discrete_or_continuous" = discrete_or_continuous,
                    "annotation" = annotation)
    return(results)
  } else if (length(args) == 9){
    test                   <- FALSE
    phenotype              <- read_in_tsv_matrix(args[1])
    tree                   <- read.tree(args[2])
    if (!is.rooted(tree)){
      tree <- midpoint(tree)
    }
    genotype               <- read_in_tsv_matrix(args[3])
    output_name            <- args[4] # Ex: log_toxin_snp_stop
    output_dir             <- args[5] # Directory in which all output files will be saved
    perm                   <- as.numeric(args[6]) #typically 10,000
    alpha                  <- as.numeric(args[7])
    bootstrap_cutoff       <- as.numeric(args[8]) # typically 0.70
    annotation             <- read_in_tsv_matrix(args[9])
    discrete_or_continuous <- check_input_format(phenotype, tree, genotype, output_name, output_dir, perm, alpha, annotation)
    results <- list("test" = test, "tree" = tree, "phenotype" = phenotype,
                    "genotype" = genotype, "output_name" = output_name,
                    "output_dir" = output_dir, "perm" = perm, "alpha" = alpha,
                    "discrete_or_continuous" = discrete_or_continuous,
                    "annotation" = annotation, "bootstrap_cutoff" = bootstrap_cutoff)
    return(results)
  } else if (length(args) == 8){
    test                   <- FALSE
    phenotype              <- read_in_tsv_matrix(args[1])
    tree                   <- read.tree(args[2])
    # added is.rooted if statement on 2018-09-25 to deal with odd midpoint rooting issue for PSM dataset.
    if (!is.rooted(tree)){
      tree <- midpoint(tree)
    }
    # End of section added on 2018-09-25
    genotype               <- read_in_tsv_matrix(args[3])
    output_name            <- args[4] # Ex: log_toxin_snp_stop
    output_dir             <- args[5] # Directory in which all output files will be saved
    perm                   <- as.numeric(args[6]) #typically 10,000
    alpha                  <- as.numeric(args[7])
    bootstrap_cutoff       <- as.numeric(args[8])
    discrete_or_continuous <- check_input_format(phenotype, tree, genotype, output_name, output_dir, perm, alpha, NULL)
    results <- list("test" = test, "tree" = tree, "phenotype" = phenotype,
                    "genotype" = genotype, "output_name" = output_name,
                    "output_dir" = output_dir, "perm" = perm, "alpha" = alpha,
                    "discrete_or_continuous" = discrete_or_continuous,
                    "annotation" = NULL, "bootstrap_cutoff" = bootstrap_cutoff)
    return(results)
  } else {
    stop("2, 8 or 9 inputs required. \n
         Either: 1. test 2. output directory or \n
         1. Phenotype 2. Tree 3. Genotype 4. Output name 5. Output directory \n
         6. Number of permutations 7. Alpha 8. Bootstrap confidence threshold and optional 9. Annotation")
  }
} # end read_in_arguments()


check_input_format <- function(pheno, tr, geno, name, dir, perm, a, annot){
  # Function description -------------------------------------------------------
  # Check that all of the inputs into phyC are in the valid format.
  #
  # Input:
  # pheno. Matrix. Phenotype.
  # tr.    Phylo. Tree.
  # geno.  Matrix. Genotype.
  # name.  Character. Output name.
  # dir.   Character. Output path.
  # perm.  Number. Times to shuffle the data on the tree to create a null distribution for the permutation test.
  # a.     Number. Alpha.
  # annot. Matrix. Annotation for heatmaps.

  # Output:
  # discrete_or_continuous. Character. Either "discrete" or "continuous". Describes the input phenotype.
  #
  # Check input ----------------------------------------------------------------
  check_dimensions(geno, Ntip(tr), 2, NULL, 2) # Genotype matrix should have same rows as tr$tip.label, at least 2 genotypes in the columns
  check_dimensions(pheno, Ntip(tr), 2, 1, 1) # Phnoetype matrix should have same rows as tr$tip.label and exactly 1 column
  check_rownames(geno, tr) # Genotype rownames should be in the same order as the tr$tip.label
  check_rownames(pheno, tr) # Phenotype rownames should be in the same order as the tr$tip.label
  check_for_NA_and_inf(geno)
  check_for_NA_and_inf(pheno)
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  check_if_binary_matrix(geno)
  check_is_string(name)
  check_if_dir_exists(dir)
  check_if_permutation_num_valid(perm)
  check_if_alpha_valid(a)
  if (!is.null(annot)){
    check_dimensions(annot, Ntip(tr), 2, 2, 2)
  }

  # Function -------------------------------------------------------------------
  discrete_or_continuous <- assign_pheno_type(pheno)

  # Check and return output ----------------------------------------------------
  check_is_string(discrete_or_continuous)
  return(discrete_or_continuous)
} # end check_input_format()

