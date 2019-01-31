# Katie Saund
# 2018-08-29
# Library of functiosn to plot GWAS hits. 

library(ape)
library(ComplexHeatmap)

print(sessionInfo())

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
  # the first section here is to improve speed of loading my larger genotype matrices: 
  sub_sample <- read.table(mat,
                           sep = "\t",
                           row.names = 1, 
                           header = TRUE,
                           stringsAsFactors = FALSE, 
                           nrows = 2)
  classes <- sapply(sub_sample, class)
  
  # now actually read in the matrix
  temp <- read.table(mat,
                     sep = "\t",
                     row.names = 1, 
                     header = TRUE,
                     stringsAsFactors = FALSE, 
                     colClasses = classes)
  temp <- as.matrix(temp)
  # Check and return output ----------------------------------------------------
  if (class(temp) != "matrix"){stop("Ouput is incorrectly formatted")}
  return(temp)
} # end read_in_tsv_matrix()


check_is_string <- function(char){
  # Function description ------------------------------------------------------- 
  # Check that the input is a character string. 
  #
  # Input: 
  # char. Character.
  #
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (!is.character(char)){
    stop("Object must be a character string.")
  }
} # end check_is_string()


check_file_exists <- function(file_name){
  # Function description ------------------------------------------------------- 
  # Check that the file exists. 
  #
  # Input: 
  # file_name. Character.   
  #
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (!file.exists(file_name)){
    stop("File does not exist")
  }
} # end check_file_exists()


create_heatmap_compatible_tree <- function(tree){
  heatmap_tree <- tree
  heatmap_tree$edge.length[which(heatmap_tree$edge.length == 0)] <- 0.00001
  heatmap_tree <- chronopl(heatmap_tree,
                           lambda = 0.1,
                           tol = 0)
  heatmap_tree <- as.dendrogram(as.hclust.phylo(heatmap_tree))
  return(heatmap_tree)
} # end create_heatmap_compatible_tree()


assign_pheno_type <- function(mat){
  # Function description ------------------------------------------------------- 
  # Determine if the matrix is discrete or continuous. 
  #
  # Input: 
  # mat: Matrix. Should be a phenotype with 1 column. 
  #
  # Outputs: 
  # type. Character. Either "discrete" or "continuous". 
  
  # Check input ----------------------------------------------------------------
  check_dimensions(mat, NULL, 1, 1, 1)
  
  # Function -------------------------------------------------------------------
  type <- "discrete"
  if (sum(!(mat %in% c(0, 1))) > 0){
    type <- "continuous"
  } 
  
  # Check and return output ----------------------------------------------------
  check_is_string(type)
  return(type)
} # end assign_pheno_type()


check_dimensions <- function(mat, exact_rows, min_rows, exact_cols, min_cols){
  # Function description ------------------------------------------------------- 
  # Check that the matrix is of the specified dimensions. 
  #
  # Input: 
  # mat:        Matrix. 
  # exact_rows. Numeric. Describes expected number of rows in matrix. Can be NULL. 
  # min_rows.   Numeric. Describes minimum number of rows in matrix. Must be specified.  
  # exact_cols. Numeric. Describes expected number of columns in matrix. Can be NULL. 
  # min_cols.   Numeric. Describes minimum number of columns in matrix. Must be specified.  
  #
  # Output: 
  # None. 
  #
  # Check input ----------------------------------------------------------------
  if (!is.null(exact_rows)){check_is_number(exact_rows)}
  check_is_number(min_rows)
  if (!is.null(exact_cols)){check_is_number(exact_cols)}
  check_is_number(min_cols)
  if (class(mat) != "matrix"){stop("input must be a matrix")}
  
  # Function -------------------------------------------------------------------
  if (nrow(mat) < min_rows){
    stop("matrix has too few rows")
  }
  if (ncol(mat) < min_cols){
    stop("matrix has too few columns")
  }
  if (!(is.null(exact_rows))){
    if (nrow(mat) != exact_rows){
      stop("matrix has wrong number of rows")
    }
  }
  if (!(is.null(exact_cols))){
    if (ncol(mat) != exact_cols){
      stop("matrix has wrong number of columns")
    }
  }
} # end check_dimensions()


check_is_number <- function(num){
  # Function description ------------------------------------------------------- 
  # Check that input is some type of number. 
  #
  # Input: 
  # num. Number. Could be numeric, double, or integer.  
  # 
  # Output: 
  # None. 
  #
  # Check input & function -----------------------------------------------------
  if (!is.numeric(num)){
    if (!is.integer(num)){
      if (!is.double(num)){
        stop("Must be a number")
      }
    }
  }
} # end check_is_number()

check_if_dir_exists <- function(dir){
  # Function description ------------------------------------------------------- 
  # Check that output directory exists so data can be saved in it. 
  #
  # Input: 
  # dir. Character. Path to output directory. 
  #
  # Output: 
  # None. 
  #
  # Check input ----------------------------------------------------------------
  check_is_string(dir)
  
  # Function -------------------------------------------------------------------
  if (!dir.exists(dir)){
    stop("The output directory indicated does not exist.")
  }
} # end check_if_dir_exists()

load_rda <- function(file_name){
  #loads an RData file, and returns it
  load(file_name)
  get(ls()[ls() != "file_name"])
} # end load_rda

format_genotype_name <- function(geno_mat){
  colnames(geno_mat) <- gsub("[ ]", ".", colnames(geno_mat))
  colnames(geno_mat) <- gsub("[;]", ".", colnames(geno_mat))
  colnames(geno_mat) <- gsub("[>]", ".", colnames(geno_mat))
  colnames(geno_mat) <- gsub("[=]", ".", colnames(geno_mat))
  colnames(geno_mat) <- gsub("[|]", ".", colnames(geno_mat))
  colnames(geno_mat) <- gsub("[/]", ".", colnames(geno_mat))
  colnames(geno_mat) <- gsub("[*]", ".", colnames(geno_mat))
  colnames(geno_mat) <- gsub("^X", "", colnames(geno_mat))
  return(geno_mat)
}


# END OF SCRIPT ----------------------------------------------------------------
