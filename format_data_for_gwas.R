# Katie Saund

# Goal is to have genotype, phenotype, and tree data all saved in formats
# that are compatible with treeWAs and phyC.

# INPUTS into this script: 
# 1. Parsed snp mat (optional)
# 2. Parsed gene data (optional)
# 3. Roary pan genome matrix (optional)
# 4. Phenotypes (discrete and/or continuous) (required)
# 6. Tree (required)

# OUTPUTS from this script: 
# For each phenotype: 
#   1. SNP mats (NS, HIGH, STOPS) (.tsv and .rda)
#   2. Gene mats (NS, HIGH, STOPS) (.tsv and .rda)
#   3. roary mat (.tsv)
#   4. phenotype vector (.tsv)
#   5. Tree (.tree)

# Naming conventions: 
#   Ex: toxin
#       toxin_snp_stop.tsv
#       toxin_snp_high.tsv
#       toxin_snp_ns.tsv
#       toxin_gene_stop.tsv
#       toxin_gene_high.tsv
#       toxin_gene_ns.tsv
#       toxin_roary_pan_genome.tsv
#       toxin_pheno.tsv
#       toxin_tree.tree

# FUNCTIONS -------------------------------------------------------------------#
format_phenotypes <- function(pheno_mat){
  # break pheno_mat into a list of matrices, each entry in the list is now a matrix corresponding to one phenotype
  pheno_list <- rep(list(matrix(NA)), ncol(pheno_mat))
  for (c in 1:ncol(pheno_mat)){
    pheno_list[[c]] <- pheno_mat[ , c, drop = FALSE]
  }
  for (l in 1:length(pheno_list)){
    pheno_list[[l]] <- remove_na_from_phenotype(pheno_list[[l]])
    pheno_list[[l]] <- remove_inf_from_phenotype(pheno_list[[l]])
  }
  return(pheno_list)
} # end format_phenotypes()

remove_na_from_phenotype <- function(phenotype){
  temp <- phenotype[!(is.na(phenotype)), , drop = FALSE]
  return(temp)
} # end remove_na_from_phenotype()

remove_inf_from_phenotype <- function(phenotype){
  temp <- phenotype[!(is.infinite(phenotype[ , 1])), , drop = FALSE]
  return(temp)
} # end remove_inf_from_phenotype()

subset_tree_to_phenotype <- function(pheno_list, tr){
  tree_list <- rep(list(tr), length(pheno_list))
  
  for (l in 1:length(pheno_list)){
    tree_list[[l]] <- drop.tip(phy = tr, tip = tr$tip.label[!(tr$tip.label %in% row.names(pheno_list[[l]]))])
  }
  return(tree_list)
} # end subset_tree_to_phenotype()

order_phenotype_to_tree <- function(geno_list, tr_list){
  if (length(tr_list) != length(geno_list)){
    stop("Different lengths")
  }
  for (l in 1:length(tr_list)){
    geno_list[[l]] <- geno_list[[l]][match(tr_list[[l]]$tip.label, row.names(geno_list[[l]])), , drop = FALSE]
  }
  return(geno_list)
} # end order_phenotype_to_tree

create_phenotype_specific_genotype <- function(genotype, tr_list){
  if (sum(row.names(genotype) %in% tr_list[[1]]$tip.label) < 1){ # transpose to get isolates in rows and loci in columns
    genotype <- t(genotype)
  }
  if (sum(row.names(genotype) %in% tr_list[[1]]$tip.label) < 1){
    stop("Sample IDs in genotype do not match any tree tip labels")
  }
  geno_list <- rep(list(matrix(NA)), length(tr_list))
  for (l in 1:length(tr_list)){
    geno_list[[l]] <- genotype[match(tr_list[[l]]$tip.label, row.names(genotype)), , drop = FALSE]
  }
  return(geno_list)
} # end create_phenotype_specific_genotype

save_results <- function(pheno_list, tr_list, s_ns, s_stop, s_high, s_ex, 
                         g_ns, g_ns_ex, g_stop, g_stop_ex, g_high, g_high_ex, 
                         roary, roary_ex, output_dir){
  
  for (l in 1:length(pheno_list)){
    file_name <- colnames(pheno_list[[l]])[1]
    write.table(pheno_list[[l]],   file = paste(output_dir, "/", file_name, "_pheno.tsv", sep = ""), sep = "\t", quote = TRUE, row.names = TRUE)
    write.tree(phy = tr_list[[l]], file = paste(output_dir, "/", file_name, ".tree",      sep = ""))
    
    if (s_ex){
      temp_ns   <- s_ns[[l]]
      temp_stop <- s_stop[[l]]
      temp_high <- s_high[[l]]
      write.table(temp_ns,    file = paste(output_dir, "/", file_name, "_snp_ns.tsv",   sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
      write.table(temp_stop,  file = paste(output_dir,  "/", file_name, "_snp_stop.tsv", sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
      write.table(temp_high,  file = paste(output_dir,  "/", file_name, "_snp_high.tsv", sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
      save(temp_ns,           file = paste(output_dir,  "/", file_name, "_snp_ns.rda",   sep = ""))
      save(temp_stop,         file = paste(output_dir,  "/", file_name, "_snp_stop.rda", sep = ""))
      save(temp_high,         file = paste(output_dir,  "/", file_name, "_snp_high.rda", sep = ""))
    }
    if (g_ns_ex){
      temp_ns   <- g_ns[[l]]
      write.table(temp_ns,   file = paste(output_dir, "/",  file_name, "_gene_ns.tsv", sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
      save(temp_ns,          file = paste(output_dir,  "/", file_name, "_gene_ns.rda",          sep = ""))
    }
    if (g_stop_ex){
      temp_stop <- g_stop[[l]]
      write.table(temp_stop, file = paste(output_dir,  "/", file_name, "_gene_stop.tsv", sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
      save(temp_stop,        file = paste(output_dir,  "/", file_name, "_gene_stop.rda",        sep = ""))
    }
    if (g_high_ex){
      temp_high <- g_high[[l]]
      write.table(temp_high, file = paste(output_dir,  "/", file_name, "_gene_high.tsv", sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
      save(temp_high,        file = paste(output_dir,  "/", file_name, "_gene_high.rda",        sep = ""))
    }
    if (roary_ex){
      temp_roary <- roary[[l]]
      write.table(temp_roary, file = paste(output_dir,  "/", file_name, "_roary_pan_genome.tsv", sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
      save(temp_roary,        file = paste(output_dir,  "/", file_name, "_roary_pan_genome.rda", sep = ""))
    }
  }
} # end save_results()

format_data <- function(opt){
  # PHENOTYPES
  if (is.null(opt$phenotype)){
    stop("Must include phenotype matrix")
  } else {
    phenotype_matrix <- read.csv(file = opt$phenotype, header = TRUE, sep = "\t", row.names = 1) 
  }
  
  # TREE
  if (is.null(opt$tree)){
    stop("Must include tree")
  } else {
    tree <- read.tree(file = opt$tree)
  }
  
  # GENOTYPES
  # Set defaults to don't exist (FALSE), but update to TRUE if files are given. 
  roary_exists <- gene_high_exists <- gene_stop_exists <- gene_ns_exists <- snp_exists <- FALSE
  
  if (!is.null(opt$roary)){
    if (file.exists(opt$roary)){
      roary_pan_genome <- read.csv(file = opt$roary, sep = "\t", header = TRUE, row.names = 1)
      roary_exists <- TRUE
    }
  }
  
  if (!is.null(opt$genehigh)){
    if (file.exists(opt$genehigh)){
      gene_high <- local(get(load(opt$genehigh)))
      gene_high_exists <- TRUE
    }
  }
  
  if (!is.null(opt$genestop)){
    if (file.exists(opt$genestop)){
      gene_stop <- local(get(load(opt$genestop)))
      gene_stop_exists <- TRUE
    }
  }
  
  if (!is.null(opt$genens)){
    if (file.exists(opt$genens)){
      gene_ns <-local(get(load(opt$genens)))
      gene_ns_exists <- TRUE
    }
  }
  
  if (!is.null(opt$snp)){
    if (file.exists(opt$snp)){
      snp      <- local(get(load(opt$snp)))
      snp_ns   <- snp$snpmat[snp$ns_mut,      , drop = FALSE]
      snp_stop <- snp$snpmat[snp$stop,        , drop = FALSE]
      snp_high <- snp$snpmat[snp$snpeff_high, , drop = FALSE]
      snp_exists <- TRUE
    }
  }
  
  # SUBSET PHENOTYPES TO HAVE NO NAs & NO INF/-INF ------------------------------#
  phenotype_list <- format_phenotypes(phenotype_matrix)
  
  # MIDPOINT ROOT TREE IF NOT ALREADY ROOTED ------------------------------------#
  if (!is.rooted(tree)){
    tree <- midpoint.root(tree)
  }
  
  # CREATE PHENOTYPE SPECIFIC TREES ---------------------------------------------#
  tree_list <- subset_tree_to_phenotype(phenotype_list, tree)
  
  # REORDER PHENOTYPES TO MATCH TREE ORDER --------------------------------------#
  phenotype_list <- order_phenotype_to_tree(phenotype_list, tree_list)
  
  # CREATE PHENOTYPE SPECIFIC GENOTYPES -----------------------------------------#
  if (roary_exists){
    roary_list <- create_phenotype_specific_genotype(roary_pan_genome, tree_list)
  }
  if (gene_high_exists){
    gene_high_list <- create_phenotype_specific_genotype(gene_high, tree_list)
  }
  if (gene_stop_exists){
    gene_stop_list <- create_phenotype_specific_genotype(gene_stop, tree_list)
  }
  if (gene_ns_exists){
    gene_ns_list <- create_phenotype_specific_genotype(gene_ns, tree_list)
  }
  if (snp_exists){
    snp_ns_list   <- create_phenotype_specific_genotype(snp_ns,   tree_list)
    snp_stop_list <- create_phenotype_specific_genotype(snp_stop, tree_list)
    snp_high_list <- create_phenotype_specific_genotype(snp_high, tree_list)
  }
  
  # OUTPUT ----------------------------------------------------------------------#
  save_results(phenotype_list, tree_list, snp_ns_list, snp_stop_list, 
               snp_high_list, snp_exists, gene_ns_list, gene_ns_exists, 
               gene_stop_list, gene_stop_exists, gene_high_list, gene_high_exists, 
               roary_list, roary_exists, opt$out)
} #end format_data()

# LIBRARIES -------------------------------------------------------------------#
library(ape)      # Phylogenetic trees
library(optparse) # Read in command line arguments with flags

# INPUTS ----------------------------------------------------------------------#
# Set up arguments
inputs     <- list(make_option(c("-p", "--phenotype"), type = "character", default = NULL, help = "path to phenotype matrix .tsv file",         metavar = "character"), 
                   make_option(c("-t", "--tree"),      type = "character", default = NULL, help = "path to phylogenetic tree .tre file",        metavar = "character"), 
                   make_option(c("-r", "--roary"),     type = "character", default = NULL, help = "path to roary pangenome .Rtab file",         metavar = "character"), 
                   make_option(c("-q", "--genestop"),  type = "character", default = NULL, help = "path to SNPeff high gene matrix .rda file",  metavar = "character"), 
                   make_option(c("-i", "--genehigh"),  type = "character", default = NULL, help = "path to stop gene matrix .rda file",         metavar = "character"), 
                   make_option(c("-n", "--genens"),    type = "character", default = NULL, help = "path to gene nonsynonymous .rda file",       metavar = "character"), 
                   make_option(c("-s", "--snp"),       type = "character", default = NULL, help = "path to simple SNP matrix .rda file",        metavar = "character"), 
                   make_option(c("-o", "--out"),       type = "character", default = NULL, help = "output directory",                           metavar = "character"))
opt_parser <- OptionParser(option_list=inputs)
opt        <- parse_args(opt_parser)
format_data(opt)

# END OF SCRIPT ---------------------------------------------------------------#