# Katie Saund
# TODO
# Change it from taking in one phenotype at a time to taking in a matrix of 1 or more phenotypes
# Remove any paths that are hard coded and change to command line arguments
# Generalize things that are Cdif specific
# Remove functions from this file and make a lib? 
# Turn this file from a script to a function. 

# Goal is to have genotype, phenotype, and tree data all saved in formats
# that are compatible with treeWAs and phyC.

# INPUTS into this script: 
# 1. Parsed snp mat
# 2. Parsed gene data
# 3. Roary pan genome matrix
# 4. Look up table defining nomenclature and included samples
# 5. Phenotypes (discrete and continuous; continuous will be log transformed)
# 6. Tree

# OUTPUTS from this script: 
# For each phenotype: 
#   1. SNP mats (NS, HIGH, STOPS) (.tsv and .rda)
#   2. Gene mats (NS, HIGH, STOPS) (.tsv and .rda)
#   3. roary mat (.tsv)
#   4. phenotype vector (.tsv)
#   5. Tree (.tree)
#   6. Optional: annotation (.tsv)

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

# FUNCTIONS --------------------------------------------------------------------
read_in_phenotype <- function(path){
  phenotype <- read.csv(path, 
                        row.names = 1, 
                        sep = "\t")
  return(phenotype)
} # end read_in_phenotype()

convert_PH_to_CD <- function(matrix, lookup){
  for (i in 1:nrow(matrix)){
    for (j in 1:nrow(lookup)){
      if (row.names(matrix)[i] == lookup$Unique_ID[j]){
        row.names(matrix)[i] <- lookup$Hanna_ID[j]
      }
    }
  }
  return(matrix)
} # end convert_PH_to_CD()

remove_na_from_phenotype <- function(phenotype){
  temp <- phenotype[!(is.na(phenotype)), , drop = FALSE]
  return(temp)
} # end remove_na_from_phenotype()

remove_inf_from_phenotype <- function(phenotype){
  temp <- phenotype[!(is.infinite(phenotype[ , 1])), , drop = FALSE]
  return(temp)
} # end remove_inf_from_phenotype()

remove_bad_samples_from_phenotype <- function(phenotype, lookup){
  temp <- phenotype[row.names(phenotype) %in% lookup$Hanna_ID, , drop = FALSE]
  return(temp)
} # end remove_bad_samples_from_phenotype()

subset_tree_to_phenotype <- function(phenotype, original_tree){
  new_tree <- drop.tip(phy = original_tree, tip = original_tree$tip.label[!(original_tree$tip.label %in% row.names(phenotype))])
  return(new_tree)
} # end subset_tree_to_phenotype()

create_phenotype_specific_genotype <- function(genotype, phenotype_tree){
  temp <- genotype[match(phenotype_tree$tip.label, row.names(genotype)), , drop = FALSE]
  return(temp)
} # end create_phenotype_specific_genotype

save_results <- function(file_name, phenotype, phenotype_tree, pheno_snp_ns, pheno_snp_stop, pheno_snp_high, 
                         pheno_gene_ns, pheno_gene_stop, pheno_gene_high, pheno_roary_pan_genome, 
                         ribo_annotation, output_dir){
  path <- output_dir
  
  # save treeWAS compatible version of the phenotype: 
  temp <- unlist(phenotype)
  names(temp) <- row.names(phenotype)
  write.table(temp, file = paste(path, file_name, "_pheno.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)
  
  # Save the rest as is
  write.tree(phy = phenotype_tree,    file = paste(path, file_name, ".tree", sep = ""))
  write.table(phenotype,              file = paste(path, file_name, "_pheno.tsv",            sep = ""), sep = "\t", quote = TRUE, row.names = TRUE)
  write.table(pheno_snp_ns,           file = paste(path, file_name, "_snp_ns.tsv",           sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
  write.table(pheno_snp_stop,         file = paste(path, file_name, "_snp_stop.tsv",         sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
  write.table(pheno_snp_high,         file = paste(path, file_name, "_snp_high.tsv",         sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
  write.table(pheno_gene_ns,          file = paste(path, file_name, "_gene_ns.tsv",          sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
  write.table(pheno_gene_stop,        file = paste(path, file_name, "_gene_stop.tsv",        sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
  write.table(pheno_gene_high,        file = paste(path, file_name, "_gene_high.tsv",        sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
  write.table(pheno_roary_pan_genome, file = paste(path, file_name, "_roary_pan_genome.tsv", sep = ""), sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
  write.table(ribo_annotation,        file = paste(path, file_name, "_annotation.tsv", sep = ""),       sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
  
  # And save the genotypes as a .rda file 
  save(pheno_snp_ns,           file = paste(path, file_name, "_snp_ns.rda",           sep = ""))
  save(pheno_snp_stop,         file = paste(path, file_name, "_snp_stop.rda",         sep = ""))
  save(pheno_snp_high,         file = paste(path, file_name, "_snp_high.rda",         sep = ""))
  save(pheno_gene_ns,          file = paste(path, file_name, "_gene_ns.rda",          sep = ""))
  save(pheno_gene_stop,        file = paste(path, file_name, "_gene_stop.rda",        sep = ""))
  save(pheno_gene_high,        file = paste(path, file_name, "_gene_high.rda",        sep = ""))
  save(pheno_roary_pan_genome, file = paste(path, file_name, "_roary_pan_genome.rda", sep = ""))
  
  
} # end save_results()

simplify_ribotype <- function(phenotype_data, tree){
  ribotypes       <- phenotype_data[ , 1, drop = FALSE]
  ribotypes       <- ribotypes[match(tree$tip.label, row.names(ribotypes)), ]
  ribotypes[!(ribotypes %in% c("014-020", "027", "078-126", "053-163"))]<- "other"
  ribo_color <- ribotypes
  ribo_color[ribo_color == "014-020"] <- "purple"
  ribo_color[ribo_color == "027"]     <- "blue"
  ribo_color[ribo_color == "078-126"] <- "green"
  ribo_color[ribo_color == "053-163"] <- "orange"
  ribo_color[ribo_color == "other"] <- "grey"
  ribo_mat <- matrix(NA, nrow = Ntip(tree), ncol = 2)
  colnames(ribo_mat) <- c("ribotype", "color")
  row.names(ribo_mat) <- tree$tip.label
  ribo_mat[ , 1] <- ribotypes
  ribo_mat[ , 2] <- ribo_color
  return(ribo_mat)
} # end simplify_ribotype()


format_tree <- function(tree, sample_lookup){
  tree$tip.label = gsub('_genome', "", tree$tip.label); #remove '_genome' at the end of the sample ID
  tree$tip.label = gsub('^Cdif_', "", tree$tip.label); #remove 'Cdif_' at the start of the sample ID
  hanna_label <- vector(mode="character", length(tree$tip.label))
  for (i in 1:length(tree$tip.label)){
    for (j in 1:nrow(sample_lookup)){
      if (tree$tip.label[i] == sample_lookup$Unique_ID[j]){
        hanna_label[i] <- as.character(sample_lookup$Hanna_ID[j])
      }
    }
  }
  tree$tip.label <- hanna_label
  return(tree)
}

# LIBRARIES --------------------------------------------------------------------
library(ape)

# INPUTS -----------------------------------------------------------------------
# PHENOTYPES
log_cfe             <- read_in_phenotype("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/data/Hanna_in_vitro_data_continuous_phenotypes/log_cfe")
log_germ_tc         <- read_in_phenotype("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/data/Hanna_in_vitro_data_continuous_phenotypes/log_germ_tc")
log_germ_tc_and_gly <- read_in_phenotype("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/data/Hanna_in_vitro_data_continuous_phenotypes/log_germ_tc_and_gly")
log_growth          <- read_in_phenotype("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/data/Hanna_in_vitro_data_continuous_phenotypes/log_growth")
log_sporulation     <- read_in_phenotype("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/data/Hanna_in_vitro_data_continuous_phenotypes/log_sporulation")
log_toxin           <- read_in_phenotype("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/data/Hanna_in_vitro_data_continuous_phenotypes/log_toxin")
fqR                 <- read.csv("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2017-08-18_comprehensive_phenotype_binning/data/all_bins/unsorted/gyrA_res_mut", sep = "", row.names = 1, header = TRUE)
severity            <- read.csv("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2017-08-18_comprehensive_phenotype_binning/data/all_bins/unsorted/severity", sep = "", row.names = 1, header = TRUE)

# TREE
tree <- read.tree("/nfs/esnitkin/Project_Cdiff/Sequence_data/Project_Hanna_collections/consensus/2019_01_11_09_42_32_core_results/gubbins/raxml_results/RAxML_bipartitionsINSERT REST OF FILE NAME HERE")
# this particular version of the tree has a trailing "_" in some tip.labels. Remove here: 
tree$tip.label <- gsub("_$", "", tree$tip.label)

# GENOTYPES
roary_pan_genome <- read.csv("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-05-23_roary_id_70_split_paralogs/data/gene_presence_absence.Rtab", 
                             sep = "\t", 
                             header = TRUE, 
                             row.names = 1)
 
load("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-01-22_snpeff_high_impact_gene_matrix.rda")
gene_high <- gene_matrix
gene_matrix <- NULL


load("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-01-22_stop_gene_matrix.rda")
gene_stop <- gene_matrix
gene_matrix <- NULL

load("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-01-22_nonsense_gene_matrix.rda")
gene_ns <- gene_matrix
gene_matrix <- NULL

load("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-01-22_parsed_simple_code_snpmat.rda")
snp_ns   <- parsed_simple_code_snpmat$snpmat[parsed_simple_code_snpmat$ns_mut,      , drop = FALSE]
snp_stop <- parsed_simple_code_snpmat$snpmat[parsed_simple_code_snpmat$stop,        , drop = FALSE]
snp_high <- parsed_simple_code_snpmat$snpmat[parsed_simple_code_snpmat$snpeff_high, , drop = FALSE]

# lookup table
lookup <- read.table("/nfs/esnitkin/Project_Cdiff/Sequence_data/Project_Hanna_collections/meta_data/high_confidence_Hanna_genomes_without_Cd160_controls_nor_Cd025_FOBT142", 
                     header = TRUE, 
                     stringsAsFactors = FALSE)

# ORIENT ALL GENOTYPE MATRICES TO BE ISOLATES IN ROWS AND LOCI IN COLUMNS ------
snp_ns   <- t(snp_ns)
snp_stop <- t(snp_stop)
snp_high <- t(snp_high)
gene_ns   <- t(gene_ns)
gene_stop <- t(gene_stop)
gene_high <- t(gene_high)
roary_pan_genome <- t(roary_pan_genome)

# USE LOOKUP TABLE TO RENAME ISOLATES TO CD/DA INSTEAD OF PH -------------------
snp_ns   <- convert_PH_to_CD(snp_ns, lookup)
snp_stop <- convert_PH_to_CD(snp_stop, lookup)
snp_high <- convert_PH_to_CD(snp_high, lookup)

gene_ns   <- convert_PH_to_CD(gene_ns, lookup)
gene_stop <- convert_PH_to_CD(gene_stop, lookup)
gene_high <- convert_PH_to_CD(gene_high, lookup)

roary_pan_genome <- convert_PH_to_CD(roary_pan_genome, lookup)

# SUBSET PHENOTYPES TO HAVE NO NAs ---------------------------------------------
log_cfe             <- remove_na_from_phenotype(log_cfe)            
log_germ_tc         <- remove_na_from_phenotype(log_germ_tc)      
log_germ_tc_and_gly <- remove_na_from_phenotype(log_germ_tc_and_gly)
log_growth          <- remove_na_from_phenotype(log_growth)        
log_sporulation     <- remove_na_from_phenotype(log_sporulation)   
log_toxin           <- remove_na_from_phenotype(log_toxin)        
fqR                 <- remove_na_from_phenotype(fqR)            
severity            <- remove_na_from_phenotype(severity)       

# SUBSET PHENOTYPES TO HAVE NO INF/-INF ----------------------------------------
log_cfe             <- remove_inf_from_phenotype(log_cfe)            
log_germ_tc         <- remove_inf_from_phenotype(log_germ_tc)      
log_germ_tc_and_gly <- remove_inf_from_phenotype(log_germ_tc_and_gly)
log_growth          <- remove_inf_from_phenotype(log_growth)        
log_sporulation     <- remove_inf_from_phenotype(log_sporulation)   
log_toxin           <- remove_inf_from_phenotype(log_toxin)        
fqR                 <- remove_inf_from_phenotype(fqR)            
severity            <- remove_inf_from_phenotype(severity)  

# SUBSET PHENOTYPES TO ONLY HAVE GOOD SAMPLES ----------------------------------
log_cfe             <- remove_bad_samples_from_phenotype(log_cfe, lookup)            
log_germ_tc         <- remove_bad_samples_from_phenotype(log_germ_tc, lookup)      
log_germ_tc_and_gly <- remove_bad_samples_from_phenotype(log_germ_tc_and_gly, lookup)
log_growth          <- remove_bad_samples_from_phenotype(log_growth, lookup)        
log_sporulation     <- remove_bad_samples_from_phenotype(log_sporulation, lookup)   
log_toxin           <- remove_bad_samples_from_phenotype(log_toxin, lookup)        
fqR                 <- remove_bad_samples_from_phenotype(fqR, lookup)            
severity            <- remove_bad_samples_from_phenotype(severity, lookup)      


# CONVERT TREE$TIP.LABEL FROM PH TO CD & DROP BAD SAMPLES ----------------------
tree <- format_tree(tree, lookup)
tree <- drop.tip(tree, c(1:Ntip(tree))[!(tree$tip.label != "")])

# MIDPOINT ROOT TREE -----------------------------------------------------------
tree <- midpoint.root(tree)

# CREATE PHENOTYPE SPECIFIC TREES ----------------------------------------------
log_cfe_tree             <- subset_tree_to_phenotype(log_cfe, tree)
log_germ_tc_tree         <- subset_tree_to_phenotype(log_germ_tc, tree)
log_germ_tc_and_gly_tree <- subset_tree_to_phenotype(log_germ_tc_and_gly, tree)
log_growth_tree          <- subset_tree_to_phenotype(log_growth, tree)
log_sporulation_tree     <- subset_tree_to_phenotype(log_sporulation, tree)
log_toxin_tree           <- subset_tree_to_phenotype(log_toxin, tree)
fqR_tree                 <- subset_tree_to_phenotype(fqR, tree)
severity_tree            <- subset_tree_to_phenotype(severity, tree)

# REORDER PHENOTYPES TO MATCH TREE ORDER ---------------------------------------
log_cfe             <- create_phenotype_specific_genotype(log_cfe,             log_cfe_tree)
log_germ_tc         <- create_phenotype_specific_genotype(log_germ_tc,         log_germ_tc_tree)
log_germ_tc_and_gly <- create_phenotype_specific_genotype(log_germ_tc_and_gly, log_germ_tc_and_gly_tree)
log_growth          <- create_phenotype_specific_genotype(log_growth,          log_growth_tree)
log_sporulation     <- create_phenotype_specific_genotype(log_sporulation,     log_sporulation_tree)
log_toxin           <- create_phenotype_specific_genotype(log_toxin,           log_toxin_tree)
fqR                 <- create_phenotype_specific_genotype(fqR,                 fqR_tree)
severity            <- create_phenotype_specific_genotype(severity,            severity_tree)

# CREATE PHENOTYPE SPECIFIC GENOTYPES ------------------------------------------

# CFE
cfe_snp_ns   <- create_phenotype_specific_genotype(snp_ns, log_cfe_tree)
cfe_snp_stop <- create_phenotype_specific_genotype(snp_stop, log_cfe_tree)
cfe_snp_high <- create_phenotype_specific_genotype(snp_high, log_cfe_tree)

cfe_gene_ns   <- create_phenotype_specific_genotype(gene_ns, log_cfe_tree)
cfe_gene_stop <- create_phenotype_specific_genotype(gene_stop, log_cfe_tree)
cfe_gene_high <- create_phenotype_specific_genotype(gene_high, log_cfe_tree)

cfe_roary_pan_genome <- create_phenotype_specific_genotype(roary_pan_genome, log_cfe_tree)

# GERM_TC
germ_tc_snp_ns   <- create_phenotype_specific_genotype(snp_ns,   log_germ_tc_tree)
germ_tc_snp_stop <- create_phenotype_specific_genotype(snp_stop, log_germ_tc_tree)
germ_tc_snp_high <- create_phenotype_specific_genotype(snp_high, log_germ_tc_tree)

germ_tc_gene_ns   <- create_phenotype_specific_genotype(gene_ns,   log_germ_tc_tree)
germ_tc_gene_stop <- create_phenotype_specific_genotype(gene_stop, log_germ_tc_tree)
germ_tc_gene_high <- create_phenotype_specific_genotype(gene_high, log_germ_tc_tree)

germ_tc_roary_pan_genome <- create_phenotype_specific_genotype(roary_pan_genome, log_germ_tc_tree)
  
# GERM_TC_AND_GLY
germ_tc_and_gly_snp_ns   <- create_phenotype_specific_genotype(snp_ns,   log_germ_tc_and_gly_tree)
germ_tc_and_gly_snp_stop <- create_phenotype_specific_genotype(snp_stop, log_germ_tc_and_gly_tree)
germ_tc_and_gly_snp_high <- create_phenotype_specific_genotype(snp_high, log_germ_tc_and_gly_tree)

germ_tc_and_gly_gene_ns   <- create_phenotype_specific_genotype(gene_ns,   log_germ_tc_and_gly_tree)
germ_tc_and_gly_gene_stop <- create_phenotype_specific_genotype(gene_stop, log_germ_tc_and_gly_tree)
germ_tc_and_gly_gene_high <- create_phenotype_specific_genotype(gene_high, log_germ_tc_and_gly_tree)

germ_tc_and_gly_roary_pan_genome <- create_phenotype_specific_genotype(roary_pan_genome, log_germ_tc_and_gly_tree)

# GROWTH
growth_snp_ns   <- create_phenotype_specific_genotype(snp_ns,   log_growth_tree)
growth_snp_stop <- create_phenotype_specific_genotype(snp_stop, log_growth_tree)
growth_snp_high <- create_phenotype_specific_genotype(snp_high, log_growth_tree)

growth_gene_ns   <- create_phenotype_specific_genotype(gene_ns,   log_growth_tree)
growth_gene_stop <- create_phenotype_specific_genotype(gene_stop, log_growth_tree)
growth_gene_high <- create_phenotype_specific_genotype(gene_high, log_growth_tree)

growth_roary_pan_genome <- create_phenotype_specific_genotype(roary_pan_genome, log_growth_tree)

# SPORULATION
sporulation_snp_ns   <- create_phenotype_specific_genotype(snp_ns,   log_sporulation_tree)
sporulation_snp_stop <- create_phenotype_specific_genotype(snp_stop, log_sporulation_tree)
sporulation_snp_high <- create_phenotype_specific_genotype(snp_high, log_sporulation_tree)

sporulation_gene_ns   <- create_phenotype_specific_genotype(gene_ns,   log_sporulation_tree)
sporulation_gene_stop <- create_phenotype_specific_genotype(gene_stop, log_sporulation_tree)
sporulation_gene_high <- create_phenotype_specific_genotype(gene_high, log_sporulation_tree)

sporulation_roary_pan_genome <- create_phenotype_specific_genotype(roary_pan_genome, log_sporulation_tree)

# TOXIN
toxin_snp_ns   <- create_phenotype_specific_genotype(snp_ns,   log_toxin_tree)
toxin_snp_stop <- create_phenotype_specific_genotype(snp_stop, log_toxin_tree)
toxin_snp_high <- create_phenotype_specific_genotype(snp_high, log_toxin_tree)

toxin_gene_ns   <- create_phenotype_specific_genotype(gene_ns,   log_toxin_tree)
toxin_gene_stop <- create_phenotype_specific_genotype(gene_stop, log_toxin_tree)
toxin_gene_high <- create_phenotype_specific_genotype(gene_high, log_toxin_tree)

toxin_roary_pan_genome <- create_phenotype_specific_genotype(roary_pan_genome, log_toxin_tree)

# FQR
fqR_snp_ns   <- create_phenotype_specific_genotype(snp_ns,   fqR_tree)
fqR_snp_stop <- create_phenotype_specific_genotype(snp_stop, fqR_tree)
fqR_snp_high <- create_phenotype_specific_genotype(snp_high, fqR_tree)

fqR_gene_ns   <- create_phenotype_specific_genotype(gene_ns,   fqR_tree)
fqR_gene_stop <- create_phenotype_specific_genotype(gene_stop, fqR_tree)
fqR_gene_high <- create_phenotype_specific_genotype(gene_high, fqR_tree)

fqR_roary_pan_genome <- create_phenotype_specific_genotype(roary_pan_genome, fqR_tree)

# SEVERITY
severity_snp_ns   <- create_phenotype_specific_genotype(snp_ns,   severity_tree)
severity_snp_stop <- create_phenotype_specific_genotype(snp_stop, severity_tree)
severity_snp_high <- create_phenotype_specific_genotype(snp_high, severity_tree)

severity_gene_ns   <- create_phenotype_specific_genotype(gene_ns,   severity_tree)
severity_gene_stop <- create_phenotype_specific_genotype(gene_stop, severity_tree)
severity_gene_high <- create_phenotype_specific_genotype(gene_high, severity_tree)

severity_roary_pan_genome <- create_phenotype_specific_genotype(roary_pan_genome, severity_tree)

# ribotype annotation
raw_phenotypes <- read.csv("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/data/Hanna_in_vitro_data.txt", 
                           row.names = 1, 
                           header = TRUE, 
                           sep = ",", 
                           stringsAsFactors = FALSE)
raw_phenotypes[raw_phenotypes == "ND"] <- NA

log_cfe_ribo             <- simplify_ribotype(raw_phenotypes, log_cfe_tree)
log_germ_tc_ribo         <- simplify_ribotype(raw_phenotypes, log_germ_tc_tree)
log_germ_tc_and_gly_ribo <- simplify_ribotype(raw_phenotypes, log_germ_tc_and_gly_tree)
log_growth_ribo          <- simplify_ribotype(raw_phenotypes, log_growth_tree)
log_sporulation_ribo     <- simplify_ribotype(raw_phenotypes, log_sporulation_tree)
log_toxin_ribo           <- simplify_ribotype(raw_phenotypes, log_toxin_tree)
fqR_ribo                 <- simplify_ribotype(raw_phenotypes, fqR_tree)
severity_ribo            <- simplify_ribotype(raw_phenotypes, severity_tree)

# SAVE EVERYTHING --------------------------------------------------------------
save_results("log_cfe", log_cfe, log_cfe_tree, cfe_snp_ns, cfe_snp_stop, cfe_snp_high, cfe_gene_ns, cfe_gene_stop, cfe_gene_high, cfe_roary_pan_genome, log_cfe_ribo, output_dir)
save_results("log_germ_tc", log_germ_tc, log_germ_tc_tree, germ_tc_snp_ns, germ_tc_snp_stop, germ_tc_snp_high, germ_tc_snp_del, germ_tc_gene_ns, germ_tc_gene_stop, germ_tc_gene_high, germ_tc_gene_del, germ_tc_roary_pan_genome, log_germ_tc_ribo, output_dir)
save_results("log_germ_tc_and_gly", log_germ_tc_and_gly, log_germ_tc_and_gly_tree, germ_tc_and_gly_snp_ns, germ_tc_and_gly_snp_stop, germ_tc_and_gly_snp_high, germ_tc_and_gly_snp_del, germ_tc_and_gly_gene_ns, germ_tc_and_gly_gene_stop, germ_tc_and_gly_gene_high, germ_tc_and_gly_gene_del, germ_tc_and_gly_roary_pan_genome, log_germ_tc_and_gly_ribo, output_dir)
save_results("log_growth", log_growth, log_growth_tree, growth_snp_ns, growth_snp_stop, growth_snp_high, growth_snp_del, growth_gene_ns, growth_gene_stop, growth_gene_high, growth_gene_del, growth_roary_pan_genome, log_growth_ribo, output_dir)
save_results("log_sporulation", log_sporulation, log_sporulation_tree, sporulation_snp_ns, sporulation_snp_stop, sporulation_snp_high, sporulation_snp_del, sporulation_gene_ns, sporulation_gene_stop, sporulation_gene_high, sporulation_gene_del, sporulation_roary_pan_genome, log_sporulation_ribo, output_dir)
save_results("log_toxin", log_toxin, log_toxin_tree, toxin_snp_ns, toxin_snp_stop, toxin_snp_high, toxin_snp_del, toxin_gene_ns, toxin_gene_stop, toxin_gene_high, toxin_gene_del, toxin_roary_pan_genome, log_toxin_ribo, output_dir)
save_results("fqR", fqR, fqR_tree, fqR_snp_ns, fqR_snp_stop, fqR_snp_high, fqR_snp_del, fqR_gene_ns, fqR_gene_stop, fqR_gene_high, fqR_gene_del, fqR_roary_pan_genome, fqR_ribo, output_dir)
save_results("severity", severity, severity_tree, severity_snp_ns, severity_snp_stop, severity_snp_high, severity_snp_del, severity_gene_ns, severity_gene_stop, severity_gene_high, severity_gene_del, severity_roary_pan_genome, severity_ribo, output_dir)

# END OF SCRIPT ----------------------------------------------------------------