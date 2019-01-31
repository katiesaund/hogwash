# Katie
# 2019-01-31
# need to figure out something to deal with this cdif specific formatting nonsense. 
#simplify_ribotype <- function(phenotype_data, tree){
#  ribotypes       <- phenotype_data[ , 1, drop = FALSE]
#  ribotypes       <- ribotypes[match(tree$tip.label, row.names(ribotypes)), ]
#  ribotypes[!(ribotypes %in% c("014-020", "027", "078-126", "053-163"))]<- "other"
#  ribo_color <- ribotypes
#  ribo_color[ribo_color == "014-020"] <- "purple"
#  ribo_color[ribo_color == "027"]     <- "blue"
#  ribo_color[ribo_color == "078-126"] <- "green"
#  ribo_color[ribo_color == "053-163"] <- "orange"
#  ribo_color[ribo_color == "other"] <- "grey"
#  ribo_mat <- matrix(NA, nrow = Ntip(tree), ncol = 2)
#  colnames(ribo_mat) <- c("ribotype", "color")
#  row.names(ribo_mat) <- tree$tip.label
#  ribo_mat[ , 1] <- ribotypes
#  ribo_mat[ , 2] <- ribo_color
#  return(ribo_mat)
#} # end simplify_ribotype()


# format_tree <- function(tree, sample_lookup){
#   tree$tip.label = gsub('_genome', "", tree$tip.label); #remove '_genome' at the end of the sample ID
#   tree$tip.label = gsub('^Cdif_', "", tree$tip.label); #remove 'Cdif_' at the start of the sample ID
#   hanna_label <- vector(mode="character", length(tree$tip.label))
#   for (i in 1:length(tree$tip.label)){
#     for (j in 1:nrow(sample_lookup)){
#       if (tree$tip.label[i] == sample_lookup$Unique_ID[j]){
#         hanna_label[i] <- as.character(sample_lookup$Hanna_ID[j])
#       }
#     }
#   }
#   tree$tip.label <- hanna_label
#   tree <- drop.tip(tree, c(1:Ntip(tree))[!(tree$tip.label != "")])
#   return(tree)
# }


# lookup table
#if (!is.null(opt$key) & file.exists(opt$key)){
#  key <- read.table(file = opt$key, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
#  key_exists <- TRUE
#}
#lookup <- read.table("/nfs/esnitkin/Project_Cdiff/Sequence_data/Project_Hanna_collections/meta_data/high_confidence_Hanna_genomes_without_Cd160_controls_nor_Cd025_FOBT142", 
#                     header = TRUE, 
#                     stringsAsFactors = FALSE)


# USE LOOKUP TABLE TO RENAME ISOLATES TO CD/DA INSTEAD OF PH -------------------
# snp_ns   <- convert_PH_to_CD(snp_ns, lookup)
# snp_stop <- convert_PH_to_CD(snp_stop, lookup)
# snp_high <- convert_PH_to_CD(snp_high, lookup)
# 
# gene_ns   <- convert_PH_to_CD(gene_ns, lookup)
# gene_stop <- convert_PH_to_CD(gene_stop, lookup)
# gene_high <- convert_PH_to_CD(gene_high, lookup)
# 
# roary_pan_genome <- convert_PH_to_CD(roary_pan_genome, lookup)


#convert_PH_to_CD <- function(matrix, lookup){
#  for (i in 1:nrow(matrix)){
#    for (j in 1:nrow(lookup)){
#      if (row.names(matrix)[i] == lookup$Unique_ID[j]){
#        row.names(matrix)[i] <- lookup$Hanna_ID[j]
#      }
#    }
#  }
#  return(matrix)
#} # end convert_PH_to_CD()


#remove_bad_samples_from_phenotype <- function(phenotype, lookup){
#  temp <- phenotype[row.names(phenotype) %in% lookup$Hanna_ID, , drop = FALSE]
#  return(temp)
#} # end remove_bad_samples_from_phenotype()


# SUBSET PHENOTYPES TO ONLY HAVE GOOD SAMPLES ----------------------------------
# log_cfe             <- remove_bad_samples_from_phenotype(log_cfe, lookup)            
# CONVERT TREE$TIP.LABEL FROM PH TO CD & DROP BAD SAMPLES --------------------
# tree <- format_tree(tree, lookup)


# # ribotype annotation
# raw_phenotypes <- read.csv("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/data/Hanna_in_vitro_data.txt", 
#                            row.names = 1, 
#                            header = TRUE, 
#                            sep = ",", 
#                            stringsAsFactors = FALSE)
# raw_phenotypes[raw_phenotypes == "ND"] <- NA
# 
# log_cfe_ribo             <- simplify_ribotype(raw_phenotypes, log_cfe_tree)
# log_germ_tc_ribo         <- simplify_ribotype(raw_phenotypes, log_germ_tc_tree)
# log_germ_tc_and_gly_ribo <- simplify_ribotype(raw_phenotypes, log_germ_tc_and_gly_tree)
# log_growth_ribo          <- simplify_ribotype(raw_phenotypes, log_growth_tree)
# log_sporulation_ribo     <- simplify_ribotype(raw_phenotypes, log_sporulation_tree)
# log_toxin_ribo           <- simplify_ribotype(raw_phenotypes, log_toxin_tree)
# fqR_ribo                 <- simplify_ribotype(raw_phenotypes, fqR_tree)
# severity_ribo            <- simplify_ribotype(raw_phenotypes, severity_tree)
