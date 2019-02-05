# Katie Saund
# 
# Write PBS scripts for treeWAS and transition phyC methods using a standard set of inputs: 
# 
# Standard GWAS inputs: 
# Genotype: 
#   Matrix. 
#   Each row is a sample. 
#   Each column is a genotype. 
#   Entries are binary. 
#
# Phenotype: 
#   Matrix. 
#   Each row is a sample.
#   One column, which is the phenotype. 
#   Can be discrete or continuous. 
#
# Tree
#   Phylo. 
#   Tree$tip.label are the samples. 
#   Phenotype, genotype, and annotation matrix are in the same order as tree$tip.label. 
#   
# Annot
#   Matrix.
#   Each row is a sample.
#   Column 1 is the annotation (ex: ribotype)
#   Column 2 is the color you want the annotation to be. 
#   This is an optional input. 


library(optparse) # Read in command line arguments with flags
# Set up arguments
inputs     <- list(
  make_option(c("-a", "--alpha"),                    type = "double",    default = 0.05,   help = "Value of alpha",                                               metavar = "character"),
  make_option(c("-b", "--reconstruction"),           type = "character", default = "ML",   help = "Ancestral reconstruction method for treewas: ML or parsimony", metavar = "character"),
  make_option(c("-c", "--multiple_test_correction"), type = "character", default = "fdr",  help = "Type of multiple test correction for treewas: fdr of bonf",    metavar = "character"),
  make_option(c("-d", "--bootstrap"),                type = "double",    default = 0.70,   help = "Confidence threshold for bootstrap for phyc",                  metavar = "character"),
  make_option(c("-f", "--generic_filename"),         type = "character", default = NULL,   help = "Name of generic genotype (no extension)",                      metavar = "character"), 
  make_option(c("-g", "--generic_genotype"),         type = "logical",   default = FALSE,  help = "Whether to run GWAS on a generic genotype mat",                metavar = "character"), 
  make_option(c("-i", "--gene_high"),                type = "logical",   default = FALSE,  help = "Whether to run GWAS on stop gene mat",                         metavar = "character"), 
  make_option(c("-k", "--input_dir"),                type = "character", default = NULL,   help = "dir where all formatted data is stored",                       metavar = "character"), 
  make_option(c("-n", "--gene_ns"),                  type = "logical",   default = FALSE,  help = "Whether to run GWAS on gene nonsynonymous mat",                metavar = "character"), 
  make_option(c("-o", "--out"),                      type = "character", default = ".",    help = "Output directory",                                             metavar = "character"),
  make_option(c("-p", "--phenotype"),                type = "character", default = NULL,   help = "phenotype matrix .tsv filename",                               metavar = "character"),
  make_option(c("-q", "--gene_stop"),                type = "logical",   default = FALSE,  help = "Whether to run GWAS on SNPeff high gene mat",                  metavar = "character"), 
  make_option(c("-r", "--roary"),                    type = "logical",   default = FALSE,  help = "Whether to run GWAS on roary",                                 metavar = "character"), 
  make_option(c("-s", "--snp"),                      type = "logical",   default = FALSE,  help = "Whether to run GWAS on SNP mats",                              metavar = "character"), 
  make_option(c("-t", "--treewas"),                  type = "logical",   default = FALSE,  help = "Whether or not to run treewas",                                metavar = "character"),
  make_option(c("-u", "--username"),                 type = "character", default = NULL,   help = "UMich ID",                                                     metavar = "character"),
  make_option(c("-x", "--permutations"),             type = "integer",   default = 10000,  help = "Number of permutations to run in phyc",                        metavar = "character"),
  make_option(c("-y", "--phyc"),                     type = "logical",   default = TRUE,   help = "Whether or not to run phyc",                                   metavar = "character"))

opt_parser <- OptionParser(option_list=inputs)
opt        <- parse_args(opt_parser)


phenotype_mat <-  read.table(opt$phenotype, 
                             sep = "\t",
                             row.names = 1, 
                             header = TRUE,
                             stringsAsFactors = FALSE)
phenotypes <- colnames(phenotype_mat)

genotypes <- NULL
if (!is.null(opt$snp)){
  genotypes <- c(genotypes, "_snp_stop", "_snp_ns", "_snp_high")
}
if (!is.null(opt$gene_stop)){
  genotypes <- c(genotypes, "_gene_stop")
}
if (!is.null(opt$gene_high)){
  genotypes <- c(genotypes, "_gene_high")
}
if (!is.null(opt$gene_ns)){
  genotypes <- c(genotypes, "_gene_ns")
}
if (!is.null(opt$roary)){
  genotypes <- c(genotypes, "_roary_pan_genome")
}
if (!is.null(opt$generic_genotype)){
  genotypes <- c(genotypes, opt$generic_filename)
}

# PHYC COMMAND LINE INPUTS
# 1. Phenotype
# 2. Tree 
# 3. Genotype
# 4. Output name
# 5. Output directory 
# 6. Number of permutations 
# 7. Alpha
# 8. Boostrap confidence threshold 
# 9. Optional: Annotation for heatmap

if (opt$phyc){
  for (p in 1:length(phenotypes)){
    for (g in 1:length(genotypes)){
      if (file.exists(paste(opt$input_dir, phenotypes[p], "_annotation.tsv", sep = ""))){
        annotation <- paste(opt$input_dir, phenotypes[p], "_annotation.tsv", sep = "")
      } else {
        annotation <- ""
      }
      
      
      command <- paste("Rscript /nfs/esnitkin/bin_group/pipeline/Github/gwas/phyc_run.R",
                       paste(opt$input_dir, phenotypes[p], "_pheno.tsv", sep = ""), 
                       paste(opt$input_dir, phenotypes[p], ".tree", sep = ""),
                       paste(opt$input_dir, phenotypes[p], genotypes[g], ".tsv", sep = ""), 
                       paste(phenotypes[p], genotypes[g], sep = ""),
                       getwd(),
                       opt$permutations, 
                       opt$alpha,
                       opt$bootstrap, 
                       annotation,
                       sep = " ")
      fname <- paste(opt$out, "/", "phyc_", phenotypes[p], genotypes[g], ".pbs", sep = "")
      writeLines(c("#!/bin/sh","####  PBS preamble",
                   paste("#PBS -N phyc_", phenotypes[p], genotypes[g], sep = ""),
                   paste("#PBS -M ", opt$username, sep = ""),  
                   "#PBS -m a",
                   "#PBS -l nodes=1:ppn=4,mem=40gb,walltime=150:00:00",
                   "#PBS -V",
                   "#PBS -j oe",
                   "#PBS -A esnitkin_flux",
                   "#PBS -q flux",
                   "#PBS -l qos=flux",
                   "cd $PBS_O_WORKDIR",
                   command),
                 fname,
                 sep = "\n")  
    }
  }
}


# TREEWAS COMMAND LINE INPUTS
# 1. Phenotype
# 2. Tree 
# 3. Genotype
# 4. Reconstruction method
# 5. Multiple test correction method
# 6. Output name
# 7. Alpha value

if (opt$treewas){
  for (p in 1:length(phenotypes)){
    for (g in 1:length(genotypes)){
      command <- paste("Rscript /nfs/esnitkin/bin_group/pipeline/Github/gwas/treewas_run.R",
                       paste(opt$input_dir, phenotypes[p], "_pheno.tsv", sep = ""), 
                       paste(opt$input_dir, phenotypes[p], ".tree", sep = ""),
                       paste(opt$input_dir, phenotypes[p], genotypes[g], ".tsv", sep = ""), 
                       opt$reconstruction,
                       opt$multiple_test_correction,
                       paste("treewas_", phenotypes[p], genotypes[g], sep = ""),
                       opt$alpha, 
                       sep = " ")
      fname <- paste(opt$out, "/", "treewas_", phenotypes[p], genotypes[g], ".pbs", sep = "")
      writeLines(c("#!/bin/sh","####  PBS preamble",
                   paste("#PBS -N treewas_", phenotypes[p], genotypes[g], sep = ""),
                   paste("#PBS -M ", opt$username, sep = ""),  
                   "#PBS -m a",
                   "#PBS -l nodes=1:ppn=4,mem=40gb,walltime=150:00:00",
                   "#PBS -V",
                   "#PBS -j oe",
                   "#PBS -A esnitkin_flux",
                   "#PBS -q flux",
                   "#PBS -l qos=flux",
                   "cd $PBS_O_WORKDIR",
                   command),
                 fname,
                 sep = "\n")  
    }
  }
}

# END SCRIPT -------------------------------------------------------------------
