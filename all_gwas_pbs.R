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

require(tools)
# Read in arguments
args <- commandArgs(trailingOnly = TRUE) # Grab arguments from the PBS script

# args[1] is path to script that runs phyC ex: /nfs/esnitkin/bin_group/pipeline/Github/gwas/phyc_run.R
# args[2] is alpha value for phyC
# args[3] is number of permutations for phyC
# args[4] is the bootstrap confidence threshold (values below this are called low confidence edges)

# args[5] is path to script that runs treewas ex: /nfs/esnitkin/bin_group/pipeline/Github/gwas/treewas_run.R
# args[6] is type of multiple test correction ex: "fdr" or "bonf"
# args[7] is type of ancestral reconstruction method ex: "ML" or "parsimony" 
# args[8] is the alpha value for treewas

# args[9] is path to formatted data on which you want to run GWAS ex: /nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/
#           phenotypes and genotype nomenclature will need to be dealt with still at a later time

# args[10] is your umich email address 

# PHYC SPECIFIC
run_phyC  <- args[1]
alpha     <- args[2]
num_perm  <- args[3]
boot_conf <- args[4]

# TREEWAS SPECIFIC
run_treewas <- args[5]
test_corr   <- args[6]
recon       <- args[7]
alpha_val   <- args[8]

# FOR BOTH: 
path <- args[9]
phenotypes <- c("log_cfe", "log_germ_tc", "log_germ_tc_and_gly", "log_growth", "log_sporulation", "log_toxin", "fqR", "severity") 
genotypes  <- c("_snp_stop", "_snp_ns", "_snp_high", "_gene_stop", "_gene_ns", "_gene_high", "_roary_pan_genome")

username <- args[10]

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

for (p in 1:length(phenotypes)){
  for (g in 1:length(genotypes)){
    if (file.exists(paste(path, phenotypes[p], "_annotation.tsv", sep = ""))){
      annotation <- paste(path, phenotypes[p], "_annotation.tsv", sep = "")
    } else {
      annotation <- ""
    }

    
    command <- paste("Rscript",
                     run_phyC,
                     paste(path, phenotypes[p], "_pheno.tsv", sep = ""), 
                     paste(path, phenotypes[p], ".tree", sep = ""),
                     paste(path, phenotypes[p], genotypes[g], ".tsv", sep = ""), 
                     paste(phenotypes[p], genotypes[g], sep = ""),
                     getwd(),
                     num_perm, 
                     alpha,
                     boot_conf, 
                     annotation,
                     sep = " ")
    fname <- paste(getwd(), "/", "phyc_", phenotypes[p], genotypes[g], ".pbs", sep = "")
    writeLines(c("#!/bin/sh","####  PBS preamble",
                 paste("#PBS -N phyc_", phenotypes[p], genotypes[g], sep = ""),
                 paste("#PBS -M ", username, sep = ""),  
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

# TREEWAS COMMAND LINE INPUTS
# 1. Phenotype
# 2. Tree 
# 3. Genotype
# 4. Reconstruction method
# 5. Multiple test correction method
# 6. Output name
# 7. Alpha value

for (p in 1:length(phenotypes)){
  for (g in 1:length(genotypes)){
    command <- paste("Rscript",
                     run_treewas,
                     paste(path, phenotypes[p], "_pheno.tsv", sep = ""), 
                     paste(path, phenotypes[p], ".tree", sep = ""),
                     paste(path, phenotypes[p], genotypes[g], ".tsv", sep = ""), 
                     recon,
                     test_corr,
                     paste("treewas_", phenotypes[p], genotypes[g], sep = ""),
                     alpha_val, 
                     sep = " ")
    fname <- paste(getwd(), "/", "treewas_", phenotypes[p], genotypes[g], ".pbs", sep = "")
    writeLines(c("#!/bin/sh","####  PBS preamble",
                 paste("#PBS -N treewas_", phenotypes[p], genotypes[g], sep = ""),
                 paste("#PBS -M ", username, sep = ""),  
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
# END SCRIPT -------------------------------------------------------------------
