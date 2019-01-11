# 2018-08-26
# Katie Saund
# 
# Write PBS scripts for all 3 GWAS methods using a standard set of inputs: 
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
# FOR ALL 3: 
path       <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-22_format_data_for_treewas/data/2018-08-23_formatted_data/"
phenotypes <- c("log_cfe", "log_germ_tc", "log_germ_tc_and_gly", "log_growth", "log_sporulation", "log_toxin", "fqR", "severity") # order matters
genotypes  <- c("_snp_stop", "_snp_ns", "_snp_del", "_snp_high", "_gene_stop", "_gene_ns", "_gene_del", "_gene_high", "_pilon_sv", "_roary_pan_genome")

# PHYC & GEE
alpha <- 0.05

# GEE SPECIFIC 
run_GEE <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-26_gee_with_treewas_inputs/lib/2018-08-26_GEE_run.R" 

# PHYC SPECIFIC
run_phyC <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-23_phyC_with_treewas_inputs/lib/2018-08-24_phyC_run.R" 
num_perm <- 10000 # run with X permutations 

# TREEWAS SPECIFIC
run_treewas <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-24_treewas/lib/2018-08-24_run_treewas.R" 
test_corr   <- "fdr" # "bonf" or "fdr"
recon       <- "ML" # "parsimony" or "ML" 

# GEE COMMANDLINE INPUTS: 
# 1. Phenotype
# 2. Tree 
# 3. Genotype
# 4. Output name
# 5. Output directory 
# 6. Alpha
# 7. Annotation for heatmap # if you have no annotation, make this NULL

# PHYC COMMAND LINE INPUTS
# 1. Phenotype
# 2. Tree 
# 3. Genotype
# 4. Output name
# 5. Output directory 
# 6. Number of permutations 
# 7. Alpha
# 8. Annotation for heatmap # if you have no annotation, make this NULL

# TREEWAS COMMAND LINE INPUTS
# 1. Phenotype
# 2. Tree 
# 3. Genotype
# 4. Reconstruction method
# 5. Multiple test correction method
# 6. Output name

for (p in 1:length(phenotypes)){
  for (g in 1:length(genotypes)){
    # GEE segment
    gee_command <- paste("Rscript",
                     run_GEE,
                     paste(path, phenotypes[p], "_treewas_pheno.tsv", sep = ""), 
                     paste(path, phenotypes[p], ".tree", sep = ""),
                     paste(path, phenotypes[p], genotypes[g], ".tsv", sep = ""), 
                     paste(phenotypes[p], genotypes[g], sep = ""),
                     getwd(),
                     alpha,
                     paste(path, phenotypes[p], "_annotation.tsv", sep = ""),
                     sep = " ")
    gee_fname <- paste(getwd(), "/", "gee_", phenotypes[p], genotypes[g], ".pbs", sep = "")
    writeLines(c("#!/bin/sh","####  PBS preamble",
                 paste("#PBS -N GEE_", phenotypes[p], genotypes[g], sep = ""),
                 "#PBS -M katiephd@umich.edu", 
                 "#PBS -m a",
                 "#PBS -l nodes=1:ppn=4,mem=20gb,walltime=50:00:00",
                 "#PBS -V",
                 "#PBS -j oe",
                 "#PBS -A esnitkin_fluxod",
                 "#PBS -q fluxod",
                 "#PBS -l qos=flux",
                 "cd $PBS_O_WORKDIR",
                 gee_command),
               gee_fname,
               sep = "\n")  
    
    # phyC
    phyc_command <- paste("Rscript",
                     run_phyC,
                     paste(path, phenotypes[p], "_treewas_pheno.tsv", sep = ""), 
                     paste(path, phenotypes[p], ".tree", sep = ""),
                     paste(path, phenotypes[p], genotypes[g], ".tsv", sep = ""), 
                     paste(phenotypes[p], genotypes[g], sep = ""),
                     getwd(),
                     num_perm, 
                     alpha,
                     paste(path, phenotypes[p], "_annotation.tsv", sep = ""),
                     sep = " ")
    phyc_fname <- paste(getwd(), "/", "phyc_", phenotypes[p], genotypes[g], ".pbs", sep = "")
    writeLines(c("#!/bin/sh","####  PBS preamble",
                 paste("#PBS -N phyc_", phenotypes[p], genotypes[g], sep = ""),
                 "#PBS -M katiephd@umich.edu", 
                 "#PBS -m a",
                 "#PBS -l nodes=1:ppn=4,mem=20gb,walltime=50:00:00",
                 "#PBS -V",
                 "#PBS -j oe",
                 "#PBS -A esnitkin_fluxod",
                 "#PBS -q fluxod",
                 "#PBS -l qos=flux",
                 "cd $PBS_O_WORKDIR",
                 phyc_command),
               phyc_fname,
               sep = "\n")  
    
    # treewas
    treewas_command <- paste("Rscript",
                     run_treewas,
                     paste(path, phenotypes[p], "_treewas_pheno.tsv", sep = ""), 
                     paste(path, phenotypes[p], ".tree", sep = ""),
                     paste(path, phenotypes[p], genotypes[g], ".tsv", sep = ""), 
                     recon,
                     test_corr,
                     paste(format(Sys.time(), "%Y-%m-%d_"), "treewas_", phenotypes[p], genotypes[g], sep = ""),
                     sep = " ")
    treewas_fname <- paste(getwd(), "/", "treewas_", phenotypes[p], genotypes[g], ".pbs", sep = "")
    writeLines(c("#!/bin/sh","####  PBS preamble",
                 paste("#PBS -N treewas_", phenotypes[p], genotypes[g], sep = ""),
                 "#PBS -M katiephd@umich.edu", 
                 "#PBS -m a",
                 "#PBS -l nodes=1:ppn=4,mem=20gb,walltime=50:00:00",
                 "#PBS -V",
                 "#PBS -j oe",
                 "#PBS -A esnitkin_fluxod",
                 "#PBS -q fluxod",
                 "#PBS -l qos=flux",
                 "cd $PBS_O_WORKDIR",
                 treewas_command),
               treewas_fname,
               sep = "\n")  
    
    # gather all results pbs
    
    
    # make heatmaps of the results pbs
  }
}

# END SCRIPT -------------------------------------------------------------------