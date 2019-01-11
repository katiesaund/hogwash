# 2018-12-14 forked from 2018-11-12
# Katie Saund

# Script will create pbs script to run phyC for all genetic variants, phenotypes, and tree available. 

require(tools)

# PHYC COMMANDLINE INPUTS: 
# 1. Phenotype
# 2. Tree 
# 3. Genotype
# 4. Output name
# 5. Output directory 
# 6. Number of permutations 
# 7. Alpha

path <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-04_format_data_for_gwas/data/2018-09-05_formatted_data_for_gwas/"

run_phyc <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-12-14_phyc_fix_delta_pheno/lib/2018-12-14_phyc_run.R" 

phenotypes <- c("log_cfe", "log_germ_tc", "log_germ_tc_and_gly", "log_growth", "log_sporulation", "log_toxin", "fqR", "severity") 

genotypes <- c("_snp_stop", "_snp_ns", "_snp_del", "_snp_high", "_gene_stop", "_gene_ns", "_gene_del", "_gene_high", "_roary_pan_genome") #_pilon_sv

num_perm <- 10000 # run with X permutations 

alpha <- 0.1

# 1. Phenotype
# 2. Tree 
# 3. Genotype
# 4. Output name
# 5. Output directory 
# 6. Number of permutations 
# 7. Alpha
# 8. Annotation for heatmap # if you have no annotation, make this NULL

for (p in 1:length(phenotypes)){
  for (g in 1:length(genotypes)){
    command <- paste("Rscript",
                     run_phyc,
                     paste(path, phenotypes[p], "_pheno.tsv", sep = ""), 
                     paste(path, phenotypes[p], ".tree", sep = ""),
                     paste(path, phenotypes[p], genotypes[g], ".tsv", sep = ""), 
                     paste(phenotypes[p], genotypes[g], sep = ""),
                     getwd(),
                     num_perm, 
                     alpha,
                     paste(path, phenotypes[p], "_annotation.tsv", sep = ""),
                     sep = " ")
    fname <- paste(getwd(), "/", phenotypes[p], genotypes[g], ".pbs", sep = "")
    writeLines(c("#!/bin/sh","####  PBS preamble",
                 paste("#PBS -N phyc_", phenotypes[p], genotypes[g], sep = ""),
                 "#PBS -M katiephd@umich.edu", 
                "#PBS -m a",
                 "#PBS -l nodes=1:ppn=4,mem=20gb,walltime=100:00:00",
                 "#PBS -V",
                 "#PBS -j oe",
                "#PBS -A esnitkin_fluxod",
                 "#PBS -q fluxod",
                "#PBS -l qos=flux",
                 "cd $PBS_O_WORKDIR",
                 command),
              fname,
              sep = "\n")  
  }
}



