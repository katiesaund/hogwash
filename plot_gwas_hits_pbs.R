# 2018-08-29
# Katie Saund
# 
# Write PBS script to plot heatmaps of the 3 method intersection hits. 
# 

results_dir <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-27_run_all_gwas/data/2018-08-27_all_pheno_no_dates_for_gather_results/compare_methods/"
gwas_format_path <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-22_format_data_for_treewas/data/2018-08-23_formatted_data/"
gwas_intersection_results <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-27_run_all_gwas/data/2018-08-27_all_pheno_no_dates_for_gather_results/compare_methods/all_methods_intersection_per_test_sig_hits.rda"
make_heatmaps <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-27_run_all_gwas/lib/2018-08-29_plot_GWAS_hits.R"

command <- paste("Rscript",
                 make_heatmaps,
                 results_dir, 
                 gwas_format_path,
                 gwas_intersection_results, 
                 sep = " ")
fname <- paste(getwd(), "/", "plot_GWAS_hits.pbs", sep = "")
writeLines(c("#!/bin/sh","####  PBS preamble",
             "#PBS -N plot_GWAS_hits",
             "#PBS -M katiephd@umich.edu", 
             "#PBS -m abe",
             "#PBS -l nodes=1:ppn=4,mem=20gb,walltime=50:00:00",
             "#PBS -V",
             "#PBS -j oe",
             "#PBS -A esnitkin_fluxod",
             "#PBS -q fluxod",
             "#PBS -l qos=flux",
             "cd $PBS_O_WORKDIR",
             command),
           fname,
           sep = "\n")  


