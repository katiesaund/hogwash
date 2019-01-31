# Katie Saund
# 2018-10-02
# Update summary figures

summary_count <- read.csv("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/summary_hits_per_method.tsv", 
                          sep = "\t")

summary_count[summary_count == 0] <- 1
summary_count <- log(summary_count)
print(summary_count)

count_per_method <- read.csv("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/hits_per_method.tsv", sep = "\t")
count_per_method[count_per_method == 0] <- 1
count_per_method <- log(count_per_method)
print(count_per_method)



white_to_black = colorRampPalette(c("white", "black")) # Set up presence/absence colors
num_color = length(min(summary_count):max(summary_count)) # We want a unique shade of grayscale for each level (should be 2 for presence/absence)
htmp_col = white_to_black(num_color) # Generate a character vector of colors with each shade of grayscale
htmp <- ComplexHeatmap::Heatmap(
  matrix = summary_count, 
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  row_names_side = "left", 
  col = htmp_col, 
  show_column_names = TRUE, 
  show_heatmap_legend = TRUE, 
  row_title = "GWAS methods",
  column_title = "Log(Number of significant hits)", column_names_gp = gpar(cex = 0.6))

draw(htmp)
pdf("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/2018-10-02_log_summary_hits_per_method_heatmap.pdf")
draw(htmp)
dev.off()

write.table(summary_count, file = "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/2018-10-02_log_summary_hits_per_method.tsv", sep = "\t")



white_to_black = colorRampPalette(c("white", "black")) # Set up presence/absence colors
num_color = length(min(count_per_method):max(count_per_method)) # We want a unique shade of grayscale for each level (should be 2 for presence/absence)
htmp_col = white_to_black(num_color) # Generate a character vector of colors with each shade of grayscale
htmp <- ComplexHeatmap::Heatmap(
  matrix = count_per_method, 
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  row_names_side = "left", 
  col = htmp_col, 
  show_column_names = TRUE, 
  show_heatmap_legend = TRUE, 
  row_title = "GWAS methods", 
  column_title = "Log(Number of significant hits)", column_names_gp = gpar(cex = 0.6))

draw(htmp)
pdf("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/2018-10-02_log_hits_per_method_heatmap.pdf")
draw(htmp)
dev.off()

write.table(count_per_method, file = "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-09-05_run_all_gwas_for_LAMGC_and_U01_conferences/data/2018-09-06_all_results/2018-10-02_log_hits_per_method.tsv", sep = "\t")


# END OF SCRIPT ----------------------------------------------------------------
