# 2019-03-05
# Katie Saund
load("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/2019-01-23_snpmat_parsed.RData")
str(parsed)

create_snp_subset_gene_lookup_mat <- function(subset_logical, name, parsed_data, output_dir){
  lookup <- matrix("", nrow = sum(subset_logical), ncol = 2)
  colnames(lookup) <- c(name, "gene")
  lookup[ , 1] <- row.names(parsed_data$mat)[subset_logical]
  lookup[ , 2] <- parsed_data$genes[subset_logical]
  write.table(x = lookup, file = paste(output_dir, format(Sys.time(), "%Y-%m-%d_"), name, "_gene_lookup.tsv", sep = ""), sep = "\t")
}

dir <- "/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2019-01-22_parse_code_snpmat/data/"

create_snp_subset_gene_lookup_mat(parsed$stop, "stop_snp", parsed, dir)
create_snp_subset_gene_lookup_mat(parsed$snpeff_high, "snpeff_high_snp", parsed, dir)
create_snp_subset_gene_lookup_mat(parsed$ns_mut, "nonsynonymous_snp", parsed, dir)
