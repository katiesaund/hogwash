library(microbeGWAS)
context("Non-plot outputs") #----------------------------------------------#

# test save_data_table ---------------------------------------------------------

# save_data_table <- function(matrix, output_dir, pheno_geno_name, extension){
#   # Function description ------------------------------------------------------#
#   # Save a table of information to file.
#   #
#   # Inputs:
#   # TODO
#   # matrix. ?
#   # output_dir. String.
#   # pheno_geno_name. String.
#   # extension. String.
#   #
#   # Output:
#   # A data file.
#   #
#   # Check inputs ---------------------------------------------------------------
#   check_is_string(output_dir)
#   check_if_dir_exists(output_dir)
#   check_is_string(pheno_geno_name)
#   check_is_string(extension)
#
#   # Function/Output ------------------------------------------------------------
#   write.table(matrix,
#               file = paste0(output_dir, "/phyc_", pheno_geno_name, extension),
#               sep = "\t",
#               quote = FALSE,
#               row.names = TRUE,
#               col.names = TRUE)
# } # end save_data_table()
#


# test save_results_as_r_object ------------------------------------------------

# save_results_as_r_object <- function(dir, name, object){
#   # Function description ------------------------------------------------------#
#   # Save a data object as a .rda file.
#   #
#   # Inputs:
#   # TODO
#   # dir. String.
#   # name. String.
#   # object. ?
#   #
#   # Output:
#   # A .rda file.
#   #
#   # Check inputs ---------------------------------------------------------------
#   check_is_string(dir)
#   check_if_dir_exists(dir)
#   check_is_string(name)
#   # TODO check on object?
#
#   # Function/Output ------------------------------------------------------------
#   save(object, file = paste0(dir, "/phyc_", name, ".rda"))
# } # end save_results_as_r_object()
#
#

# test save_hits ---------------------------------------------------------------


# save_hits <- function(hits, output_dir, output_name, pval_name){
#   # Function description -------------------------------------------------------
#   # Create a file name and save results to that file name.
#   #
#   # Inputs:
#   # hits.        Vector. Pvals.
#   # output_dir.  Character.
#   # output_name. Character.
#   # pval_name.   Character.
#   #
#   # Output:
#   # None
#   #
#   # Check inputs ---------------------------------------------------------------
#   check_is_string(pval_name)
#   check_is_string(output_name)
#   check_if_dir_exists(output_dir)
#
#   # Function/Output ------------------------------------------------------------
#   if (nrow(hits) > 0){
#     file_name <- create_file_name(output_dir, output_name, pval_name)
#     hit_and_rank <- cbind(hits, rank(hits))
#     colnames(hit_and_rank)[2] <- "rank"
#     write.table(x = hit_and_rank, file = paste(file_name, ".tsv", sep = ""), sep = "\t", quote = FALSE, eol = "\n", row.names = TRUE, col.names = TRUE)
#   } else { # added 2018-11-13
#     file_name <- create_file_name(output_dir, output_name, pval_name)
#     empty <- matrix(NA, nrow = 1, ncol = 1)
#     colnames(empty) <- "no_sig_hits"
#     row.names(empty) <- "no_sig_hits"
#     write.table(x = empty, file = paste(file_name, ".tsv", sep = ""), sep = "\t", quote = FALSE, eol = "\n", row.names = TRUE, col.names = TRUE)
#   }
# } #end save_hits()


# END SCRIPT -------------------------------------------------------------------
