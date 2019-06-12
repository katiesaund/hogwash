save_results_as_r_object <- function(dir, name, object){
  # Function description -------------------------------------------------------
  # Save all of the non-plot outputs in a .rda file.
  #
  # Inputs:
  # dir. String. Path to directory where file should be saved.
  # name. String. Name of file without suffix.
  # object. List of various inputs to save.
  #
  # Output:
  # A .rda file.
  #
  # Check inputs ---------------------------------------------------------------
  check_is_string(dir)
  check_if_dir_exists(dir)
  check_is_string(name)
  # TODO check on object?

  # Function & output ----------------------------------------------------------
  save(object, file = paste0(dir, "/phyc_", name, ".rda"))
} # end save_results_as_r_object()

# End of script ----------------------------------------------------------------
