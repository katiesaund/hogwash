save_results_as_r_object <- function(dir, name, object, prefix, group_logical){
  # Function description -------------------------------------------------------
  # Save all of the non-plot outputs in a .rda file.
  #
  # Inputs:
  # dir. String. Path to directory where file should be saved.
  # name. String. Name of file without suffix.
  # object. List of various inputs to save.
  # prefix. Character. Test type (continuous, synchronous, convergence)
  # group_logical. Logical. Whether or not genotypes were grouped.
  #
  # Output:
  # .rda file.
  #
  # Check inputs ---------------------------------------------------------------
  check_is_string(dir)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_is_string(prefix)
  # TODO check on object?

  # Function & output ----------------------------------------------------------
  if (group_logical) {
    save(object, file = paste0(dir, "/hogwash_", prefix, "_grouped_", name, ".rda"))
  } else {
    save(object, file = paste0(dir, "/hogwash_", prefix, "_", name, ".rda"))
  }

} # end save_results_as_r_object()

# End of script ----------------------------------------------------------------
