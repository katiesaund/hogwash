
save_continuous <- function(hogwash_continuous, file_name){
  save(hogwash_continuous, file = file_name)
}
save_synchronous <- function(hogwash_synchronous, file_name){
  save(hogwash_synchronous, file = file_name)
}
save_convergence <- function(hogwash_phyc, file_name){
  save(hogwash_phyc, file = file_name)
}

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
    fname <- paste0(dir, "/hogwash_", prefix, "_grouped_", name, ".rda")
    # save(object, file = paste0(dir, "/hogwash_", prefix, "_grouped_", name, ".rda"))
  } else {
    fname <- paste0(dir, "/hogwash_", prefix, "_", name, ".rda")
    # save(object, file = paste0(dir, "/hogwash_", prefix, "_", name, ".rda"))
  }

  if (prefix == "continuous") {
    save_continuous(object, fname)
  } else if (prefix == "synchronous") {
    save_synchronous(object, fname)
  } else if (prefix == "convergence") {
    save_convergence(object, fname)
  } else {
    stop("test name incorrect")
  }

} # end save_results_as_r_object()


# End of script ----------------------------------------------------------------
