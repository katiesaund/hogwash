#' save_continuous
#'
#' @description Save the results of the Continuous Test in an rdata object
#'  called hogwash_continuous.
#' @param hogwash_continuous List of objects.
#' @param file_name Character.
#'
#' @noRd
save_continuous <- function(hogwash_continuous, file_name){
  save(hogwash_continuous, file = file_name)
} # end save_continuous()

#' save_synchronous
#'
#' @description Save the results of the Synchronous Test in an rdata object
#'  called hogwash_synchronous
#'
#' @param hogwash_synchronous List of objects.
#' @param file_name Character.
#'
#' @noRd
save_synchronous <- function(hogwash_synchronous, file_name){
  save(hogwash_synchronous, file = file_name)
} # end save_synchronous()

#' save_phyc
#'
#' @description Save the results of the PhyC Test in an rdata object
#'  called hogwash_synchronous
#'
#' @param hogwash_phyc List of objects.
#' @param file_name Character.
#'
#' @noRd
save_phyc <- function(hogwash_phyc, file_name){
  save(hogwash_phyc, file = file_name)
} # end save_phyc()

#' save_results_as_r_object
#'
#' @description Save all of the non-plot outputs in a .rda file.
#'
#' @param dir String. Path to directory where file should be saved.
#' @param name String. Name of file without suffix.
#' @param object  List of various inputs to save.
#' @param prefix Character. Test type (continuous, synchronous, convergence)
#' @param group_logical  Logical. Whether or not genotypes were grouped.
#'
#' @noRd
save_results_as_r_object <- function(dir, name, object, prefix, group_logical){
  # Check inputs ---------------------------------------------------------------
  check_is_string(dir)
  check_if_dir_exists(dir)
  check_is_string(name)
  check_is_string(prefix)

  # Function & output ----------------------------------------------------------
  if (group_logical) {
    fname <- paste0(dir, "/hogwash_", prefix, "_grouped_", name, ".rda")
  } else {
    fname <- paste0(dir, "/hogwash_", prefix, "_", name, ".rda")
  }

  if (prefix == "continuous") {
    save_continuous(object, fname)
  } else if (prefix == "synchronous") {
    save_synchronous(object, fname)
  } else if (prefix == "phyc") {
    save_phyc(object, fname)
  } else {
    stop("test name incorrect")
  }
} # end save_results_as_r_object()
