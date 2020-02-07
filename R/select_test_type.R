#' Select which GWAS to run
#'
#' @description Given the type of input to the GWAS program, run either the
#'   continuous test or both of the discrete tests (Synchronous and PhyC).
#' @param args Object with all of the user provided inputs necessary to run the
#'   gwas.
#'
#' @noRd
select_test_type <- function(args){
  # Check inputs ---------------------------------------------------------------
  check_str_is_discrete_or_continuous(args$discrete_or_continuous)

  # Function -------------------------------------------------------------------
  if (args$discrete_or_continuous == "continuous") {
    run_continuous(args)
  } else {
    run_synchronous(args)
    run_phyc(args)
  }
}
