select_test_type <- function(args){
  # Function description -------------------------------------------------------
  # Given the type of input to the GWAS program, run either the continuous test
  # or both of the discrete tests.
  #
  # Inputs:
  # args$discrete_or_continuous: Character string. Either "discrete" or
  #                              "continuous".
  #
  # Output:
  # Saved plots and an .Rdata file.
  #
  # Check inputs -------------------------------------------------------------
  check_str_is_discrete_or_continuous(args$discrete_or_continuous)

  # Function -----------------------------------------------------------------
  if (args$discrete_or_continuous == "continuous") {
    run_continuous(args)
  } else {
    run_binary_transition(args)
    run_binary_original(args)
  }

} # end select_test_type()

