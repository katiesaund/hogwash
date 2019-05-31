select_test_type <- function(args){
  check_str_is_discrete_or_continuous(args$discrete_or_continuous)

  if(args$discrete_or_continuous == "continuous"){
    run_continuous(args)
  } else {
    run_binary_transition(args)
    run_binary_original(args)
  }

} # end select_test_type()

