context("Non-plot outputs") #----------------------------------------------#

test_that("save_results_as_r_object error for invalid input", {
  temp_dir <- "."
  temp_name <- "dummy"
  temp_object <- letters[1:3]
  temp_prefix <- "fake_test"
  temp_group_logical <- FALSE
  expect_error(save_results_as_r_object(temp_dir,
                                        temp_name,
                                        temp_object,
                                        temp_prefix,
                                        temp_group_logical))
})

test_that("save_results_as_r_object works valid input", {
  temp_dir <- "."
  temp_name <- "dummy"
  temp_object <- letters[1:3]
  temp_prefix <- "phyc"
  temp_group_logical <- FALSE
  expect_error(save_results_as_r_object(temp_dir,
                                        temp_name,
                                        temp_object,
                                        temp_prefix,
                                        temp_group_logical),
               NA)
})

test_that("save_results_as_r_object works valid input", {
  temp_dir <- "."
  temp_name <- "dummy"
  temp_object <- letters[1:3]
  temp_prefix <- "synchronous"
  temp_group_logical <- FALSE
  expect_error(save_results_as_r_object(temp_dir,
                                        temp_name,
                                        temp_object,
                                        temp_prefix,
                                        temp_group_logical),
               NA)
})

test_that("save_results_as_r_object works valid input", {
  temp_dir <- "."
  temp_name <- "dummy"
  temp_object <- letters[1:3]
  temp_prefix <- "continuous"
  temp_group_logical <- FALSE
  expect_error(save_results_as_r_object(temp_dir,
                                        temp_name,
                                        temp_object,
                                        temp_prefix,
                                        temp_group_logical),
               NA)
})
