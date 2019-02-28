library(microbeGWAS)

context("first test")
test_that("check_tree_is_valid works", {
  for (i in 2:10){
    expect_true(check_tree_is_valid(rtree(i)))
  }
  invalid_tree <- rtree(10)
  for (j in 1:Nedge(invalid_tree)){
    if (invalid_tree$edge[j, 2] == 1){
      invalid_tree$edge[j, 2] <- 10000
      break
    }
  }

  expect_that(check_tree_is_valid(invalid_tree), throws_error())
})
