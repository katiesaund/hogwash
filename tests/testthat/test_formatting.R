context("Formatting") #--------------------------------------------------------#

# Test format tree
test_that("format_tree gives error for invalid input", {
  set.seed(10)
  temp_tree <- ape::rtree(20)
  expect_error(format_tree(temp_tree))
})

test_that("format_tree works for valid input", {
  set.seed(10)
  temp_tree <- ape::rtree(20)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  expect_error(format_tree(temp_tree), NA)
})

test_that("format_tree works for valid input", {
  set.seed(10)
  temp_tree <- ape::rtree(20)
  temp_tree$node.label <- c("", rep(100, ape::Nnode(temp_tree) - 1))
  expect_error(format_tree(temp_tree), NA)
})

test_that("format_tree works for valid input", {
    set.seed(10)
    temp_tree <- ape::rtree(20)
    temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
    expect_error(format_tree(temp_tree), NA)
})

test_that("format_tree roots unrooted tree", {
    set.seed(10)
    temp_tree <- ape::rtree(20)
    temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
    temp_tree <- ape::unroot(temp_tree)
    expect_error(format_tree(temp_tree), NA)
})
