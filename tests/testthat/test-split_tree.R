test_that("returns list of trees", {
  expect_equal(divide_tree_into_subtrees(tree)[[1]], tree)
})
