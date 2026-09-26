test_that("rev_comp reverses and complements", {
  skip("unskip when you implement rev_comp()")
  expect_equal(bioseqkit::rev_comp("ATGC"), "GCAT")
  expect_equal(bioseqkit::rev_comp("A"), "T")
  expect_equal(bioseqkit::rev_comp("AAAAAC"), "GTTTTT")
})

test_that("rev_comp is case-insensitive, output uppercase", {
  skip("unskip when you implement rev_comp()")
  expect_equal(bioseqkit::rev_comp("atgc"), "GCAT")
})

test_that("rev_comp handles N", {
  skip("unskip when you implement rev_comp()")
  expect_equal(bioseqkit::rev_comp("ATGN"), "NCAT")
})

test_that("rev_comp guards its inputs", {
  skip("unskip when you implement rev_comp()")
  expect_error(bioseqkit::rev_comp(c("AT", "GC")), "single string")
})
