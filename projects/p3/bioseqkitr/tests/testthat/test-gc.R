test_that("gc_content returns percentages", {
  expect_equal(bioseqkit::gc_content("ATGC"), 50)
  expect_equal(bioseqkit::gc_content("GGCC"), 100)
  expect_equal(bioseqkit::gc_content("ATAT"), 0)
})

test_that("gc_content is case-insensitive", {
  expect_equal(bioseqkit::gc_content("atgcatgc"), 50)
})

test_that("gc_content guards its inputs", {
  expect_error(bioseqkit::gc_content(c("ATGC", "ATGC")), "single string")
  expect_error(bioseqkit::gc_content(42), "single string")
})

test_that("gc_content handles N runs and empty strings", {
  expect_equal(bioseqkit::gc_content("NNNN"), 0)
  expect_equal(bioseqkit::gc_content(""), 0)
})
