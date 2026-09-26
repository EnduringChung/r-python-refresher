test_that("kmer_counts counts all k-mers", {
  skip("unskip when you implement kmer_counts()")
  expect_equal(
    bioseqkit::kmer_counts("ATATA", 2),
    c(AT = 2L, TA = 2L)
  )
})

test_that("kmer_counts with k=1 is base composition", {
  skip("unskip when you implement kmer_counts()")
  expect_equal(
    bioseqkit::kmer_counts("AAGCT", 1),
    c(A = 2L, C = 1L, G = 1L, T = 1L)
  )
})

test_that("kmer_counts sorts output by name", {
  skip("unskip when you implement kmer_counts()")
  expect_true(!is.unsorted(names(bioseqkit::kmer_counts("ATCGATCG", 3))))
})

test_that("kmer_counts rejects k > nchar(seq)", {
  skip("unskip when you implement kmer_counts()")
  expect_error(bioseqkit::kmer_counts("AT", 3), "k")
})
