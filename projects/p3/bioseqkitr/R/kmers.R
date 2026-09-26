# Counts of all k-mers (substrings of length k) in a sequence.
#
# @param seq character, length 1
# @param k integer, 1 <= k <= nchar(seq)
# @return named integer vector, names are the k-mers, sorted by name
# @examples
# kmer_counts("ATATA", 2)   # AT=2 TA=2
kmer_counts <- function(seq, k) {
  stopifnot(
    "seq must be a single string" = is.character(seq) && length(seq) == 1,
    "k must be a positive integer" = is.numeric(k) && k >= 1 && k == round(k)
  )
  # TODO (learner): implement — see tests/test-kmers.R for the contract
  stop("kmer_counts() not implemented yet")
}
