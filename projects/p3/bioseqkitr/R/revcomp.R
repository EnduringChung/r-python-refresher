# Reverse complement of a DNA string.
#
# @param seq character, length 1, letters A C G T N (case-insensitive)
# @return character, same length, reversed and complemented
# @examples
# rev_comp("ATGC")   # "GCAT"
rev_comp <- function(seq) {
  stopifnot(
    "seq must be a single string" = is.character(seq) && length(seq) == 1
  )
  # TODO (learner): implement — see tests/test-revcomp.R for the contract
  stop("rev_comp() not implemented yet")
}
