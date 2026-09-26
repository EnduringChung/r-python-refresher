# GC content as a percentage (0-100).
#
# @param seq character, length 1
# @return double in [0, 100], rounded to 2 decimals
# @examples
# gc_content("ATGC")   # 50
# gc_content("atgcatgc")  # 50
gc_content <- function(seq) {
  stopifnot(
    "seq must be a single string" = is.character(seq) && length(seq) == 1
  )
  bases <- strsplit(toupper(seq), "")[[1]]
  if (length(bases) == 0) return(0)
  round(sum(bases %in% c("G", "C")) / length(bases) * 100, 2)
}
