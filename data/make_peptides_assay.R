# Generates the deliberately messy Project 1 dataset.
# Run once: Rscript data/make_peptides_assay.R
set.seed(42)

n <- 120
samples <- sprintf("sample_%03d", 1:(n / 3)) |>
  rep(each = 3)

treatment <- rep(c("Control", "Treated"), length.out = n) |>
  (\(x) sample(x))()

peptides <- rep(c("AngII", "SubstanceP", "Bradykinin", "Oxytocin", "Vasopressin", "Neurotensin"),
                length.out = n)
replicate <- rep(1:3, length.out = n)
batch <- rep(c("B1", "B2", "B3"), length.out = n)

dates <- seq(as.Date("2024-03-01"), by = "day", length.out = 14)
collection_date <- sample(dates, n, replace = TRUE)

concentration <- round(rlnorm(n, meanlog = log(25), sdlog = 0.6), 2)

messy <- tibble::tibble(
  ` Sample ID ` = samples,
  Treatment = treatment,
  Peptide = peptides,
  `Concentration (ng/uL)` = concentration,
  Replicate = replicate,
  Batch = batch,
  `Collection Date` = collection_date
)

# --- mess it up ---
messy <- messy |>
  dplyr::mutate(
    Treatment = dplyr::case_when(
      Treatment == "Control" & replicate == 1 ~ "control",
      Treatment == "Control" & replicate == 2 ~ "CTRL",
      Treatment == "Treated" & replicate == 2 ~ "treated",
      .default = Treatment
    ),
    `Concentration (ng/uL)` = dplyr::case_when(
      runif(n) < 0.05 ~ -999,
      .default = `Concentration (ng/uL)`
    )
  )

# mixed date formats
messy$`Collection Date` <- as.character(messy$`Collection Date`)
idx2 <- seq(2, n, by = 3)
messy$`Collection Date`[idx2] <- format(as.Date(messy$`Collection Date`[idx2]), "%m/%d/%Y")
idx3 <- seq(3, n, by = 3)
messy$`Collection Date`[idx3] <- format(as.Date(messy$`Collection Date`[idx3]), "%B %d, %Y")

# string NAs in some concentrations
na_idx <- sample(n, 8)
messy$`Concentration (ng/uL)`[na_idx[1:4]] <- "n/a"
messy$`Concentration (ng/uL)`[na_idx[5:8]] <- ""

# whitespace + case chaos in sample ids
messy$` Sample ID ` <- ifelse(runif(n) < 0.1, paste0(" ", toupper(messy$` Sample ID `), " "),
                              messy$` Sample ID `)

# duplicate a few rows entirely
dupes <- messy[sample(n, 5), ]
messy <- rbind(messy, dupes)

messy <- messy[sample(nrow(messy)), ]   # shuffle row order

readr::write_csv(messy, "data/peptides_assay.csv")
cat("wrote", nrow(messy), "messy rows to data/peptides_assay.csv\n")
