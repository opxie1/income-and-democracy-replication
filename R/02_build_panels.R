# ---------------------------------------------------------------------------
# 02_build_panels.R
# Trim each raw panel to the variables the IV tables use, coerce types, attach
# a label to every column, and write one compressed parquet per panel.
# No lags are stored here: they are deterministic transforms that the table
# scripts build with add_lags(), so the parquet holds only source variables.
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

panels_raw <- readRDS(FILE_RAW)

INT_COLS <- c("code_numeric", "year", "year_numeric", "sample",
              "samplebalancefe", "samplebalancegmm", "socialist",
              "noextrapolation", "madid")

build_panel <- function(df, keep) {
  df <- df[, keep, drop = FALSE]
  for (cc in intersect(INT_COLS, names(df))) df[[cc]] <- as.integer(round(df[[cc]]))
  apply_labels(df)
}

panels <- list(
  list(df = build_panel(panels_raw$p5yr,   COLS_5YR),  file = FILE_P5,  name = "5-year"),
  list(df = build_panel(panels_raw$annual, COLS_FREQ), file = FILE_PA,  name = "annual"),
  list(df = build_panel(panels_raw$p10yr,  COLS_FREQ), file = FILE_P10, name = "10-year"),
  list(df = build_panel(panels_raw$p20yr,  COLS_FREQ), file = FILE_P20, name = "20-year"),
  list(df = build_panel(panels_raw$p25yr,  COLS_MAD),  file = FILE_P25, name = "25-year"),
  list(df = build_panel(panels_raw$p50yr,  COLS_MAD),  file = FILE_P50, name = "50-year")
)

for (p in panels) {
  # every kept column must carry a label
  miss <- names(p$df)[vapply(p$df, function(x) is.null(attr(x, "label")), logical(1))]
  if (length(miss)) stop("Missing labels in ", p$name, ": ", paste(miss, collapse = ", "))
  write_parquet(p$df, p$file, compression = "zstd", compression_level = 9)
  cat(sprintf("Wrote %-22s %5d rows x %2d cols  (%.0f KB)\n",
              basename(p$file), nrow(p$df), ncol(p$df),
              file.info(p$file)$size / 1024))
}

cat("Built", length(panels), "labelled panels.\n")
