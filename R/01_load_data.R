source(here::here("R", "00_setup.R"))

cat("Reading sheets from", basename(FILE_XLS), "...\n")

read_sheet <- function(sheet) {
  suppressMessages(read_excel(FILE_XLS, sheet = sheet))
}

panels_raw <- list(
  p5yr    = read_sheet("5 Year Panel"),
  annual  = read_sheet("Annual Panel"),
  p10yr   = read_sheet("10 Year Panel"),
  p20yr   = read_sheet("20 Year Panel"),
  p25yr   = read_sheet("25 Year Panel"),
  p50yr   = read_sheet("50 Year Panel")
)

for (nm in names(panels_raw)) {
  cat(sprintf("  %-8s %5d rows x %2d cols\n",
              nm, nrow(panels_raw[[nm]]), ncol(panels_raw[[nm]])))
}

saveRDS(panels_raw, FILE_RAW)
cat("Saved raw panels to", basename(FILE_RAW), "\n")
