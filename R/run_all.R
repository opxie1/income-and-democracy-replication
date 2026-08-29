source(here::here("R", "00_setup.R"))

cat("Running the Income and Democracy IV replication...\n\n")

# order
scripts <- c("01_load_data.R", "02_build_panels.R",
             "03_table2.R", "04_table3.R", "05_table4.R",
             "06_table5.R", "07_table6.R", "08_table7.R",
             "09_verify.R", "10_alternatives.R", "12_instruments.R",
             "13_aggregation.R", "14_weakiv.R", "16_montecarlo.R")
for (s in scripts) {
  cat("====", s, "====\n")
  source(here::here("R", s), local = new.env())
  cat("\n")
}

# session
sink(file.path(PATH_DOCS, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done. Tables are in output/, the parquet panels in data/, the checks in docs/.\n")
if (.Platform$OS.type == "windows") {
  cat("On Windows, R sometimes prints a non-zero exit code after this line.\n",
      "It comes from the arrow package as the package closes.\n",
      "Every file is already written at this point.\n", sep = "")
}
