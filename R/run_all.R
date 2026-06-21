source(here::here("R", "00_setup.R"))

cat("Running the Income and Democracy IV replication...\n\n")

scripts <- c("01_load_data.R", "02_build_panels.R",
             "03_table2.R", "04_table3.R", "05_table4.R",
             "06_table5.R", "07_table6.R", "08_table7.R",
             "09_verify.R", "10_alternatives.R")
for (s in scripts) {
  cat("====", s, "====\n")
  source(here::here("R", s), local = new.env())
  cat("\n")
}

sink(file.path(PATH_DOCS, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done. Tables are in output/, the parquet panels in data/, the checks in docs/.\n")
if (.Platform$OS.type == "windows") {
  cat("(On Windows, R may print a non-zero exit code after this line from an arrow\n",
      "shutdown bug; every file is already written by this point.)\n", sep = "")
}
