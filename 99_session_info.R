# ============================================================
# File: 99_session_info.R
# Purpose: Save R session information for reproducibility
# ============================================================

source(file.path("R", "00_setup.R"))

writeLines(
  capture.output(sessionInfo()),
  con = output_table_file("session_info.txt")
)

cat("Session information saved to outputs/tables/session_info.txt\n")
