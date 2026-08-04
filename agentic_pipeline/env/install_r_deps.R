#!/usr/bin/env Rscript
# Install (or audit) the pinned R dependencies of the PEPPER / DisPo pipeline.
#
#   Rscript agentic_pipeline/env/install_r_deps.R          # audit only
#   Rscript agentic_pipeline/env/install_r_deps.R --install # install what is missing
#
# Exact versions are installed via remotes::install_version so that a figure
# regenerated years from now still uses the code it was validated against.

args <- commandArgs(trailingOnly = TRUE)
do_install <- "--install" %in% args

req_file <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                      "r-requirements.txt")
if (!file.exists(req_file)) stop("r-requirements.txt not found: ", req_file)

lines <- readLines(req_file, warn = FALSE)
lines <- trimws(lines)
lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
parts <- strsplit(lines, "[[:space:]]+")
req <- setNames(vapply(parts, `[`, "", 2L), vapply(parts, `[`, "", 1L))

status <- character(0)
missing <- character(0)
mismatch <- character(0)

for (pkg in names(req)) {
  want <- req[[pkg]]
  have <- tryCatch(as.character(packageVersion(pkg)), error = function(e) NA_character_)
  if (is.na(have)) {
    status <- c(status, sprintf("%-14s %-10s MISSING", pkg, want))
    missing <- c(missing, pkg)
  } else if (identical(have, want)) {
    status <- c(status, sprintf("%-14s %-10s OK", pkg, want))
  } else {
    status <- c(status, sprintf("%-14s %-10s INSTALLED: %s", pkg, want, have))
    mismatch <- c(mismatch, pkg)
  }
}

cat(paste(status, collapse = "\n"), "\n\n", sep = "")

if (length(missing) == 0 && length(mismatch) == 0) {
  cat("All R dependencies are at the expected version.\n")
  quit(status = 0)
}

cat(sprintf("%d absente(s), %d en version differente.\n", length(missing), length(mismatch)))

if (!do_install) {
  cat("Rerun with --install to install the pinned versions.\n")
  # A version mismatch is not fatal for exploration, but the regression tests
  # may legitimately fail on it, so signal it to the caller.
  quit(status = if (length(missing)) 1L else 0L)
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

for (pkg in c(missing, mismatch)) {
  cat("\n>>> installation de ", pkg, " ", req[[pkg]], "\n", sep = "")
  remotes::install_version(pkg, version = req[[pkg]],
                           repos = "https://cloud.r-project.org",
                           upgrade = "never")
}
