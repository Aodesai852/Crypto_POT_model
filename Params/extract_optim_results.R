rm(list = ls())
# ---------- user settings ----------
setwd("D:/MyFiles/EVT/Crypto_POT")
# =========================
# Inputs
# =========================
input_dir_down <- "./result_prior/down_tail/optim_result"
input_dir_up   <- "./result_prior/up_tail/optim_result"

# =========================
# Outputs
# =========================
output_file_down <- "./code/Params/est_param_down.R"
output_file_up   <- "./code/Params/est_param_up.R"

# =========================
# Helper: numeric vector -> "c(a, b, ...)" string
# =========================
format_c <- function(x, digits = 15) {
  x <- as.numeric(x)
  s <- formatC(x, digits = digits, format = "fg", flag = "#")
  s <- sub("\\.$", "", s)
  paste0("c(", paste(s, collapse = ", "), ")")
}

# =========================
# Read one .RData file and extract c(est_para[[2]], est_para[[3]])
# Assumption: each .RData contains exactly ONE object (e.g., neg_optim_result)
# =========================
extract_est <- function(rdata_file) {
  e <- new.env(parent = emptyenv())
  obj_names <- load(rdata_file, envir = e)
  
  if (length(obj_names) != 1) {
    stop(sprintf("File %s: expected exactly 1 object, got %d.",
                 basename(rdata_file), length(obj_names)))
  }
  
  obj <- e[[obj_names]]
  
  if (!is.list(obj) || is.null(obj[["est_para"]])) {
    stop(sprintf("File %s: object has no element named 'est_para'.",
                 basename(rdata_file)))
  }
  
  est_para <- obj[["est_para"]]
  if (!is.list(est_para) || length(est_para) < 3) {
    stop(sprintf("File %s: 'est_para' is not a list with length >= 3.",
                 basename(rdata_file)))
  }
  
  c(est_para[[2]], est_para[[3]])
}

# =========================
# Process one folder -> write one output .R file
# =========================
process_folder <- function(input_dir, output_file) {
  if (!dir.exists(input_dir)) stop("Folder not found: ", input_dir)
  
  files <- list.files(input_dir, pattern = "\\.RData$", full.names = TRUE)
  if (length(files) == 0) stop("No .RData files found in: ", input_dir)
  
  # Names like "AAVE" from "AAVE.RData"
  keys <- tools::file_path_sans_ext(basename(files))
  
  # Extract
  vals <- lapply(files, extract_est)
  names(vals) <- keys
  
  # Write output script
  lines <- character(0)
  lines <- c(lines, "est.param <- list(")
  
  for (i in seq_along(vals)) {
    nm <- names(vals)[i]
    vec_str <- format_c(vals[[i]])
    comma <- if (i < length(vals)) "," else ""
    lines <- c(lines, sprintf("  %s = %s%s", nm, vec_str, comma))
  }
  
  lines <- c(lines, ")")
  writeLines(lines, con = output_file)
  
  message("Wrote: ", output_file)
}

# =========================
# Run for both folders
# =========================
process_folder(input_dir_down, output_file_down)
process_folder(input_dir_up,   output_file_up)