rm(list=ls())
setwd("D:/MyFiles/EVT/BTC_POT/code")


# ---- return descriptive statistics ----

path_desc <- "./result/ret_desc_stats"
files <- list.files(path_desc,
                    pattern = "\\.RData$",
                    full.names = TRUE)

desc_get_one_row <- function(f) {
  # load into a temporary environment to avoid name clashes
  env <- new.env()
  load(f, envir = env)   # object name: desc.stats
  
  temp <- env$ret.desc.stats
  
  coin <- tools::file_path_sans_ext(basename(f))
  
  # function to generate stars from p-value
  stars_fun <- function(p) {
    if (is.na(p)) return("")
    if (p <= 0.01) return("\\text{***}")
    if (p <= 0.05) return("\\text{***}")
    if (p <= 0.10) return("\\text{***}")
    ""
  }
  
  # build JB and ADF strings: stat(pvalue)***
  jb_star  <- stars_fun(temp$JB_pvalue)
  adf_star <- stars_fun(temp$adf_pvalue)
  
  JB_test <- paste0(
    gsub(
      "([0-9.]+)e[+]?(0*)([0-9]+)",
      "\\1\\\\times 10^{\\3}",
      formatC(temp$JB_stat, format = "e", digits = 2)
    ),
    "(",
    sprintf("%.2f", temp$JB_pvalue),
    ")",
    jb_star
  )
  
  ADF_test <- sprintf("%.2f(%.2f)%s",
                      temp$adf_stat, temp$adf_pvalue, adf_star)
  
  data.frame(
    coin      = coin,
    Mean    = temp$mean,
    Med     = temp$median,
    Min     = temp$min,
    Max     = temp$max,
    Std_dev = temp$std_dev,
    Skew    = temp$skew,
    Kurt    = temp$kurt,
    JB_test = JB_test,
    ADF_test = ADF_test,
    stringsAsFactors = FALSE
  )
}

tab_raw <- do.call(rbind, lapply(files, desc_get_one_row))

## (optional) sort rows by coin name
tab_raw <- tab_raw[order(tab_raw$coin), ]

## Round numeric columns to 2 decimals for printing
num_cols <- c("Mean", "Med", "Min", "Max", "Std_dev", "Skew", "Kurt")
tab_print <- tab_raw
tab_print[num_cols] <- lapply(tab_print[num_cols],
                              function(x) sprintf("%.2f", x))

# make Latex table
desc_make_latex_table <- function(df,
                             caption = "Descriptive statistics of cryptocurrency returns",
                             label   = "tab:coin_descriptive") {
  
  # Convert everything to character first
  df_chr <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
  
  header <- paste0(
    "\\begin{table}[htbp]
\\centering
\\caption{", caption, "}
\\label{", label, "}
\\begin{adjustbox}{width=\\textwidth}
\\begin{tabular}{lccccccccc}
\\toprule
Coin & Mean & Med & Min & Max & Std. Dev. & Skew. & Kurt. & JB Test & ADF Test \\\\
\\midrule"
  )
  
  # Build rows:
  # column 1 = FX name (no $)
  # columns 2–10 = wrap entire cell in $ ... $
  rows <- apply(df_chr, 1, function(r) {
    r <- trimws(r)                     # remove leading/trailing spaces
    r[-1] <- paste0("$", r[-1], "$")   # wrap all numeric/stat fields
    paste(r, collapse = " & ")
  })
  
  rows <- paste0(rows, " \\\\")
  
  footer <-
    "\\bottomrule
\\end{tabular}
\\end{adjustbox}
\\begin{flushleft}
{\\footnotesize
\\textit{Notes:} This table shows daily retrun statistics for different cryptocurrencies. Descriptive statistics are expressed as percentages. The sample period spans from 2020-10-18 to 2026-02-03. The ADF test (Dickey \\& Fuller, 1979, 1981) is used to examine stationarity; a statistically significant test statistic indicates rejection of the null hypothesis of a unit root, implying that the series is stationary. The Jarque--Bera test is used to assess normality; a statistically significant result indicates rejection of the null hypothesis of normality, suggesting that the distribution exhibits non-normal skewness and/or kurtosis. ***, **, * indicate significance at the 1\\%, 5\\%, and 10\\% levels.}
\\end{flushleft}
\\end{table}"
  
  cat(header, "\n",
      paste(rows, collapse = "\n"), "\n",
      footer, sep = "")
}


sink("./table_txt/desc_table.txt")   # start writing file
desc_make_latex_table(tab_print)
sink()     # close


# ---- optim result ----

path_up_optim_result <- "./result/up_tail/optim_result"
path_down_optim_result <- "./result/down_tail/optim_result"

latex_sci <- function(x, digits = 2) {
  if (is.na(x)) return("")
  
  sci <- formatC(x, format = "e", digits = digits)
  parts <- strsplit(sci, "e")[[1]]
  base <- parts[1]
  expo <- as.integer(parts[2])
  
  if (expo == 0) {
    return(base)  # e.g. "0.97"
  } else {
    return(paste0(base, "\\times 10^{", expo, "}"))
  }
}


# for test: file = "./result/up_tail/optim_result/BTC.RData"
# this helper builds the 2 LaTeX lines for one coin
optim_result_build_row_from_file <- function(file) {
  env <- new.env()
  load(file, envir = env)          # object name: pos_optim_result or neg_optim_result
  
  if (exists("pos_optim_result", envir = env)) {
    env$optim_result <- env$pos_optim_result
    rm("pos_optim_result", envir = env)
    
  } else if (exists("neg_optim_result", envir = env)) {
    env$optim_result <- env$neg_optim_result
    rm("neg_optim_result", envir = env)
  }
  
  
  obj <- env$optim_result
  
  # FX name: USD_GBP -> USD/GBP
  coin <- tools::file_path_sans_ext(basename(file))
  
  # parameters
  param <- obj$est_para
  phi  <- as.numeric(param[[2]])   # phi_0, phi_1, phi_2
  psi  <- as.numeric(param[[3]])   # psi_0, psi_1, psi_2
  
  # standard errors and p-values
  se    <- as.numeric(obj$stderr)        # length 6
  pvals <- as.numeric(obj$p_values)      # length 6
  
  # log-likelihood
  loglik <- -as.numeric(obj$neg_loglik)
  
  # correlations
  corr  <- as.numeric(obj$corr)
  rho1  <- corr[1]
  rho2  <- corr[2]
  
  # stars
  star_fun <- function(p) {
    if (is.na(p)) return("")
    if (p <= 0.01) return("\\text{***}")
    if (p <= 0.05) return("\\text{***}")
    if (p <= 0.10) return("\\text{***}")
    ""
  }
  
  # format coefficient + star and std.error
  fmt_coef <- function(beta, p, digits = 2) {
    if (is.na(beta)) return("")
    inside <- paste0(latex_sci(beta, digits), star_fun(p))
    paste0("$", inside, "$")
  }
  
  fmt_se <- function(s, digits = 2) {
    if (is.na(s)) return("")
    inside <- latex_sci(s, digits)
    paste0("$(", inside, ")$")
  }
  
  # six parameters in order: phi0,phi1,phi2,psi0,psi1,psi2
  coef_vec <- c(
    fmt_coef(phi[1], pvals[1]),
    fmt_coef(phi[2], pvals[2]),
    fmt_coef(phi[3], pvals[3]),
    fmt_coef(psi[1], pvals[4]),
    fmt_coef(psi[2], pvals[5]),
    fmt_coef(psi[3], pvals[6])
  )
  se_vec <- c(
    fmt_se(se[1]),
    fmt_se(se[2]),
    fmt_se(se[3]),
    fmt_se(se[4]),
    fmt_se(se[5]),
    fmt_se(se[6])
  )
  
  fmt_plain <- function(x, digits = 3, sci = FALSE) {
    if (is.na(x)) return("")
    inside <- if (sci) latex_sci(x, digits) else sprintf(paste0("%.", digits, "f"), x)
    paste0("$", inside, "$")
  }
  
  # first line: coin + coefficients + loglik + correlations
  line1 <- paste(
    c(
      coin,
      coef_vec,
      fmt_plain(loglik, digits = 2, sci = TRUE),  # or sci = FALSE
      fmt_plain(rho1,   digits = 3, sci = FALSE),
      fmt_plain(rho2,   digits = 3, sci = FALSE)
    ),
    collapse = " & "
  )
  
  
  # second line: empty country + std.errors + blanks for last 3 columns
  line2 <- paste(
    c("",
      se_vec,
      "", "", ""),
    collapse = " & "
  )
  
  c(
    paste0(line1, " \\\\"),
    paste0(line2, " \\\\")
  )
}


optim_result_make_latex_table <- function(
    folder,
    caption = "Estimation results of the uPoT model",
    label   = "tab:tail_uPot"
) {
  files <- list.files(folder, pattern = "\\.RData$", full.names = TRUE)
  files <- sort(files)  # optional: sort by name
  
  # build LaTeX lines for each file (two lines per asset)
  row_lines <- unlist(lapply(files, optim_result_build_row_from_file))
  
  header <- paste0(
    "\\begin{landscape}
\\begin{longtable}{lccccccccc}
\\caption{", caption, "}\\\\
\\label{", label, "} \\\\
\\toprule
& \\multicolumn{3}{c}{Shape parameter (tail index) $\\xi_t$}
& \\multicolumn{3}{c}{Scale parameter $\\alpha_t$}
& & & \\\\
\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}
Coin
& $\\phi_0$ & $\\phi_1$ & $\\phi_2$
& $\\psi_0$ & $\\psi_1$ & $\\psi_2$
& $\\log L$
& $\\rho(\\hat{\\alpha}_{it},\\tilde{\\alpha}_{it})$
& $\\rho(\\Delta\\hat{\\alpha}_{it},\\Delta\\tilde{\\alpha}_{it})$ \\\\
\\midrule
\\endfirsthead

\\multicolumn{10}{r}{\\textit{(continued from previous page)}}\\\\
\\toprule
& \\multicolumn{3}{c}{Shape parameter (tail index) $\\xi_t$}
& \\multicolumn{3}{c}{Scale parameter $\\alpha_t$}
& & & \\\\
\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}
Coin
& $\\phi_0$ & $\\phi_1$ & $\\phi_2$
& $\\psi_0$ & $\\psi_1$ & $\\psi_2$
& $\\log L$
& $\\rho(\\hat{\\alpha}_{it},\\tilde{\\alpha}_{it})$
& $\\rho(\\Delta\\hat{\\alpha}_{it},\\Delta\\tilde{\\alpha}_{it})$ \\\\
\\midrule
\\endhead

\\midrule
\\multicolumn{10}{r}{\\textit{(continued on next page)}}\\\\
\\endfoot

\\bottomrule
\\endlastfoot
"
  )
  
  
  footer <- 
    "\\end{longtable}

\\begin{flushleft}
{\\footnotesize
\\textit{Notes:} The sample period spans from 2020-10-18 to 2026-02-03. Standard errors are reported in parentheses. $\\log L$ denotes the maximized log-likelihood value. $\\rho(\\hat{\\alpha}_{it},\\tilde{\\alpha}_{it})$ and $\\rho(\\Delta\\hat{\\alpha}_{it},\\Delta\\tilde{\\alpha}_{it})$ are Pearson correlation coefficients between the estimated risk measures and their benchmarks. The key hyperparameter of the uPOT model is the threshold, which is set at the 95\\% quantile. As part of robustness checks, we also estimate the uPOT model using the 85\\% and 90\\% quantiles as alternative thresholds. The results indicate that the estimates are robust to the choice of the threshold, and the corresponding tables are available from the authors upon request. $^{***}$, $^{**}$, and $^{*}$ indicate significance at the 1\\%, 5\\%, and 10\\% levels, respectively.
}
\\end{flushleft}

\\end{landscape}
"


cat(header, paste(row_lines, collapse = "\n"), "\n", footer, sep = "")
}

sink("./table_txt/up_optim_result_table.txt")
optim_result_make_latex_table(path_up_optim_result)
sink()

sink("./table_txt/down_optim_result_table.txt")
optim_result_make_latex_table(path_down_optim_result)
sink()

# ---- the DY table ----

dca_table_to_latex <- function(dca,
                               digits  = 2,
                               caption = "Volatility spillover table.",
                               label   = "tab:spillover",
                               from_col_name = "FROM",
                               to_row_name   = "TO",
                               inc_row_name  = "Inc.Own",
                               net_row_name  = "NET") {
  if (is.null(dca$TABLE)) stop("dca$TABLE is NULL.")
  
  tbl_raw <- dca$TABLE
  if (is.data.frame(tbl_raw)) tbl_raw <- as.matrix(tbl_raw)
  
  rn <- rownames(tbl_raw)
  if (is.null(rn)) stop("dca$TABLE must have rownames.")
  cn <- colnames(tbl_raw)
  if (is.null(cn)) stop("dca$TABLE must have colnames.")
  
  # --- helpers ---
  clean_cell <- function(x) {
    x <- trimws(x)
    x <- gsub('"', "", x, fixed = TRUE)
    x
  }
  
  latex_escape <- function(x) {
    x <- gsub("\\\\", "\\\\textbackslash{}", x, perl = TRUE)
    x <- gsub("([_%&#$])", "\\\\\\1", x, perl = TRUE)
    x
  }
  
  # Main-body numbers
  fmt_num <- function(x) {
    ifelse(
      is.na(x),
      "",
      paste0("$", formatC(x, format = "f", digits = digits), "$")
    )
  }
  
  # Bottom summary numbers
  fmt_num_small <- function(x) {
    ifelse(
      is.na(x),
      "",
      paste0("$", formatC(x, format = "f", digits = digits), "$")
    )
  }
  
  fmt_colname <- function(x) {
    x <- latex_escape(x)
    if (grepl("/", x, fixed = TRUE)) {
      parts <- strsplit(x, "/", fixed = TRUE)[[1]]
      paste0("\\shortstack{", parts[1], "\\\\", parts[2], "}")
    } else {
      x
    }
  }
  
  # Clean strings
  tbl_chr <- apply(tbl_raw, c(1, 2), clean_cell)
  
  # Identify columns: assets + FROM
  has_from <- from_col_name %in% cn
  asset_cols <- if (has_from) setdiff(cn, from_col_name) else cn
  if (!has_from) {
    tbl_chr <- cbind(tbl_chr, FROM = "")
    cn <- c(cn, from_col_name)
  }
  
  # Find special rows
  to_idx  <- match(to_row_name, rn)
  inc_idx <- match(inc_row_name, rn)
  net_idx <- match(net_row_name, rn)  # optional, may be NA
  
  if (is.na(to_idx) || is.na(inc_idx)) {
    stop(sprintf(
      "Cannot find required rows '%s' and/or '%s' in rownames(dca$TABLE).",
      to_row_name, inc_row_name
    ))
  }
  
  main_end <- min(to_idx, inc_idx) - 1
  if (main_end < 1) stop("Main block has no rows (check TO/Inc.Own positions).")
  
  # Numeric version
  suppressWarnings({
    num_mat <- apply(tbl_chr, c(1, 2), as.numeric)
  })
  rownames(num_mat) <- rn
  colnames(num_mat) <- cn
  
  # Main block
  main_mat <- num_mat[seq_len(main_end), c(asset_cols, from_col_name), drop = FALSE]
  main_fmt <- apply(main_mat, c(1, 2), fmt_num)
  
  # Summary rows
  to_num  <- num_mat[to_idx,  asset_cols]
  inc_num <- num_mat[inc_idx, asset_cols]
  
  to_vals  <- fmt_num_small(to_num)
  inc_vals <- fmt_num_small(inc_num)
  
  # Net row: use existing row if present, otherwise compute TO - FROM
  if (!is.na(net_idx)) {
    net_num <- num_mat[net_idx, asset_cols]
  } else {
    # FROM vector by asset row name (preferred)
    from_by_asset <- rep(NA_real_, length(asset_cols))
    names(from_by_asset) <- asset_cols
    
    common_assets <- intersect(asset_cols, rn[seq_len(main_end)])
    if (length(common_assets) > 0) {
      from_by_asset[common_assets] <- num_mat[common_assets, from_col_name]
    } else {
      # fallback by position
      k <- min(length(asset_cols), main_end)
      from_by_asset[seq_len(k)] <- num_mat[seq_len(k), from_col_name]
    }
    
    net_num <- to_num - from_by_asset
  }
  net_vals <- fmt_num_small(net_num)
  
  # Total spillover index note
  N <- length(asset_cols)
  tsi_num <- sum(num_mat[seq_len(main_end), from_col_name], na.rm = TRUE)
  tsi_denom <- 100 * N
  tsi_pct <- if (is.finite(tsi_num) && is.finite(tsi_denom) && tsi_denom > 0) {
    100 * tsi_num / tsi_denom
  } else {
    NA_real_
  }
  
  note <- if (is.finite(tsi_pct)) {
    sprintf(
      "\\hfill\\shortstack[r]{Total spillover index\\\\ (%.1f/%d): %.1f\\%%}",
      tsi_num, tsi_denom, tsi_pct
    )
  } else {
    ""
  }
  
  header_cells <- c(
    "",
    vapply(asset_cols, fmt_colname, character(1)),
    "\\shortstack{Directional\\\\FROM others}"
  )
  header_line <- paste(header_cells, collapse = " & ")
  
  col_spec <- paste0("@{}l", paste(rep("c", N + 1), collapse = ""), "@{}")
  
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\setlength{\\tabcolsep}{2pt}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{adjustbox}{width=\\linewidth}",
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\toprule",
    paste0(header_line, " \\\\"),
    "\\midrule"
  )
  
  # Main rows
  for (i in seq_len(nrow(main_fmt))) {
    row_label <- latex_escape(rn[i])
    row_line <- paste(c(row_label, main_fmt[i, ]), collapse = " & ")
    lines <- c(lines, paste0(row_line, " \\\\"))
  }
  
  # Summary rows + Net spillover row
  lines <- c(
    lines,
    "\\midrule",
    paste0(
      "\\shortstack{Directional\\\\TO others} & ",
      paste(to_vals, collapse = " & "),
      " & \\\\[4pt]"
    ),
    paste0(
      "\\shortstack{Directional\\\\including\\\\own} & ",
      paste(inc_vals, collapse = " & "),
      " & \\\\[4pt]"
    ),
    paste0(
      "Net spillover & ",
      paste(net_vals, collapse = " & "),
      " & ",
      note,
      " \\\\"
    ),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{adjustbox}",
    "\\begin{flushleft}
    \\raggedright
    \\footnotesize
    \\textit{Notes:} This table reports lower-tail risk spillovers across currency markets. The sample period spans from 2020-10-18 to 2026-02-03. The $(i,j)$-th entry represents the estimated contribution of innovations in market $j$ to the forecast error variance of market $i$. The off-diagonal column sums (labeled \\emph{Directional TO others}) and row sums (labeled \\emph{Directional FROM others}) measure directional spillovers transmitted to and received from other markets, respectively, while the difference between the two (\\emph{TO minus FROM}) captures net risk spillovers. The total risk spillover index, reported in the lower-right corner of the table, is defined as the ratio of the sum of all off-diagonal elements to the total sum of all elements (including diagonals), expressed as a percentage. The spillover table thus provides an approximate input--output decomposition of the total spillover index. Risk spillovers are estimated using a generalized vector autoregressive (VAR) framework following Diebold and Yilmaz (2012), based on forecast error variance decompositions of the reciprocal lower-tail index.
    \\end{flushleft}",
    "\\end{table}"
  )
  
  cat(paste(lines, collapse = "\n"))
}

load("./result/DY/dca.RData")

sink("./table_txt/spillover_down_table.txt")
dca_table_to_latex(dca_down)
sink()

sink("./table_txt/spillover_up_table.txt")
dca_table_to_latex(dca_up)
sink()




