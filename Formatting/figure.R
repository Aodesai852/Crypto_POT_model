rm(list=ls())
setwd("D:/MyFiles/EVT/BTC_POT/code/table")

generate_latex_images <- function(folder) {
  
  # Escape LaTeX special characters (for captions only)
  escape_latex_text <- function(x) {
    x <- gsub("\\\\", "\\\\textbackslash{}", x)
    x <- gsub("([{}$&#_^%])", "\\\\\\1", x, perl = TRUE)
    x
  }
  
  # List images
  files <- list.files(
    path = folder,
    pattern = "\\.(png|PNG|pdf|PDF)$",
    full.names = FALSE
  )
  files <- sort(files)
  
  for (i in seq_along(files)) {
    fname <- files[i]
    base  <- tools::file_path_sans_ext(fname)
    caption <- escape_latex_text(base)
    
    # Every row begins with \noindent (prevents indentation)
    if ((i - 1) %% 4 == 0) {
      cat("% ---- New row ----\n")
      cat("\\noindent\n")
    }
    
    # minipage for each image
    cat("\\begin{minipage}{0.24\\textwidth}\n")
    cat("  \\centering\n")
    cat(sprintf("  \\includegraphics[width=\\textwidth]{%s/%s}\n", folder, fname))
    cat(sprintf("  {\\small %s}\n", caption))
    cat("\\end{minipage}%\n")   # <-- % prevents line break!
    
    # Insert spacing except after the 4th image
    if (i %% 4 != 0) {
      cat("\\hfill%\n")
    } else {
      cat("\n\n")  # end of row
    }
  }
}



sink("../table_txt/tail_index_images.txt")
generate_latex_images("../result/down_tail/P_xi_exc")
generate_latex_images("../result/down_tail/P_std_sigma")
generate_latex_images("../result/up_tail/P_xi_exc")
generate_latex_images("../result/up_tail/P_std_sigma")
sink()

sink("../table_txt/up_NET_plot.txt")
generate_latex_images("../result/up_tail/NET_plots")
sink()

sink("../table_txt/down_NET_plot.txt")
generate_latex_images("../result/down_tail/NET_plots")
sink()

sink("../table_txt/price_plot.txt")
generate_latex_images("../result/price_plot")
sink()

sink("../table_txt/return_plot.txt")
generate_latex_images("../result/return_plot")
sink()

