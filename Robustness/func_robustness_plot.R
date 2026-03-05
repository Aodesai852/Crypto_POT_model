vecs_to_dt <- function(vectors) {
  stopifnot(is.list(vectors), length(vectors) > 0)
  ok <- vapply(vectors, function(x) !is.null(names(x)), logical(1))
  if (!all(ok)) stop("Every vector must be a *named* vector (names(x) cannot be NULL).")
  cols <- names(vectors)
  all_names <- unique(unlist(lapply(vectors, names), use.names = FALSE))
  DT <- data.table(name = all_names)
  for (col in cols) {
    x <- vectors[[col]]
    DT[, (col) := x[name]]
  }
  setkey(DT, name)
  DT
}

matrices_to_dt_by_time <- function(lst) {
  stopifnot(is.list(lst), length(lst) > 0)
  
  if (is.null(names(lst)) || any(names(lst) == "")) {
    stop("The list must have non-empty names to be used as column name prefixes.")
  }
  
  dts <- Map(function(mat, prefix) {
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    
    rn <- rownames(mat)
    if (is.null(rn)) stop("Each matrix must have rownames representing dates.")
    
    time_vec <- as.Date(rn)
    if (anyNA(time_vec)) stop("Date conversion failed. rownames(mat) must be coercible to Date.")
    
    if (!is.numeric(mat)) storage.mode(mat) <- "numeric"
    
    dt <- data.table(time = time_vec)
    dt <- cbind(dt, as.data.table(mat))
    
    # Ensure column names exist
    mcols <- names(dt)[-1]
    if (is.null(mcols) || any(mcols == "")) {
      setnames(dt, old = mcols, new = paste0("V", seq_along(mcols)))
      mcols <- names(dt)[-1]
    }
    
    # Prefix columns with list item name to avoid collisions
    setnames(dt, old = mcols, new = paste0(prefix, "_", mcols))
    
    dt
  }, lst, names(lst))
  
  out <- Reduce(function(a, b) merge(a, b, by = "time", all = TRUE), dts)
  setorder(out, time)
  out
}

rank45_plot <- function(dt) {
  stopifnot(data.table::is.data.table(dt))
  if (ncol(dt) < 3) stop("Input must have at least three columns: name, benchmark, alternative.")
  
  name_col <- names(dt)[1]
  y_col    <- names(dt)[2]
  x_col    <- names(dt)[3]
  
  d <- data.table::copy(dt)
  d[[name_col]] <- as.factor(d[[name_col]])
  
  rng <- range(c(d[[x_col]], d[[y_col]]), na.rm = TRUE)
  
  ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[name_col]]
    )
  ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "solid") +
    ggplot2::geom_point(size = 3, shape = 16, alpha = 0.9) +
    ggplot2::coord_equal(xlim = rng, ylim = rng, expand = TRUE) +
    ggplot2::labs(x = x_col, y = y_col, color = "Asset") +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 3, byrow = TRUE)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 9),
      legend.key.height = grid::unit(0.6, "lines"),
      legend.spacing.y = grid::unit(0.15, "lines"),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

plot_multi_series <- function(dt) {
  stopifnot(data.table::is.data.table(dt))
  
  if (ncol(dt) < 2) {
    stop("Input data.table must have at least two columns: time and one numeric series.")
  }
  
  time_col <- names(dt)[1]
  value_cols <- names(dt)[-1]
  
  d <- data.table::copy(dt)
  
  d_long <- data.table::melt(
    d,
    id.vars = time_col,
    measure.vars = value_cols,
    variable.name = "series",
    value.name = "value"
  )
  
  p <- ggplot2::ggplot(
    d_long,
    ggplot2::aes(
      x = .data[[time_col]],
      y = value,
      color = series
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.3,
      linetype = "solid",
      lineend  = "butt"
    ) +
    ggplot2::labs(
      x = "Time",
      y = "Value",
      color = "Series"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}
