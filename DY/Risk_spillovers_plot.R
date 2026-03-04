#' @title Dynamic from total directional connectedness plot
#' @description Visualize dynamic from total directional connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @return Return connectedness plot
#' @export
PlotFROM = function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), ...) {
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$config$approach)) {
    ca = list(ca)
  }
  x = dca$FROM
  date = as.Date(rownames(x))
  t = length(date)
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  if (is.null(ylim[1])) {
    lower = min(x)
  }
  if (is.null(ylim[2])) {
    upper = max(x)
  }
  k_row = ceiling(sqrt(k))
  k_col = ceiling(k/k_row)
  lower = ylim[1]
  upper = ylim[2]
  
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/FROM.pdf"), width=10, height=7)
  par(mfcol=c(k_row,k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  if (length(dim(dca$NET))>2) {
    for (i in 1:k) {
      x_ = x[,i,]
      if (is.null(lower)) {
        lower = min(x)
      }
      if (is.null(upper)) {
        upper = max(apply(x,1:2,sum))
      }
      plot(date, x_[,1], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper))#, ...)
      grid(NA, NULL, lty=2)
      for (j in 1:ncol(x_)) {
        polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,j])),col=j, border=j)
      }
      for (j in 1:ncol(x_)) {
        lines(date, x_[,j],col=j)
      }
      abline(h=0, lty=3)
      legend("topleft", colnames(x_), fill=c(1:(ncol(x_)+1)), bty="n")
      box()
    }
  } else {
    if (is.null(lower)) {
      lower = min(x)
    }
    if (is.null(upper)) {
      upper = max(x)
    }
    for (i in 1:k) {
      plot(date, x[,i], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
      grid(NA, NULL, lty=2)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[,i])),col=1, border=1)
      if (!is.null(ca)) {
        for (il in 1:length(ca)) {
          lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$FROM[,i], col=il+1)
        }
      }
      abline(h=0, lty=3)
      box()
    }
  }
  if (!is.null(path)) dev.off()
}

#' @title Dynamic influence connectedness plot
#' @description Visualize dynamic influence connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param selection Indidcator of the illustrated series
#' @param ... Arguments to be passed to methods, such as graphidcal parameters (see par).
#' @return Return connectedness plot
#' @export
PlotINF = function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), selection=NULL, ...) {
  message("The influence connectedness index is implemented according to:\n Gabauer, D. (2021). Dynamic measures of asymmetric & pairwise connectedness within an optimal currency area: Evidence from the ERM I system. Journal of Multinational Financial Management, 60, 100680.")
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$config$approach)) {
    ca = list(ca)
  }
  x = dca$INFLUENCE
  if (is.null(x)) {
    stop(paste(dca$config$approach, "has no INFLUENCE values."))
  }
  date = as.Date(dimnames(x)[[3]])
  t = length(date)
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  lower = ylim[1]
  upper = ylim[2]
  
  kk = k*(k-1)/2
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/INFLUENCE.pdf"), width=10, height=7)
  if (is.null(selection)) {
    k_row = ceiling(sqrt(kk))
    k_col = ceiling(kk/k_row)
    par(mfcol=c(k_row, k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  } else {
    k_row = ceiling(sqrt(k))
    k_col = ceiling(k/k_row)
    par(mfcol=c(k_row, k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  }
  if (length(dim(dca$NET))>2) {
    for (j in 1:k) {
      for (i in 1:k) {
        if (i>j) {
          if (i==selection || j==selection || is.null(selection)) {
            x_ = x[i,j,,]
            if (is.null(lower)) {
              lower = min(x)
            }
            if (is.null(upper)) {
              upper = max(x)
            }
            plot(date, x_[,1], type="l", main=paste(NAMES[j],"-",NAMES[i]), las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper))#, ...)
            grid(NA, NULL, lty=2)
            polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,1])),col=2, border=2)
            for (l in ncol(x_):1) {
              polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,l])),col=l+1, border=l+1)
            }
            for (l in 1:ncol(x_)) {
              lines(date, x_[,l],col=l+1)
            }
            abline(h=0, lty=3)
            legend("topleft", colnames(x_), fill=2:(1+ncol(x_)), bty="n")
            box()
          }
        }
      }
    }
  } else {
    if (is.null(lower)) {
      lower = min(x)
    }
    if (is.null(upper)) {
      upper = max(x)
    }
    for (j in 1:k) {
      for (i in 1:k) {
        if (i>j) {
          if (i==selection || j==selection || is.null(selection)) {
            plot(date, x[i,j,], type="l", main=paste(NAMES[j],"-",NAMES[i]), las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
            grid(NA, NULL, lty=2)
            polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[i,j,])),col=1, border=1)
            if (!is.null(ca)) {
              for (il in 1:length(ca)) {
                lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$INFLUENCE[i,j,], col=il+1)
              }
            }
            abline(h=0, lty=3)
            box()
          }
        }
      }
    }
  }
  if (!is.null(path)) dev.off()
}


# #' @title Dynamic net total directional connectedness plot
# #' @description Visualize dynamic net total directional connectedness
# #' @param dca Connectedness object
# #' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
# #' @param path Path where plots should be saved
# #' @param ylim A vector including the lower and upper limit of the y-axis
# #' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
# #' @return Return connectedness plot
# #' @export
# PlotNET = function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), ...) {
#   if (!is.null(path)) {
#     if (!dir.exists(path)) {
#       dir.create(path)
#     }
#   }
#   if (length(ca)>0 && !is.null(ca$config$approach)) {
#     ca = list(ca)
#   }
#   x = dca$NET
#   date = as.Date(rownames(x))
#   t = length(date)
#   k = ncol(x)
#   NAMES = colnames(x)
#   if (is.null(NAMES)) {
#     NAMES = 1:k
#   }
#   lower = ylim[1]
#   upper = ylim[2]
# 
#   k_row = ceiling(sqrt(k))
#   k_col = ceiling(k/k_row)
#   oldpar = par(no.readonly=TRUE)
#   on.exit(par(oldpar))
#   if (!is.null(path)) pdf(file=paste0(path, "/NET.pdf"), width=10, height=7)
#   par(mfcol=c(k_row,k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
#   if (length(dim(dca$NET))>2) {
#     for (i in 1:k) {
#       x_ = x[,i,]
#       if (is.null(lower)) {
#         lower = min(c(min(x_), min(apply(x_,1,sum))))
#       }
#       if (is.null(upper)) {
#         upper = max(c(max(x_), max(apply(x_,1,sum))))
#       }
#       plot(date, x_[,1], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper))#, ...)
#       grid(NA, NULL, lty=2)
#       for (j in 1:ncol(x_)) {
#         polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,j])),col=j, border=j)
#       }
#       for (j in 1:ncol(x_)) {
#         lines(date, x_[,j],col=j)
#       }
#       abline(h=0, lty=3)
#       legend("topleft", colnames(x_), fill=c(1:(ncol(x_)+1)), bty="n")
#       box()
#     }
#   } else {
#     if (is.null(lower)) {
#       lower = min(x)
#     }
#     if (is.null(upper)) {
#       upper = max(x)
#     }
#     for (i in 1:k) {
#       plot(date, x[,i], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper))#, ...)
#       grid(NA, NULL, lty=2)
#       polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[,i])),col=1, border=1)
#       if (!is.null(ca)) {
#         for (il in 1:length(ca)) {
#           lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$NET[,i], col=il+1)
#         }
#       }
#       abline(h=0, lty=3)
#       box()
#     }
#   }
#   if (!is.null(path)) dev.off()
# }



#' @title Dynamic net total directional connectedness plot
#' @description Visualize dynamic net total directional connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param separate If TRUE, save one PDF per series (currency pair)
#' @param file_prefix Prefix for output pdf filenames when separate=TRUE
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @return Return connectedness plot (invisibly)
#' @export
PlotNET <- function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL),
                    separate=TRUE, file_prefix="NET_", ...) {
  
  if (!is.null(path)) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  }
  
  if (length(ca) > 0 && !is.null(ca$config$approach)) ca <- list(ca)
  
  x <- dca$NET
  date <- as.Date(rownames(x))
  t <- length(date)
  k <- ncol(x)
  
  NAMES <- colnames(x)
  if (is.null(NAMES)) NAMES <- as.character(1:k)
  
  base_lower <- ylim[1]
  base_upper <- ylim[2]
  
  safe_name <- function(s) gsub("[^A-Za-z0-9_\\-]+", "_", s)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  if (!isTRUE(separate)) {
    # fallback: original multi-panel behavior (optional)
    if (!is.null(path)) pdf(file=paste0(path, "/NET.pdf"), width=10, height=7)
    on.exit(if (!is.null(path)) dev.off(), add = TRUE)
    
    k_row <- ceiling(sqrt(k))
    k_col <- ceiling(k / k_row)
    par(mfcol=c(k_row,k_col), oma=c(0,0,0,0)+0.5, mar=c(1,1,1,1)+.5, mgp=c(1,0.4,0))
    
    for (i in 1:k) {
      lower_i <- base_lower; upper_i <- base_upper
      if (is.null(lower_i)) lower_i <- min(x[,i], na.rm=TRUE)
      if (is.null(upper_i)) upper_i <- max(x[,i], na.rm=TRUE)
      
      plot(date, x[,i], type="l", main=NAMES[i], las=1, xlab="", ylab="",
           xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower_i, upper_i), ...)
      grid(NA, NULL, lty=2)
      polygon(c(date, rev(date)), c(rep(0,t), rev(x[,i])), col=1, border=1)
      abline(h=0, lty=3)
      box()
    }
    return(invisible(TRUE))
  }
  
  # --- separate = TRUE: one PDF per series ---
  for (i in 1:k) {
    
    # compute per-series y-limits (recommended)
    lower_i <- base_lower
    upper_i <- base_upper
    
    if (length(dim(dca$NET)) > 2) {
      x_ <- x[, i, , drop=FALSE]
      x_ <- x_[, 1, ]  # t x m matrix
      
      if (is.null(lower_i)) lower_i <- min(c(min(x_), min(apply(x_,1,sum))), na.rm=TRUE)
      if (is.null(upper_i)) upper_i <- max(c(max(x_), max(apply(x_,1,sum))), na.rm=TRUE)
      
      fn <- if (!is.null(path)) paste0(path, "/", file_prefix, safe_name(NAMES[i]), ".pdf") else NULL
      
      if (!is.null(fn)) pdf(file=fn, width=10, height=7)
      
      tryCatch({
        par(mfrow=c(1,1), oma=c(0,0,0,0)+0.5, mar=c(3,3,2,1)+.5, mgp=c(1,0.4,0))
        par(family="sans")  # helps on Windows if a custom font is missing
        
        plot(date, x_[,1], type="l", main=NAMES[i], las=1, xlab="", ylab="",
             xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower_i, upper_i), ...)
        grid(NA, NULL, lty=2)
        
        for (j in 1:ncol(x_)) {
          polygon(c(date, rev(date)), c(rep(0,t), rev(x_[,j])), col=j, border=j)
        }
        for (j in 1:ncol(x_)) lines(date, x_[,j], col=j)
        
        abline(h=0, lty=3)
        legend("topleft", colnames(x_), fill=1:ncol(x_), bty="n")
        box()
      }, finally = {
        if (!is.null(fn)) dev.off()
      })
      
    } else {
      if (is.null(lower_i)) lower_i <- min(x[,i], na.rm=TRUE)
      if (is.null(upper_i)) upper_i <- max(x[,i], na.rm=TRUE)
      
      fn <- if (!is.null(path)) paste0(path, "/", file_prefix, safe_name(NAMES[i]), ".pdf") else NULL
      if (!is.null(fn)) pdf(file=fn, width=10, height=7)
      
      tryCatch({
        par(mfrow=c(1,1), oma=c(0,0,0,0)+0.5, mar=c(3,3,2,1)+.5, mgp=c(1,0.4,0))
        par(family="sans")
        
        plot(date, x[,i], type="l", main=NAMES[i], las=1, xlab="", ylab="",
             xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower_i, upper_i), ...)
        grid(NA, NULL, lty=2)
        polygon(c(date, rev(date)), c(rep(0,t), rev(x[,i])), col=1, border=1)
        
        if (!is.null(ca)) {
          for (il in 1:length(ca)) {
            lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$NET[,i], col=il+1)
          }
        }
        
        abline(h=0, lty=3)
        box()
      }, finally = {
        if (!is.null(fn)) dev.off()
      })
    }
  }
  
  invisible(TRUE)
}


#' @title Dynamic net pairwise connectedness plot
#' @description Visualize dynamic net pairwise connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param selection Indicator of the illustrated series
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @return Return connectedness plot
#' @export
PlotNPDC = function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), selection=NULL, ...) {
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$config$approach)) {
    ca = list(ca)
  }
  x = dca$NPDC
  if (is.null(x)) {
    stop(paste(ca$config$approach, "has no NPDC."))
  }
  date = as.Date(dimnames(x)[[3]])
  t = length(date)
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  lower = ylim[1]
  upper = ylim[2]
  
  kk = k*(k-1)/2
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/NPDC.pdf"), width=10, height=7)
  if (is.null(selection)) {
    k_row = ceiling(sqrt(kk))
    k_col = ceiling(kk/k_row)
    par(mfcol=c(k_row, k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  } else {
    k_row = ceiling(sqrt(k))
    k_col = ceiling(k/k_row)
    par(mfcol=c(k_row, k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  }
  if (length(dim(dca$NET))>2) {
    for (j in 1:k) {
      for (i in 1:k) {
        if (i>j) {
          if (i==selection || j==selection || is.null(selection)) {
            x_ = x[i,j,,]
            if (is.null(lower)) {
              lower = min(c(min(x), min(apply(x_,1,sum))))
            }
            if (is.null(upper)) {
              upper = max(c(max(x), max(apply(x_,1,sum))))
            }
            plot(date, x_[,1], type="l", main=paste(NAMES[j],"-",NAMES[i]), las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
            grid(NA, NULL, lty=2)
            for (l in 1:ncol(x_)) {
              polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,l])),col=l, border=l)
            }
            for (l in 1:ncol(x_)) {
              lines(date, x_[,l],col=l)
            }
            abline(h=0, lty=3)
            legend("topleft", colnames(x_), fill=1:ncol(x_), bty="n")
            box()
          }
        }
      }
    }
  } else {
    if (is.null(lower)) {
      lower = min(x)
    }
    if (is.null(upper)) {
      upper = max(x)
    }
    for (j in 1:k) {
      for (i in 1:k) {
        if (i==selection || j==selection || is.null(selection)) {
          if (i>j) {
            plot(date, x[i,j,], type="l", main=paste(NAMES[j],"-",NAMES[i]), las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
            grid(NA, NULL, lty=2)
            polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[i,j,])),col=1, border=1)
            if (!is.null(ca)) {
              for (il in 1:length(ca)) {
                lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$NPDC[i,j,], col=il+1)
              }
            }
            abline(h=0, lty=3)
            box()
          }
        }
      }
    }
  }
  if (!is.null(path)) dev.off()
}

#' @title Dynamic net pairwise transmission plot
#' @description Visualize dynamic net total directional connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ... Arguments to be passed to methods, such as graphidcal parameters (see par).
#' @return Return connectedness plot
#' @export
PlotNPT = function(dca, ca=NULL, path=NULL, ...) {
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$config$approach)) {
    ca = list(ca)
  }
  x = dca$NPT
  if (is.null(x)) {
    stop(paste(dca$config$approach, "has no NPT."))
  }
  date = as.Date(rownames(x))
  t = length(date)
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  k_row = ceiling(sqrt(k))
  k_col = ceiling(k/k_row)

  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/NPT.pdf"), width=10, height=7)
  par(mfcol=c(k_row,k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  if (length(dim(dca$NET))>2) {
    for (i in 1:k) {
      plot(date, x[,i], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(0,k-1), ...)
      grid(NA, NULL, lty=2)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[,i])),col=1, border=1)
      if (!is.null(ca)) {
        for (il in 1:length(ca)) {
          lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$NPT[,i], col=il+1)
        }
      }
      abline(h=0, lty=3)
      box()
    }
  } else {
    for (i in 1:k) {
      x_ = x[,i,]
      plot(date, apply(x_,1,sum), type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(0,k-1), ...)
      grid(NA, NULL, lty=2)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(apply(x_,1,sum))),col=1, border=1)
      for (j in ncol(x_):1) {
        polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,j])),col=j+1, border=j+1)
      }
      for (j in 1:ncol(x_)) {
        lines(date, x_[,j],col=j+1)
      }
      lines(date, apply(x_,1,sum), col=1)
      abline(h=0, lty=3)
      legend("topleft", c("Total", colnames(x_)), fill=1:(ncol(x_)+1), bty="n")
      box()
    }
  }
  if (!is.null(path)) dev.off()
}


#' @title Network plot
#' @description Visualize net pairwise or pairwise connectedness measures
#' @param dca Connectedness object
#' @param path Path where plots should be saved
#' @param method Either visualizing NPDC or PCI
#' @param name_length Length of variable names in the network plot
#' @param threshold Threshold for bivariate connections between 0 and 1 (edges below are not drawn)
#' @param threshold_color Threshold for coloring edges between 0 and 1 (edges above are colored red)
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @return Return connectedness plot
#' @export
#' @import igraph
PlotNetwork = function(dca, method="NPDC", path=NULL, name_length=NULL,
                       threshold=0.25, threshold_color=0.50, ...) {

  # ---- select connectedness measure ----
  if (method=="NPDC") {
    x = dca$NPDC
  } else if (method=="PCI") {
    x = dca$PCI
  } else {
    stop("This method does not exist")
  }
  
  date = as.Date(dimnames(x)[[3]])
  t = length(date)
  k = ncol(x)
  
  # ---- variable names ----
  NAMES = dimnames(x)[[1]]
  if (is.null(NAMES)) {
    NAMES = 1:k
  } else {
    NAMES = colnames(x)
    if (!is.null(name_length)) {
      NAMES = substr(NAMES, 1, name_length)
    }
  }
  
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar))
  
  # ---- dimension handling ----
  if (length(dim(x)) > 3) {
    kk = dim(x)[4]
    k1 = ceiling(sqrt(kk))
    k2 = ceiling(kk/k1)
  } else {
    kk = k1 = k2 = 1
    x = array(x, c(k,k,t,1))
  }
  
  par(mfrow = c(k1,k2), oma = c(0,0,0,0), mar = c(0,0,0,0), mgp = c(0, 0, 0))
  
  # ---- open PDF device ----
  if (!is.null(path)) {
    pdf(file = path, width = 10, height = 10)
  }
  
  for (ijk in 1:kk) {
    
    x_ = t(apply(x[,,,ijk], 1:2, mean))
    x_ = ifelse(x_ < 0, 0, x_)
    
    colnames(x_) = rownames(x_) = NAMES
    diag(x_) = 0
    
    # ---- normalize safely ----
    x_ = x_ - min(x_, na.rm=TRUE)
    max_val <- max(x_, na.rm=TRUE)
    if (max_val > 0) {
      x_ = x_ / max_val
    }
    
    x_[x_ < threshold] = 0
    m = 5 * x_
    
    # ---- network construction ----
    if (isTRUE(all.equal(x_, t(x_)))) {
      
      gr = igraph::graph.adjacency(m, mode="undirected", weighted=TRUE)
      lo = igraph::layout_in_circle(gr)
      net = igraph::graph.adjacency(m, mode="undirected", weighted=TRUE, diag=FALSE)
      
      color = rep("#A6CEE3", k)
      nn = rep(1, k)
      
    } else {
      
      gr = igraph::graph.adjacency(m, mode="undirected", weighted=TRUE)
      net = igraph::graph.adjacency(m, mode="directed", weighted=TRUE, diag=FALSE)
      lo <- igraph::layout_in_circle(gr)
      
      color = rep("#A6CEE3", k)
      
      nn = apply(apply(x[,,,ijk],1:2,mean),1,sum)
      
      color[nn > 0] = "#FDBF6F"
      
      max_nn <- max(abs(nn))
      if (max_nn > 0) {
        nn = abs(nn / max_nn)
      } else {
        nn = rep(0, length(nn))
      }
    }
    
    # ---- edge styles ----
    edge_col = ifelse(igraph::E(net)$weight >= 5 * threshold_color, "red", "grey50")
    
    edge_arrow_size  = 0.6
    edge_arrow_width = 1.2
    edge_width = igraph::E(net)$weight * 0.7
    
    igraph::plot.igraph(
      net,
      vertex.label = igraph::V(net)$name,
      layout = lo,
      vertex.label.family = "sans",
      vertex.label.cex = 0.8,
      vertex.size = 5 + nn * 10,
      vertex.label.font = 2,
      vertex.color = color,
      vertex.frame.color = color,
      vertex.label.color = "black",
      mark.col = "steelblue4",
      edge.width = edge_width,
      edge.color = edge_col,
      edge.arrow.size = edge_arrow_size,
      edge.arrow.width = edge_arrow_width,
      edge.curved = 0.25,
      edge.lty = 1
    )
  }
  
  if (!is.null(path)) {
    dev.off()
  }
}



#' @title Dynamic pairwise connectedness plot
#' @description Visualize dynamic pairwise connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param selection Indidcator of the illustrated series
#' @param ... Arguments to be passed to methods, such as graphidcal parameters (see par).
#' @return Return connectedness plot
#' @export
PlotPCI = function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), selection=NULL, ...) {
  message("The pairwise connectedness index is implemented according to:\n Gabauer, D. (2021). Dynamic measures of asymmetric & pairwise connectedness within an optimal currency area: Evidence from the ERM I system. Journal of Multinational Financial Management, 60, 100680.")
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$config$approach)) {
    ca = list(ca)
  }
  x = dca$PCI
  if (is.null(x)) {
    stop(paste(dca$config$approach, "has no PCI."))
  }
  date = as.Date(dimnames(x)[[3]])
  t = length(date)
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  lower = ylim[1]
  upper = ylim[2]
  
  kk = k*(k-1)/2
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/PCI.pdf"), width=10, height=7)
  if (is.null(selection)) {
    k_row = ceiling(sqrt(kk))
    k_col = ceiling(kk/k_row)
    par(mfcol=c(k_row, k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  } else {
    k_row = ceiling(sqrt(k))
    k_col = ceiling(k/k_row)
    par(mfcol=c(k_row, k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  }
  if (length(dim(dca$NET))>2) {
    for (j in 1:k) {
      for (i in 1:k) {
        if (i>j) {
          if (i==selection || j==selection || is.null(selection)) {
            x_ = x[i,j,,]
            x[which(x>=99.99999,arr.ind=T)] = 0
            if (is.null(lower)) {
              lower = min(x)
            }
            if (is.null(upper)) {
              upper = max(x)
            }
            plot(date, x_[,1], type="l", main=paste(NAMES[j],"-",NAMES[i]), las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
            grid(NA, NULL, lty=2)
            polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,1])),col=2, border=2)
            for (l in ncol(x_):1) {
              polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,l])),col=l+1, border=l+1)
            }
            for (l in 1:ncol(x_)) {
              lines(date, x_[,l],col=l+1)
            }
            abline(h=0, lty=3)
            legend("topleft", colnames(x_), fill=2:(1+ncol(x_)), bty="n")
            box()
          }
        }
      }
    }
  } else {
    if (is.null(lower)) {
      lower = min(x)
    }
    if (is.null(upper)) {
      upper = max(x)
    }
    for (j in 1:k) {
      for (i in 1:k) {
        if (i>j) {
          if (i==selection || j==selection || is.null(selection)) {
            plot(date, x[i,j,], type="l", main=paste(NAMES[j],"-",NAMES[i]), las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
            grid(NA, NULL, lty=2)
            polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[i,j,])),col=1, border=1)
            if (!is.null(ca)) {
              for (il in 1:length(ca)) {
                lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$PCI[i,j,], col=il+1)
              }
            }
            abline(h=0, lty=3)
            box()
          }
        }
      }
    }
  }
  if (!is.null(path)) dev.off()
}


#' @title Dynamic total connectedness plot
#' @description Visualize dynamic total connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @return Return connectedness plot
#' @import graphics
#' @import grDevices
#' @export
PlotTCI <- function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), fill=FALSE, ...) {
  
  # normalize ca into a list if a single connectedness object is provided
  if (length(ca) > 0 && !is.null(ca$config$approach)) {
    ca <- list(ca)
  }
  
  x <- dca$TCI
  date <- as.Date(rownames(x))
  t <- length(date)
  
  # handle x shape robustly
  x_mat <- as.matrix(x)
  k <- ncol(x_mat)
  
  lower <- ylim[1]
  upper <- ylim[2]
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  # open PDF device
  if (!is.null(path)) {
    pdf(file = path, width = 10, height = 5)
    on.exit(dev.off(), add = TRUE)
  }
  
  par(mfrow=c(1,1), oma=c(0,0,0,0) + 0.5, mar=c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  
  # case 1: rolling / multi-layer (NET has >2 dims)
  if (length(dim(dca$NET)) > 2) {
    
    x_ <- x_mat
    
    if (is.null(lower)) lower <- min(x_, na.rm = TRUE)
    if (is.null(upper)) upper <- max(x_, na.rm = TRUE)
    
    plot(date, x_[,1], type="l", main="", las=1,
         xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02,
         ylim=c(lower, upper))
    
    grid(NA, NULL, lty=2)
    
    # optional fill (off by default)
    if (isTRUE(fill)) {
      for (j in 1:ncol(x_)) {
        polygon(c(date, rev(date)), c(rep(0, t), rev(x_[,j])),
                col=j, border=j)
      }
    }
    
    # draw lines
    for (j in 1:ncol(x_)) {
      lines(date, x_[,j], col=j)
    }
    
    legend("topleft", colnames(x_), col=1:ncol(x_), lty=1, bty="n")
    
    abline(h=0, lty=3)
    box()
    
  } else {
    # case 2: standard single TCI series
    
    if (is.null(lower)) lower <- 0
    if (is.null(upper)) upper <- 100
    
    y <- as.numeric(x_mat)
    
    plot(date, y, type="l", main="", las=1,
         xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02,
         ylim=c(lower, upper), ...)
    
    grid(NA, NULL, lty=2)
    
    # optional fill (off by default)
    if (isTRUE(fill)) {
      polygon(c(date, rev(date)), c(rep(0, t), rev(y)),
              col=1, border=1)
    }
    
    # optional comparisons
    if (!is.null(ca)) {
      for (il in 1:length(ca)) {
        
        ca_tci <- ca[[il]]$TCI
        ca_date <- as.Date(rownames(ca_tci))
        
        lines(ca_date, as.numeric(ca_tci), col=il+1)
        
        gTCI <- ca[[il]]$gTCI
        if (!is.null(gTCI)) {
          gTCI_mat <- as.matrix(gTCI)
          for (ij in 1:ncol(gTCI_mat)) {
            lines(ca_date, gTCI_mat[, ij], col=ij+2)
          }
        }
      }
      
      # legend behavior kept similar to original (only for single ca)
      if (length(ca) == 1) {
        approach <- ca[[1]]$config$approach
        
        gTCI <- ca[[1]]$gTCI
        gTCI_mat <- if (!is.null(gTCI)) as.matrix(gTCI) else NULL
        
        if (!is.null(approach) && (approach == "Internal" || approach == "External")) {
          if (!is.null(gTCI_mat)) {
            legend("topleft",
                   c("TCI", paste("TCI", approach), colnames(gTCI_mat)),
                   col = 1:(ncol(gTCI_mat) + 2),
                   lty = 1, bty="n")
          } else {
            legend("topleft",
                   c("TCI", paste("TCI", approach)),
                   col = 1:2, lty=1, bty="n")
          }
        } else if (!is.null(approach) && (approach == "Inclusive" || approach == "Exclusive")) {
          legend("topleft",
                 c("TCI", paste("TCI", approach)),
                 col = 1:2, lty=1, bty="n")
        }
      }
    }
    
    box()
  }
  
  invisible(TRUE)
}

#' @title Dynamic to total directional connectedness plot
#' @description Visualize dynamic to total directional connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @return Return connectedness plot
#' @export
PlotTO = function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), ...) {
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$config$approach)) {
    ca = list(ca)
  }
  x = dca$TO
  date = as.Date(rownames(x))
  t = length(date)
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  k_row = ceiling(sqrt(k))
  k_col = ceiling(k/k_row)
  lower = ylim[1]
  upper = ylim[2]
  
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/TO.pdf"), width=10, height=7)
  par(mfcol=c(k_row,k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  if (length(dim(dca$NET))>2) {
    for (i in 1:k) {
      x_ = x[,i,]
      if (is.null(lower)) {
        lower = min(x)
      }
      if (is.null(upper)) {
        upper = max(apply(x,1:2,sum))
      }
      plot(date, x_[,1], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper))#, ...)
      grid(NA, NULL, lty=2)
      for (j in 1:ncol(x_)) {
        polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,j])),col=j, border=j)
      }
      for (j in 1:ncol(x_)) {
        lines(date, x_[,j],col=j)
      }
      abline(h=0, lty=3)
      legend("topleft", colnames(x_), fill=c(1:(ncol(x_)+1)), bty="n")
      box()
    }
  } else {
    if (is.null(lower)) {
      lower = min(x)
    }
    if (is.null(upper)) {
      upper = max(x)
    }
    for (i in 1:k) {
      plot(date, x[,i], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
      grid(NA, NULL, lty=2)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[,i])),col=1, border=1)
      if (!is.null(ca)) {
        for (il in 1:length(ca)) {
          lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$TO[,i], col=il+1)
        }
      }
      abline(h=0, lty=3)
      box()
    }
  }
  if (!is.null(path)) dev.off()
}



TCIplot<-function(datainput,flag){

tci<-data.frame(date=as.Date(row.names(datainput$TCI)),TCI=datainput$TCI) 



 ggplot(tci, aes(x = as.Date(date))) +
  geom_line(aes(y = TCI),size=1) +
 
  labs(x = "Date", y = "Values") +
  theme(
    
    
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.background = element_rect(fill = "white", color = NA),
  
    axis.title = element_text(size = 14, face = "bold")
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "12 months", expand = c(0, 0.05))
 



















}


PlotNetwork_updata = function(dca, method="NPDC", path=NULL, name_length=NULL, threshold=0.2, ...) {
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (method=="NPDC") {
    x = dca
  } else if (method=="PCI") {
    x = dca$PCI
  } else {
    stop("This method does not exists")
  }
 # date = as.Date(dimnames(x)[[3]])
  t = 1
  k = ncol(x)
  
  NAMES = dimnames(x)[[1]]
  if (is.null(NAMES)) {
    NAMES = 1:k
  } else {
    NAMES = colnames(x)
    if (!is.null(name_length)) {
      NAMES = substr(NAMES, 1, name_length)
    }
  }
  
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (length(dim(x))>3) {
    kk = dim(x)[4]
    k1 = ceiling(sqrt(kk))
    k2 = ceiling(kk/k1)
  } else {
    kk = k1 = k2 = 1
    x = array(x, c(k,k,t,1))
  }
  
  par(mfrow = c(k1,k2), oma = c(0,0,0,0), mar = c(0,0,0,0), mgp = c(0, 0, 0))
  if (!is.null(path)) pdf(file=paste0(path, "/NetworkPlot.pdf"), width=10, height=10)
    for (ijk in 1:kk) {
      x_ = t(apply(x[,,,ijk], 1:2, mean))
      x_ = ifelse(x_<0, 0, x_)
      colnames(x_) = rownames(x_) = NAMES
      diag(x_) = 0
      x_ = x_ - min(x_)
      x_ = x_ / max(x_)
      x_[x_<threshold] = 0
      m = 5 * x_
      if (isTRUE(all.equal(x_, t(x_)))) {
        gr = graph.adjacency(m, mode="undirected", weighted=TRUE)
        lo = layout_in_circle(gr)
        net = graph.adjacency(m, mode="undirected", weighted=TRUE, diag=FALSE)
        color = rep("steelblue4", k)
        nn = 1
      } else {
        gr = graph.adjacency(m, mode="undirected", weighted=TRUE)	
		
        net = graph.adjacency(m, mode="directed", weighted=TRUE, diag=FALSE)	
		lo <- layout_on_sphere(gr)
		#lo <- layout_on_grid(gr)
		
        color = rep("#4D4D4D", k)
        nn = apply(apply(x[,,,ijk],1:2,mean),1,sum)
		
		
        color[nn>0] = "#BFBFBF"
		
        nn = abs(nn/max(abs(nn))) 
      }
	 
			lo[43,2]<- -0.66740211 

	  
      plot.igraph(net,vertex.label=V(net)$name,layout=lo, vertex.label.family = "SimHei",vertex.label.cex=0.8, vertex.size=5+nn*10, vertex.label.font=2,vertex.color=color, vertex.frame.color=color, vertex.label.color="black", mark.col="steelblue4",
                  edge.width=E(net)$weight*0.7	, edge.color="grey50", edge.arrow.size=0.5, edge.curved=0.1, lty=6)
    }
  if (!is.null(path)) dev.off()
}

