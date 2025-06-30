reconstruct <- function(data, p, method = "coglasso", sel_method = "xestars", cor_cutoff = 0.6) {
  if (method == "coglasso") {
    if (length(p) > 2) {
      # Until multiple datasets issue is solved for D > 2
      rec <- coglasso::bs(data, p, method = sel_method, c = 0)
    } else {
      rec <- coglasso::bs(data, p, method = sel_method)
    }
  } else if (method == "glasso") {
    rec <- coglasso::bs(data, lock_lambdas = TRUE)
  } else if (method == "correlation") {
    rec <- stats::cor(data)
    above_cutoff <- which(abs(rec) >= cor_cutoff, arr.ind = TRUE)
    rec <- rec[above_cutoff]
  }
  return(rec)
}
