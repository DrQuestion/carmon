dss_parameters <- function(x, n) {
  ## estimate means
  # Add smoothed?
  estprob <- dss_compute.mean.noSmooth(x, n)

  ## estimate dispersion
  phi <- dss_est.dispersion.BSseq(x, n, estprob)

  ## weight the counts
  wt <- 1 / (1 + (n - 1) * phi)
  wt <- wt / mean(wt)
  x.wt <- x * wt
  n.wt <- n * wt

  ## re-estimate means
  estprob <- dss_compute.mean.noSmooth(x.wt, n.wt)

  list(prob = estprob, phi = phi)
}

dss_compute.mean.noSmooth <- function(X, N) {
  p <- X / N

  const <- mean(p, na.rm = TRUE)
  p <- (rowSums(X) + const) / (rowSums(N) + 1)

  nreps <- ncol(N)
  res <- matrix(rep(p, nreps), ncol = nreps)
  return(res)
}

dss_est.dispersion.BSseq <- function(X, N, estprob) {
  prior <- dss_est.prior.BSseq.logN(X, N)
  dss_dispersion.shrinkage.BSseq(X, N, prior, estprob)
}

dss_est.prior.BSseq.logN <- function(X, N) {
  ## rowMeans = DelayedArray::rowMeans
  ## rowSums = DelayedArray::rowSums

  if (ncol(X) == 1) { ## single rep
    return(c(-3, 1))
  }

  ## keep sites with large coverage and no missing data
  ix <- rowMeans(N > 10) == 1 & rowSums(N == 0) == 0
  if (sum(ix) < 50) {
    warning("The coverages are too low. Cannot get good estimations of prior. Use arbitrary prior N(-3,1).")
    return(c(-3, 1))
  }

  X <- X[ix, , drop = FALSE]
  N <- N[ix, , drop = FALSE]
  ## compute sample mean/var
  p <- X / N
  mm <- rowMeans(p)
  mm[mm == 0] <- 1e-5
  mm[mm == 1] <- 1 - 1e-5
  vv <- stats::rowVars(p)
  phi <- vv / mm / (1 - mm)
  ## exclude those with vv==0. Those are sites with unobservable phis.
  ## But this will over estimate the prior.
  ## What will be the consequences????
  phi <- phi[vv > 0]
  lphi <- log(phi[phi > 0])
  prior.mean <- stats::median(lphi, na.rm = TRUE)
  prior.sd <- stats::IQR(lphi, na.rm = TRUE) / 1.39

  ## It seems this over-estimates the truth. Need to use the tricks in
  ## my biostat paper to remove the over-estimation. To be done later.
  c(prior.mean, prior.sd)
}

dss_dispersion.shrinkage.BSseq <- function(X, N, prior, estprob) {
  ## penalized likelihood function
  plik.logN <- function(size, X, prob, m0, tau, phi) {
    -(sum(dbetabinom(size, X, prob, exp(phi))) +
      stats::dnorm(phi, mean = m0, sd = tau, log = TRUE))
  }

  ## for CG sites with no coverage, use prior
  shrk.phi <- exp(rep(prior[1], nrow(N)))

  ## deal with estprob, make it a matrix if not.
  if (!is.matrix(estprob)) {
    estprob <- as.matrix(estprob)
  }

  ## skip those without coverage
  ix <- rowSums(N > 0) > 0
  X2 <- X[ix, , drop = FALSE]
  N2 <- N[ix, , drop = FALSE]
  estprob2 <- estprob[ix, , drop = FALSE]
  shrk.phi2 <- rep(0, nrow(X2))

  for (i in 1:nrow(X2)) {
    shrk.one <- stats::optimize(
      f = plik.logN, size = N2[i, ], X = X2[i, ], prob = estprob2[i, ],
      m0 = prior[1], tau = prior[2],
      interval = c(-5, log(0.99)), tol = 1e-3
    )
    shrk.phi2[i] <- exp(shrk.one$minimum)
  }

  shrk.phi[ix] <- shrk.phi2

  return(shrk.phi)
}

pbetabinom <- function(q, n, prob, phi, log.p = FALSE) {
  pbetabinom.ab(q,
    size = n,
    shape1 = prob * (1 - phi) / phi,
    shape2 = (1 - prob) * (1 - phi) / phi,
    # limit.prob = prob,  # 20180217, as phi = 0.
    log.p = log.p
  )
}

pbetabinom.ab <- function(q, size, shape1, shape2, log.p = FALSE) {
  LLL <- max(
    length(q), length(size), length(shape1),
    length(shape2)
  )

  ind3 <- !is.na(size) & !is.na(q) & q > size
  q[ind3] <- size[ind3] # Useful if q == Inf as it makes q finite

  ans <- q # Retains names(q)
  ans[] <- NA_real_ # Handles NAs in size, shape1, etc. hopefully

  for (ii in 1:LLL) {
    qstar <- floor(q[ii])
    qvec <- if (!is.na(qstar) && qstar >= 0) {
      0:qstar
    } else {
      NA
    }
    ans[ii] <-
      sum(dbetabinom.ab(
        q = qvec,
        size = size[ii],
        shape1 = shape1[ii],
        shape2 = shape2[ii]
      ))
  }

  ind4 <- !is.na(size) & !is.na(q) &
    !is.na(shape1) & !is.na(shape2) & q < 0
  ans[ind4] <- 0

  if (log.p) log(ans) else ans
}

dbetabinom <- function(size, x, prob, phi, log.p = TRUE) {
  ## first convert prob/phi to alpha/beta
  tmp <- 1 / phi - 1
  alpha <- prob * tmp
  beta <- tmp - alpha
  v <- lchoose(size, x) - lbeta(beta, alpha) + lbeta(size - x + beta, x + alpha)
  if (!log.p) {
    return(exp(v))
  } else {
    return(v)
  }
}

dbetabinom.ab <- function(q, size, shape1, shape2) {
  choose(size, q) * beta(q + shape1, size - q + shape2) / lbeta(shape1, shape2)
}
