# partial
partial_out <- function(dat, cols, exog, fe = FE_TERMS) {
  W <- model.matrix(reformulate(c(fe, exog)), data = dat)
  q <- qr(W)
  list(R = qr.resid(q, as.matrix(dat[, cols, drop = FALSE])), rank = q$rank)
}

cluster_scale <- function(G, n, rank, k) (G / (G - 1)) * ((n - 1) / (n - rank - k))

# moments
moment_parts <- function(dat, dep, endog, inst, exog, cluster = "code") {
  po <- partial_out(dat, c(dep, endog, inst), exog)
  yt <- po$R[, dep]; xt <- po$R[, endog]; Zt <- po$R[, inst, drop = FALSE]
  g <- dat[[cluster]]
  ZZi <- solve(crossprod(Zt))
  uy <- as.numeric(yt - Zt %*% (ZZi %*% crossprod(Zt, yt)))
  ux <- as.numeric(xt - Zt %*% (ZZi %*% crossprod(Zt, xt)))
  list(yt = yt, xt = xt, Zt = Zt, g = g, rank = po$rank, k = ncol(Zt),
       n = nrow(Zt), G = n_distinct(g),
       P = rowsum(Zt * yt, g), Q = rowsum(Zt * xt, g),
       Ps = rowsum(Zt * uy, g), Qs = rowsum(Zt * ux, g),
       nk = as.numeric(table(g)))
}

# anderson-rubin
ar_pcurve <- function(dat, dep, endog, inst, exog, grid, cluster = "code") {
  mp <- moment_parts(dat, dep, endog, inst, exog, cluster)
  cst <- cluster_scale(mp$G, mp$n, mp$rank, mp$k)
  stat_at <- function(bs) vapply(bs, function(b0) {
    s <- colSums(mp$P - b0 * mp$Q)
    as.numeric(crossprod(s, solve(crossprod(mp$Ps - b0 * mp$Qs), s))) / cst
  }, numeric(1))
  p_at <- function(bs) pf(stat_at(bs) / mp$k, mp$k, mp$G - 1, lower.tail = FALSE)
  chisq <- stat_at(grid)
  c(mp, list(cst = cst, stat_at = stat_at, p_at = p_at, chisq = chisq,
             p = pf(chisq / mp$k, mp$k, mp$G - 1, lower.tail = FALSE),
             p_chisq = pchisq(chisq, mp$k, lower.tail = FALSE)))
}

# likelihood
lr_stat <- function(QS, QT, cross) 0.5 * (QS - QT + sqrt(pmax((QS + QT)^2 - 4 * cross, 0)))

clr_null_draws <- function(k, nsim, seed) {
  if (k == 1L) return(NULL)
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  w <- rchisq(nsim, k - 1)
  out <- list(QS = rnorm(nsim)^2 + w, w = w)
  if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv)
  out
}

clr_pvalue <- function(LR, QT, k, draws) {
  if (k == 1L) return(pchisq(LR, 1, lower.tail = FALSE))
  mean(lr_stat(draws$QS, QT, QT * draws$w) >= LR)
}

# moreira
clr_pcurve <- function(dat, dep, endog, inst, exog, grid,
                       nsim = CLR_SIM, seed = CLR_SEED) {
  po <- partial_out(dat, c(dep, endog, inst), exog)
  Y <- cbind(po$R[, dep], po$R[, endog]); Zt <- po$R[, inst, drop = FALSE]
  k <- length(inst); dfree <- nrow(Y) - po$rank - k
  ZZ <- crossprod(Zt); ZY <- crossprod(Zt, Y)
  Om <- crossprod(Y - Zt %*% solve(ZZ, ZY)) / dfree
  ei <- eigen(ZZ, symmetric = TRUE)
  ZZmh <- ei$vectors %*% diag(1 / sqrt(ei$values), k) %*% t(ei$vectors)
  Omi <- solve(Om)
  draws <- clr_null_draws(k, nsim, seed)
  p_at <- function(bs) vapply(bs, function(b0) {
    b <- c(1, -b0); a <- c(b0, 1)
    S <- ZZmh %*% (ZY %*% b) / sqrt(as.numeric(t(b) %*% Om %*% b))
    Tt <- ZZmh %*% (ZY %*% (Omi %*% a)) / sqrt(as.numeric(t(a) %*% Omi %*% a))
    QS <- sum(S^2); QT <- sum(Tt^2); QST <- sum(S * Tt)
    clr_pvalue(lr_stat(QS, QT, QS * QT - QST^2), QT, k, draws)
  }, numeric(1))
  list(p = p_at(grid), p_at = p_at, k = k)
}

# kleibergen
clr_robust_pcurve <- function(dat, dep, endog, inst, exog, grid, cluster = "code",
                              nsim = CLR_SIM, seed = CLR_SEED) {
  mp <- moment_parts(dat, dep, endog, inst, exog, cluster)
  k <- mp$k
  cst <- cluster_scale(mp$G, mp$n, mp$rank, k)
  qtot <- colSums(mp$Q)
  Qc <- mp$Qs - outer(mp$nk / mp$n, colSums(mp$Qs))
  Vqq <- crossprod(Qc) * cst
  draws <- clr_null_draws(k, nsim, seed)
  parts_at <- function(b0) {
    mg <- mp$Ps - b0 * mp$Qs
    m <- colSums(mp$P - b0 * mp$Q)
    Vmi <- solve(crossprod(mg) * cst)
    Cqm <- crossprod(Qc, mg) * cst
    qh <- qtot - as.numeric(Cqm %*% (Vmi %*% m))
    AR <- as.numeric(crossprod(m, Vmi %*% m))
    rk <- as.numeric(crossprod(qh, solve(Vqq - Cqm %*% Vmi %*% t(Cqm), qh)))
    K <- as.numeric(crossprod(m, Vmi %*% qh))^2 / as.numeric(crossprod(qh, Vmi %*% qh))
    c(AR = AR, K = K, J = max(AR - K, 0), rk = rk)
  }
  p_at <- function(bs) vapply(bs, function(b0) {
    z <- parts_at(b0)
    clr_pvalue(lr_stat(z["AR"], z["rk"], z["rk"] * z["J"]), z["rk"], k, draws)
  }, numeric(1))
  list(p = p_at(grid), p_at = p_at, parts_at = parts_at, k = k, G = mp$G, cst = cst)
}

# montiel-pflueger
effective_f <- function(dat, endog, inst, exog, fs) {
  po <- partial_out(dat, c(endog, inst), exog)
  A <- crossprod(po$R[, inst, drop = FALSE])
  pi_hat <- fs$coefficients[inst]
  V <- fs$vcov[inst, inst, drop = FALSE]
  as.numeric(crossprod(pi_hat, A %*% pi_hat) / sum(diag(V %*% A)))
}

# cragg-donald
cragg_donald <- function(dat, endog, inst, exog, fe = FE_TERMS) {
  po <- partial_out(dat, c(endog, inst), exog, fe = fe)
  X <- po$R[, endog, drop = FALSE]; Zt <- po$R[, inst, drop = FALSE]
  k <- ncol(Zt)
  PX <- Zt %*% solve(crossprod(Zt), crossprod(Zt, X))
  Sig <- crossprod(X - PX) / (nrow(X) - po$rank - k)
  e <- eigen(Sig, symmetric = TRUE)
  Sih <- e$vectors %*% diag(1 / sqrt(e$values), ncol(X)) %*% t(e$vectors)
  min(eigen(Sih %*% crossprod(PX) %*% Sih / k, symmetric = TRUE, only.values = TRUE)$values)
}

# endpoints
refine_edge <- function(p_at, lo, hi, level) {
  f <- function(x) p_at(x) - level
  if (f(lo) * f(hi) > 0) return(NA_real_)
  uniroot(f, c(lo, hi), tol = 1e-8)$root
}

# sets
set_info <- function(grid, p, p_at = NULL, level = CI_LEVEL) {
  keep <- p >= level
  if (!any(keep)) return(list(text = "empty", lo = NA_real_, hi = NA_real_,
                              pieces = 0L, grid_bounded = TRUE))
  r <- rle(keep)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  idx <- which(r$values)
  step <- if (length(grid) > 1L) grid[2] - grid[1] else 0
  edge <- function(i, side) {
    j <- if (side == "lo") starts[idx[i]] else ends[idx[i]]
    if ((side == "lo" && j == 1L) || (side == "hi" && j == length(grid)) || is.null(p_at))
      return(grid[j])
    out <- if (side == "lo") refine_edge(p_at, grid[j] - step, grid[j], level)
           else refine_edge(p_at, grid[j], grid[j] + step, level)
    if (is.na(out)) grid[j] else out
  }
  los <- vapply(seq_along(idx), function(i) edge(i, "lo"), numeric(1))
  his <- vapply(seq_along(idx), function(i) edge(i, "hi"), numeric(1))
  first <- starts[idx[1]]; last <- ends[idx[length(idx)]]
  list(text = paste(vapply(seq_along(idx), function(i) sprintf("[%s, %s]",
         if (starts[idx[i]] == 1L) "-inf" else sprintf("%.2f", los[i]),
         if (ends[idx[i]] == length(grid)) "+inf" else sprintf("%.2f", his[i])),
       character(1)), collapse = " U "),
       lo = los[1], hi = his[length(his)], pieces = length(idx),
       grid_bounded = first > 1L && last < length(grid))
}
