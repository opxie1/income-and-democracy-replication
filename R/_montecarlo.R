# stationary
lyapunov <- function(A, S, tol = 1e-13, iter = 5000L) {
  V <- S
  for (i in seq_len(iter)) {
    Vn <- A %*% V %*% t(A) + S
    if (max(abs(Vn - V)) < tol) return(Vn)
    V <- Vn
  }
  V
}

match_fe_scale <- function(A, S, c_fe, target_var) {
  V <- lyapunov(A, S)
  M <- solve(diag(2) - A)
  want <- target_var - diag(V)
  if (any(want <= 0)) return(NULL)
  obj <- function(p) {
    s <- exp(p)
    Om <- matrix(c(s[1]^2, c_fe * s[1] * s[2], c_fe * s[1] * s[2], s[2]^2), 2, 2)
    sum((diag(M %*% Om %*% t(M)) - want)^2)
  }
  fit <- optim(log(sqrt(want) * abs(1 - diag(A))), obj, method = "BFGS")
  list(sd = exp(fit$par), V = V, conv = fit$convergence, obj = fit$value)
}

# calibrate
dgp_inputs <- local({
  cache <- list()
  function(dep, inc) {
    k <- paste(dep, inc, sep = "\r")
    if (is.null(cache[[k]])) {
      d <- prep_dynamic(FILE_P5, dep, inc)
      s <- filter(d, sample == 1)
      est <- complete_on(s, c("y", "dLdep", "dLinc", "L2inc"))
      m <- fit_abgmm(d, est, dep_level = dep, endog = c("dLdep", "dLinc"),
                     inst_extra = "L2inc")
      pd <- gmm_panel(s, dep, inc)
      rho <- unname(coef(pgmm(
        as.formula(sprintf("%s ~ lag(%s, 1) | lag(%s, %s)", inc, inc, inc, GMM_LAGS_ALL)),
        data = pd, effect = "twoways", model = "onestep", transformation = "d"))[1])
      z <- complete_on(s, c(dep, inc, "Ldep", "Linc"))
      sk <- d[!is.na(d[[dep]]) | !is.na(d[[inc]]), ]
      sk <- sk[, c("code", "year", "year_numeric", "period", "sample")]
      sk$miss_dep <- is.na(d[[dep]])[!is.na(d[[dep]]) | !is.na(d[[inc]])]
      sk$miss_inc <- is.na(d[[inc]])[!is.na(d[[dep]]) | !is.na(d[[inc]])]
      cache[[k]] <<- list(
        dep = dep, inc = inc, z = z, rho = rho,
        gamma = unname(ce(fit_ols(s, inc, c("Linc", "Ldep"), TRUE), "Ldep")["est"]),
        alpha_hat = unname(m$coef[["dLdep"]]), beta_hat = unname(m$coef[["dLinc"]]),
        ti = as.numeric(table(z$code)),
        target_var = c(mean(tapply(z[[dep]], z$year, var), na.rm = TRUE),
                       mean(tapply(z[[inc]], z$year, var), na.rm = TRUE)),
        skeleton = sk,
        n_country = n_distinct(sk$code), n_obs = nrow(sk),
        gmm_obs = nrow(est), gmm_ctry = n_distinct(est$code),
        gmm_periods = n_distinct(est$period),
        periods = sort(unique(sk$period)),
        fe_real = unname(ce(fit_ols(z, dep, c("Ldep", "Linc"), TRUE), "Ldep")["est"]),
        ols_real = unname(ce(fit_ols(z, dep, c("Ldep", "Linc"), FALSE), "Ldep")["est"]))
    }
    cache[[k]]
  }
})

calibrate_dgp <- function(dep, inc, alpha = NULL, beta = 0, strict = TRUE) {
  ip <- dgp_inputs(dep, inc)
  if (is.null(alpha)) alpha <- ip$alpha_hat
  z <- ip$z; rho <- ip$rho; gamma <- ip$gamma

  u <- z[[dep]] - alpha * z$Ldep - beta * z$Linc
  v <- z[[inc]] - rho * z$Linc - gamma * z$Ldep
  lam <- tapply(u, z$year, mean); phi <- tapply(v, z$year, mean)
  u <- u - lam[as.character(z$year)]; v <- v - phi[as.character(z$year)]

  eta <- tapply(u, z$code, mean); mu <- tapply(v, z$code, mean)
  eps <- u - eta[z$code]; nu <- v - mu[z$code]
  s_eps <- sd(eps); s_nu <- sd(nu); c_en <- cor(eps, nu)

  shrink <- function(v_hat, v_shock) max(v_hat - mean(v_shock / ip$ti), 1e-8)
  v_eta <- shrink(var(eta), s_eps^2)
  v_mu  <- shrink(var(mu),  s_nu^2)
  c_fe  <- (cov(eta, mu) - mean(c_en * s_eps * s_nu / ip$ti)) / sqrt(v_eta * v_mu)

  A <- matrix(c(alpha, beta, gamma, rho), 2, 2, byrow = TRUE)
  Sig <- matrix(c(s_eps^2, c_en * s_eps * s_nu, c_en * s_eps * s_nu, s_nu^2), 2, 2)
  fit <- match_fe_scale(A, Sig, c_fe, ip$target_var)
  if (is.null(fit) || fit[["obj"]] > 1e-8) {
    if (strict) stop(sprintf(
      "no feasible country-effect variance at alpha = %.3f, beta = %.3f for %s",
      alpha, beta, dep))
    return(NULL)
  }

  list(dep = dep, inc = inc, alpha = alpha, beta = beta, rho = rho, gamma = gamma,
       sd_eta = fit[["sd"]][1], sd_eps = s_eps, sd_mu = fit[["sd"]][2], sd_nu = s_nu,
       sd_eta_raw = sqrt(v_eta), sd_mu_raw = sqrt(v_mu),
       target_sd = sqrt(ip$target_var),
       alpha_hat = ip$alpha_hat, beta_hat = ip$beta_hat,
       fe_real = ip$fe_real, ols_real = ip$ols_real,
       cor_shock = c_en, cor_fe = c_fe, lambda = lam, phi = phi,
       skeleton = ip$skeleton, n_country = ip$n_country, n_obs = ip$n_obs,
       gmm_obs = ip$gmm_obs, gmm_ctry = ip$gmm_ctry, gmm_periods = ip$gmm_periods,
       periods = ip$periods)
}

# draw
sim_panel <- function(cal, alpha = cal$alpha, beta = cal$beta, design = "stationary",
                      burn = 40L, kappa = 1) {
  ids <- unique(cal$skeleton$code)
  N <- length(ids)
  pers <- cal$periods
  Tn <- length(pers)

  fe <- MASS::mvrnorm(N, c(0, 0), matrix(c(cal$sd_eta^2,
    cal$cor_fe * cal$sd_eta * cal$sd_mu, cal$cor_fe * cal$sd_eta * cal$sd_mu,
    cal$sd_mu^2), 2, 2))
  eta <- fe[, 1]; mu <- fe[, 2]

  Sig <- matrix(c(cal$sd_eps^2, cal$cor_shock * cal$sd_eps * cal$sd_nu,
                  cal$cor_shock * cal$sd_eps * cal$sd_nu, cal$sd_nu^2), 2, 2)

  A <- matrix(c(alpha, beta, cal$gamma, cal$rho), 2, 2, byrow = TRUE)
  lr <- solve(diag(2) - A, rbind(eta, mu))
  dcur <- lr[1, ]; ycur <- lr[2, ]

  if (design == "nonstationary") {
    dcur <- dcur + kappa * eta / (1 - alpha)
    ycur <- ycur + kappa * mu / (1 - cal$rho)
    burn <- 0L
  }

  step <- function(dprev, yprev, lam, phi) {
    e <- MASS::mvrnorm(N, c(0, 0), Sig)
    list(d = alpha * dprev + beta * yprev + eta + lam + e[, 1],
         y = cal$rho * yprev + cal$gamma * dprev + mu + phi + e[, 2])
  }
  for (b in seq_len(burn)) {
    nx <- step(dcur, ycur, 0, 0); dcur <- nx$d; ycur <- nx$y
  }

  dmat <- matrix(NA_real_, N, Tn); ymat <- matrix(NA_real_, N, Tn)
  for (tt in seq_len(Tn)) {
    yr <- as.character(cal$skeleton$year[match(pers[tt], cal$skeleton$period)])
    eff <- function(v) { x <- unname(v[yr]); if (length(x) != 1L || is.na(x)) 0 else x }
    nx <- step(dcur, ycur, eff(cal$lambda), eff(cal$phi))
    dcur <- nx$d; ycur <- nx$y
    dmat[, tt] <- dcur; ymat[, tt] <- ycur
  }

  key <- match(cal$skeleton$code, ids)
  pos <- match(cal$skeleton$period, pers)
  out <- cal$skeleton
  out[[cal$dep]] <- dmat[cbind(key, pos)]
  out[[cal$inc]] <- ymat[cbind(key, pos)]
  out[[cal$dep]][out$miss_dep] <- NA_real_
  out[[cal$inc]][out$miss_inc] <- NA_real_
  derive_dynamic(out, cal$dep, cal$inc)
}

# ladder
LADDER_BASE <- list(
  list(key = "ols", label = "Pooled OLS"),
  list(key = "fe",  label = "Fixed effects"),
  list(key = "ah",  label = "Anderson-Hsiao IV"),
  list(key = "abr", label = "Arellano-Bond, difference GMM (replication)"))

run_ladder <- function(d, dep, inc) {
  s <- filter(d, sample == 1)
  rows <- list()
  put <- function(key, label, a, ase, b, bse) {
    rows[[length(rows) + 1]] <<- tibble(key = key, estimator = label,
      alpha_hat = a, alpha_se = ase, beta_hat = b, beta_se = bse)
  }
  grab <- function(key, label, expr, ta, tb) {
    m <- tryCatch(expr, error = function(e) NULL)
    if (is.null(m)) return(put(key, label, NA, NA, NA, NA))
    put(key, label, ce(m, ta)["est"], ce(m, ta)["se"], ce(m, tb)["est"], ce(m, tb)["se"])
  }

  grab("ols", "Pooled OLS", fit_ols(s, dep, c("Ldep", "Linc"), FALSE), "Ldep", "Linc")
  grab("fe",  "Fixed effects", fit_ols(s, dep, c("Ldep", "Linc"), TRUE), "Ldep", "Linc")
  grab("ah",  "Anderson-Hsiao IV",
       fit_iv(s, "y", endog = c("dLdep", "dLinc"), inst = c("L2dep", "L2inc"),
              country_fe = FALSE), "dLdep", "dLinc")
  est <- complete_on(s, c("y", "dLdep", "dLinc", "L2inc"))
  grab("abr", "Arellano-Bond, difference GMM (replication)",
       fit_abgmm(d, est, dep_level = dep, endog = c("dLdep", "dLinc"),
                 inst_extra = "L2inc"), "dLdep", "dLinc")

  pd <- tryCatch(gmm_panel(s, dep, inc), error = function(e) NULL)
  for (sp in GMM_SPECS) {
    m <- if (is.null(pd)) NULL else
      tryCatch(fit_gmm_spec(pd, dep, inc, sp), error = function(e) NULL)
    if (is.null(m)) { put(sp$key, sp$label, NA, NA, NA, NA); next }
    co <- tryCatch(pgmm_co(m), error = function(e) NULL)
    if (is.null(co)) { put(sp$key, sp$label, NA, NA, NA, NA); next }
    ic <- co_row(co, inc); dc <- co_row(co, dep)
    put(sp$key, sp$label, dc[1], dc[2], ic[1], ic[2])
  }
  bind_rows(rows)
}

LADDER_ORDER <- c(vapply(LADDER_BASE, function(x) x$key, character(1)),
                  vapply(GMM_SPECS, function(x) x$key, character(1)))

# persistence
alpha_feasible <- function(dep, inc, a)
  !is.null(calibrate_dgp(dep, inc, alpha = a, beta = 0, strict = FALSE))

alpha_for_fe <- function(dep, inc, target, reps = 40L, seed = 4242L,
                         grid = seq(0.2, 0.9, by = 0.05)) {
  ok <- grid[vapply(grid, function(a) alpha_feasible(dep, inc, a), logical(1))]
  stopifnot(length(ok) >= 2L)
  mean_fe <- function(a) {
    cal <- calibrate_dgp(dep, inc, alpha = a, beta = 0)
    set.seed(seed)
    mean(vapply(seq_len(reps), function(r) {
      s <- filter(sim_panel(cal), sample == 1)
      unname(ce(fit_ols(s, dep, c("Ldep", "Linc"), TRUE), "Ldep")["est"])
    }, numeric(1)))
  }
  vals <- vapply(ok, mean_fe, numeric(1))
  i <- which(vals[-length(vals)] <= target & vals[-1] >= target)[1]
  if (is.na(i)) return(ok[which.min(abs(vals - target))])
  uniroot(function(a) mean_fe(a) - target, c(ok[i], ok[i + 1]), tol = 1e-3)[["root"]]
}

# metrics
mc_summary <- function(draws, level = CI_LEVEL) {
  crit <- qnorm(1 - level / 2)
  draws |>
    pivot_longer(c(alpha_hat, alpha_se, beta_hat, beta_se),
                 names_to = c("param", ".value"), names_sep = "_") |>
    mutate(truth = if_else(param == "alpha", alpha_true, beta_true)) |>
    group_by(measure, design, calibration, alpha_true, beta_true,
             key, estimator, param, truth) |>
    summarise(
      done     = mean(!is.na(hat)),
      mean_est = mean(hat, na.rm = TRUE),
      bias     = mean(hat - truth, na.rm = TRUE),
      med_bias = median(hat - truth, na.rm = TRUE),
      sd_est   = sd(hat, na.rm = TRUE),
      rmse     = sqrt(mean((hat - truth)^2, na.rm = TRUE)),
      mean_se  = mean(se, na.rm = TRUE),
      coverage = mean(abs(hat - truth) <= crit * se, na.rm = TRUE),
      reject0  = mean(abs(hat / se) > crit, na.rm = TRUE),
      .groups  = "drop")
}

# inflate
inflate_panel <- function(cal, times) {
  sk <- cal$skeleton
  cal$skeleton <- bind_rows(lapply(seq_len(times), function(j) {
    x <- sk; x$code <- paste0(x$code, "_", j); x
  }))
  cal$n_country <- cal$n_country * times
  cal$n_obs <- cal$n_obs * times
  cal
}
