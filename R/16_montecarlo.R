source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(plm); library(parallel) })
source(here::here("R", "_montecarlo.R"))
source(here::here("R", "_mc_report.R"))

# targets
real_fit <- function(ms) {
  z <- complete_on(filter(prep_dynamic(FILE_P5, ms$dep, ms$inc), sample == 1),
                   c(ms$dep, ms$inc, "Ldep", "Linc"))
  list(ols = unname(ce(fit_ols(z, ms$dep, c("Ldep", "Linc"), FALSE), "Ldep")["est"]),
       fe  = unname(ce(fit_ols(z, ms$dep, c("Ldep", "Linc"), TRUE),  "Ldep")["est"]))
}

# calibrations
cals <- list()
for (ms in MEASURES) {
  tg <- real_fit(ms)
  a_fe <- alpha_for_fe(ms$dep, ms$inc, tg$fe)
  base <- calibrate_dgp(ms$dep, ms$inc, beta = 0)
  cals[[ms$label]] <- list(
    ms = ms, target = tg,
    alpha = c(gmm = base$alpha_hat, fe = a_fe),
    cal = list(gmm = base,
               fe = calibrate_dgp(ms$dep, ms$inc, alpha = a_fe, beta = 0)))
  cat(sprintf("%-14s real FE %.3f  alpha(GMM) %.3f  alpha(FE-matched) %.3f\n",
              ms$label, tg$fe, base$alpha_hat, a_fe))
}

# fidelity
for (lab in names(cals)) {
  cc <- cals[[lab]]$cal$gmm
  set.seed(MC_SEED)
  sim <- filter(sim_panel(cc, beta = 0), sample == 1)
  got <- complete_on(sim, c("y", "dLdep", "dLinc", "L2inc"))
  stopifnot(nrow(got) == cc$gmm_obs, n_distinct(got$code) == cc$gmm_ctry,
            n_distinct(got$period) == cc$gmm_periods)
  cat(sprintf("%-14s simulated estimation sample matches the real one: %d obs, %d countries\n",
              lab, cc$gmm_obs, cc$gmm_ctry))
}

# workers
cl <- makeCluster(MC_WORKERS)
invisible(clusterEvalQ(cl, {
  source(here::here("R", "00_setup.R"))
  suppressPackageStartupMessages(library(plm))
  source(here::here("R", "_montecarlo.R"))
  NULL
}))

one_draw <- function(r, cal, alpha, design, seed) {
  set.seed(seed + r)
  d <- sim_panel(cal, alpha = alpha, beta = 0, design = design)
  out <- run_ladder(d, cal$dep, cal$inc)
  out$rep <- r
  out
}
environment(one_draw) <- globalenv()

# run
draws <- list()
for (lab in names(cals)) {
  cc <- cals[[lab]]
  for (kind in names(MC_CALIB)) for (dsg in names(MC_DESIGNS)) {
    t0 <- Sys.time()
    seed <- MC_SEED + 1000L * match(lab, names(cals)) +
      100L * match(kind, names(MC_CALIB)) + 10L * match(dsg, names(MC_DESIGNS))
    res <- parLapply(cl, seq_len(MC_REPS), one_draw,
                     cal = cc$cal[[kind]], alpha = cc$alpha[[kind]],
                     design = dsg, seed = seed)
    tb <- bind_rows(res)
    tb$measure <- lab; tb$design <- unname(MC_DESIGNS[dsg])
    tb$calibration <- unname(MC_CALIB[kind])
    tb$alpha_true <- cc$alpha[[kind]]; tb$beta_true <- 0
    draws[[length(draws) + 1]] <- tb
    cat(sprintf("%-14s %-40s %-32s %5.1f s\n", lab, MC_CALIB[kind], MC_DESIGNS[dsg],
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
draws <- bind_rows(draws) |>
  mutate(estimator = factor(estimator,
    levels = c(vapply(LADDER_BASE, function(x) x$label, character(1)),
               vapply(GMM_SPECS, function(x) x$label, character(1)))))
write_csv(draws, file.path(PATH_OUTPUT, "montecarlo_draws.csv.gz"))

sm <- mc_summary(draws)
write_csv(sm, file.path(PATH_OUTPUT, "montecarlo.csv"))
stopifnot(nrow(sm) == length(cals) * length(MC_CALIB) * length(MC_DESIGNS) *
            length(LADDER_ORDER) * 2L,
          all(sm$done > 0.9))
cat("Every estimator returned an estimate in at least 90 per cent of the draws.\n")
big_worker <- function(r, cal, a) {
  set.seed(MC_SEED + 77000L + r)
  run_ladder(sim_panel(cal, alpha = a, beta = 0, design = "stationary"),
             cal[["dep"]], cal[["inc"]])
}
environment(big_worker) <- globalenv()


# consistency
big_check <- function(lab, kind, times = 6L, reps = 24L) {
  cc <- cals[[lab]]
  cal <- inflate_panel(cc$cal[[kind]], times)
  out <- bind_rows(parLapply(cl, seq_len(reps), big_worker,
                             cal = cal, a = cc[["alpha"]][[kind]]))
  out |> group_by(key, estimator) |>
    summarise(bias_big = mean(beta_hat, na.rm = TRUE),
              bias_big_alpha = mean(alpha_hat - cc$alpha[[kind]], na.rm = TRUE),
              .groups = "drop") |>
    mutate(measure = lab, times = times)
}
big <- big_check(names(cals)[1], "gmm")
write_csv(big, file.path(PATH_OUTPUT, "montecarlo_largeN.csv"))

small <- filter(sm, measure == names(cals)[1], grepl("difference GMM", calibration),
                grepl("Mean", design), param == "beta") |>
  select(key, bias_small = bias)
shrink <- left_join(big, small, by = "key") |>
  mutate(shrunk = abs(bias_big) < abs(bias_small) + 1e-8)
stopifnot(mean(shrink$shrunk, na.rm = TRUE) >= 0.7)
stopCluster(cl)

cat(sprintf("With %d times the countries the income bias shrinks for %d of %d estimators.\n",
            big$times[1], sum(shrink$shrunk, na.rm = TRUE), nrow(shrink)))


# report
write_mc_report(sm, big, shrink, cals)
saveRDS(list(sm = sm, big = big, shrink = shrink, cals = cals),
        file.path(PATH_OUTPUT, "montecarlo_state.rds"))
cat("Monte Carlo writeup written.\n")
