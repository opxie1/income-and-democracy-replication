# designs
DESIGNS <- list(
  list(name = "Paper's design (lags of democracy only)", short = "paper",
       levels = function(dep, inc) dep, extra = "L2inc"),
  list(name = "Symmetric (lags of democracy and income)", short = "symmetric",
       levels = function(dep, inc) c(dep, inc), extra = character()))

# frame
est_frame <- function(dep, inc) {
  d <- prep_dynamic(FILE_P5, dep, inc)
  s <- filter(d, sample == 1)
  list(full = d, est = complete_on(s, c("y", "dLdep", "dLinc", "L2inc")),
       fe = as.numeric(ce(fit_ols(s, dep, c("Ldep", "Linc"), TRUE), "Linc")["est"]))
}

# rules
agg_sweep <- function(dep, inc, panel_label) {
  fr <- est_frame(dep, inc)
  rows <- list()
  for (dg in DESIGNS) for (sc in AGG_SCHEMES) for (L in c(LAG_WINDOW, Inf)) {
    fitc <- function(model) fit_abgmm(fr$full, fr$est, dep_level = dg$levels(dep, inc),
                                      endog = c("dLdep", "dLinc"), inst_extra = dg$extra,
                                      lag_max = L, scheme = sc, model = model)
    m <- tryCatch(fitc("onestep"), error = function(e) NULL)
    if (is.null(m)) next
    m2 <- tryCatch(fitc("twostep"), error = function(e) NULL)
    rows[[length(rows) + 1]] <- tibble(
      panel = panel_label, design = dg$name, design_short = dg$short,
      scheme = sc, scheme_label = unname(SCHEME_LABEL[sc]),
      lag_max = L, income = unname(m$coef["dLinc"]), income_se = unname(m$se["dLinc"]),
      dem = unname(m$coef["dLdep"]), n_gmm = m$n_gmm, n_inst = m$n_inst,
      n_par = m$n_par, obs = m$nobs, countries = m$n_country,
      hansen_p = if (is.null(m2)) NA_real_ else m2$hansen_p,
      hansen_df = if (is.null(m2)) NA_integer_ else m2$hansen_df)
  }
  tb <- bind_rows(rows)
  stopifnot(nrow(tb) == length(DESIGNS) * length(AGG_SCHEMES) * (length(LAG_WINDOW) + 1L),
            length(unique(tb$obs)) == 1L, length(unique(tb$countries)) == 1L)
  list(tab = tb, bench = tibble(panel = panel_label, fe = fr$fe))
}

# families
family_sweep <- function(dep, inc, panel_label) {
  fr <- est_frame(dep, inc)
  nper <- n_distinct(fr$full$period)
  fit <- function(...) fit_abgmm(fr$full, fr$est, dep_level = dep,
                                 endog = c("dLdep", "dLinc"), inst_extra = "L2inc", ...)
  grab <- function(m) c(unname(m$coef["dLinc"]), unname(m$se["dLinc"]), m$n_gmm)

  stopifnot(
    max(abs(grab(fit(scheme = "block", block_size = 1L)) - grab(fit(scheme = "none")))) < 1e-9,
    max(abs(grab(fit(scheme = "block", block_size = nper)) - grab(fit(scheme = "lag")))) < 1e-9,
    max(abs(grab(fit(scheme = "period_geom", rho = 1)) - grab(fit(scheme = "period")))) < 1e-9)

  rows <- list()
  add <- function(family, setting, label, g)
    rows[[length(rows) + 1]] <<- tibble(panel = panel_label, family = family,
      setting = setting, label = label,
      income = g[1], income_se = g[2], n_gmm = as.integer(g[3]))
  for (b in seq_len(nper))
    add("Blocking the years", b, sprintf("%d", b), grab(fit(scheme = "block", block_size = b)))
  for (r in GEOM_RHO)
    add("Fading out older lags", r, sprintf("%.1f", r), grab(fit(scheme = "period_geom", rho = r)))
  bind_rows(rows)
}

# subsets
subset_sweep <- function(dep, inc, panel_label) {
  fr <- est_frame(dep, inc)
  fit <- function(keep, model) fit_abgmm(fr$full, fr$est, dep_level = dep,
    endog = c("dLdep", "dLinc"), inst_extra = "L2inc", keep_gmm = keep, model = model)
  total <- fit(NULL, "onestep")$n_gmm
  sizes <- unique(pmin(seq(4L, total, by = 8L), total))
  set.seed(SUBSET_SEED)
  rows <- list()
  for (s in sizes) for (r in seq_len(SUBSET_DRAWS)) {
    keep <- sort(sample.int(total, s))
    m1 <- fit(keep, "onestep")
    m2 <- fit(keep, "twostep")
    rows[[length(rows) + 1]] <- tibble(panel = panel_label, size = s, draw = r,
      income = unname(m1$coef["dLinc"]), income_se = unname(m1$se["dLinc"]),
      n_inst = m1$n_inst, hansen_p = m2$hansen_p, hansen_df = m2$hansen_df)
  }
  bind_rows(rows)
}

# outputs
write_aggregation_outputs <- function(tab, bench, fam, sub) {
  txt <- c()
  for (pl in unique(tab$panel)) {
    txt <- c(txt, sprintf("== %s (fixed effects %.3f) ==", pl, filter(bench, panel == pl)$fe),
             "Each cell: estimate (SE) [lagged-level instruments, overid p]",
             "'exact id' means there is nothing left over to test")
    for (dg in unique(tab$design)) {
      txt <- c(txt, sprintf("-- %s --", dg),
               sprintf("%-38s %-30s %-30s", "Collapse of the lagged levels",
                       sprintf("Lags %d-%d", min(LAG_WINDOW), 4L), "All lags"))
      for (sc in AGG_SCHEMES) {
        cell <- function(L) {
          r <- filter(tab, panel == pl, design == dg, scheme == sc,
                      if (is.infinite(L)) is.infinite(lag_max) else lag_max == L)
          if (nrow(r) == 0) "" else sprintf("%s (%s) [%d, %s]",
            num(r$income), num(r$income_se), r$n_gmm,
            if (r$hansen_df == 0L) "exact id" else num(r$hansen_p, 2))
        }
        txt <- c(txt, sprintf("%-38s %-30s %-30s", SCHEME_LABEL[sc], cell(4), cell(Inf)))
      }
      txt <- c(txt, "")
    }
  }
  writeLines(txt, file.path(PATH_OUTPUT, "aggregation.txt"))

  plotdat <- filter(tab, is.finite(lag_max))
  fig <- ggplot(plotdat, aes(lag_max, income, color = scheme_label)) +
    geom_hline(data = bench, aes(yintercept = fe), color = "grey45", linetype = "dotted") +
    geom_line() + geom_point(size = 1.4) +
    facet_grid(design ~ panel, scales = "free_y") +
    scale_x_continuous(breaks = LAG_WINDOW) +
    labs(x = "Longest lag used to build the instruments",
         y = "Estimated effect of income", color = "Instruments built as",
         title = "How the lagged levels collapse into instruments changes the answer",
         subtitle = "Difference GMM, one-step; dotted line is the fixed-effects estimate") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")
  ggsave(file.path(PATH_OUTPUT, "aggregation.png"), fig, width = 11, height = 6, dpi = 150)

  famfig <- ggplot(fam, aes(setting, income)) +
    geom_hline(data = bench, aes(yintercept = fe), color = "grey45", linetype = "dotted") +
    geom_line(color = "steelblue") + geom_point(color = "steelblue", size = 1.6) +
    facet_grid(panel ~ family, scales = "free") +
    labs(x = "Block size (left panel) or weight on each older lag (right panel)",
         y = "Estimated effect of income",
         title = "Turning the collapse up and down by degrees",
         subtitle = paste("Left: years grouped into blocks of this size; 1 is uncollapsed and",
                          "the largest is Roodman's rule.\nRight: each older lag multiplied",
                          "by this factor before adding; 1 is the plain sum.")) +
    theme_minimal(base_size = 11)
  ggsave(file.path(PATH_OUTPUT, "aggregation_families.png"), famfig,
         width = 10, height = 5.5, dpi = 150)

  fixed_ref <- filter(tab, is.infinite(lag_max), design_short == "paper") |>
    select(panel, scheme, scheme_label, income, n_inst)
  subfig <- ggplot(sub, aes(n_inst, income)) +
    geom_point(alpha = 0.3, size = 1, color = "grey45") +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                color = "grey20", linewidth = 0.6) +
    geom_hline(data = bench, aes(yintercept = fe), color = "grey45", linetype = "dotted") +
    geom_point(data = fixed_ref, aes(n_inst, income, color = scheme_label),
               size = 3, shape = 18) +
    facet_wrap(~panel, scales = "free_y") +
    labs(x = "Instruments in the draw (lagged levels plus year dummies)",
         y = "Estimated effect of income", color = "Deliberate rule",
         title = "Roodman's footnote 7: random subsets of the instrument set",
         subtitle = paste(SUBSET_DRAWS, "draws at each size, one-step difference GMM.",
                          "Diamonds are the deliberate rules at their own",
                          "\ninstrument counts; dotted line is the fixed-effects estimate.")) +
    theme_minimal(base_size = 11)
  ggsave(file.path(PATH_OUTPUT, "aggregation_subsets.png"), subfig,
         width = 11, height = 5, dpi = 150)

  subsum <- sub |> group_by(panel, size) |>
    summarise(n_inst = first(n_inst), med = median(income), lo = min(income), hi = max(income),
              sd = sd(income), med_se = median(income_se),
              rej = mean(hansen_p < CI_LEVEL, na.rm = TRUE), .groups = "drop")
  write_csv(subsum, file.path(PATH_OUTPUT, "aggregation_subsets_summary.csv"))
  subtxt <- c("Random subsets of the uncollapsed instrument set (Roodman footnote 7).",
              sprintf("%d draws at each size; one-step difference GMM for the coefficient,",
                      SUBSET_DRAWS),
              "two-step for the Hansen test. Paper's instrument design, all lags.",
              sprintf("%-16s %6s %8s %8s %8s %8s %8s", "Panel", "Insts", "median",
                      "min", "max", "sd", "rej 5%"))
  for (i in seq_len(nrow(subsum))) {
    r <- subsum[i, ]
    subtxt <- c(subtxt, sprintf("%-16s %6d %8s %8s %8s %8s %8s", r$panel, r$n_inst,
      num(r$med), num(r$lo), num(r$hi), num(r$sd), num(r$rej, 2)))
  }
  writeLines(subtxt, file.path(PATH_OUTPUT, "aggregation_subsets.txt"))
  list(fixed_ref = fixed_ref, subsum = subsum)
}
