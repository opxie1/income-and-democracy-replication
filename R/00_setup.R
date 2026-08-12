suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(readxl)
  library(estimatr)
  library(here)
})

here::i_am("R/00_setup.R")

# paths
PATH_KIT    <- here("replication-kit")
FILE_XLS    <- file.path(PATH_KIT, "Income-and-Democracy-Data-AER-adjustment.xls")

PATH_DATA   <- here("data")
PATH_OUTPUT <- here("output")
PATH_DOCS   <- here("docs")
for (p in c(PATH_DATA, PATH_OUTPUT, PATH_DOCS)) dir.create(p, showWarnings = FALSE)

FILE_RAW <- file.path(PATH_DATA, "panels_raw.rds")

FILE_P5  <- file.path(PATH_DATA, "panel_5yr.parquet")
FILE_PA  <- file.path(PATH_DATA, "panel_annual.parquet")
FILE_P10 <- file.path(PATH_DATA, "panel_10yr.parquet")
FILE_P20 <- file.path(PATH_DATA, "panel_20yr.parquet")
FILE_P25 <- file.path(PATH_DATA, "panel_25yr.parquet")
FILE_P50 <- file.path(PATH_DATA, "panel_50yr.parquet")

if (!file.exists(FILE_XLS)) {
  stop("Can't find the data workbook. Unzip the AER replication kit (openICPSR 113251) ",
       "into a folder called 'replication-kit' in the project root (see README).")
}

# columns
COLS_5YR <- c("code", "country", "code_numeric", "year", "year_numeric",
              "sample", "samplebalancefe", "samplebalancegmm", "socialist",
              "fhpolrigaug", "polity4", "lrgdpch", "nsave", "worldincome",
              "worlddemocracy", "laborshare", "lpop", "medage", "education",
              "age_veryyoung", "age_young", "age_midage", "age_old")
COLS_FREQ <- c("code", "country", "code_numeric", "year", "year_numeric",
               "sample", "fhpolrigaug", "polity4", "lrgdpch")
COLS_MAD  <- c("code", "country", "code_numeric", "year", "year_numeric",
               "sample", "noextrapolation", "madid", "polity4", "lrgdpmad")

# labels
VAR_LABELS <- c(
  code             = "Country code (3-letter)",
  country          = "Country name",
  code_numeric     = "Country numeric code",
  year             = "Year of observation",
  year_numeric     = "Panel time index (consecutive within panel)",
  sample           = "In base estimation sample (1 = yes)",
  samplebalancefe  = "In balanced fixed-effects sample, 1970-2000 (1 = yes)",
  samplebalancegmm = "In balanced GMM sample (1 = yes)",
  socialist        = "Former Soviet bloc country (1 = yes)",
  noextrapolation  = "Income not extrapolated in Maddison data (1 = yes)",
  madid            = "Maddison aggregation id (clustering unit for long-run panels)",
  fhpolrigaug      = "Augmented Freedom House political rights index (0-1)",
  polity4          = "Polity IV democracy index (rescaled 0-1)",
  lrgdpch          = "Log real GDP per capita (Penn World Tables)",
  lrgdpmad         = "Log real GDP per capita (Maddison)",
  nsave            = "Nominal savings rate, (Y-C-G)/Y",
  worldincome      = "Trade-weighted world income (instrument)",
  worlddemocracy   = "Trade-weighted world democracy",
  laborshare       = "Labor share of gross value added (%)",
  lpop             = "Log total population (in thousands)",
  medage           = "Median age of the population",
  education        = "Average years of schooling",
  age_veryyoung    = "Share of population age 0-15",
  age_young        = "Share of population age 15-30",
  age_midage       = "Share of population age 30-45",
  age_old          = "Share of population age 45-60"
)

apply_labels <- function(df) {
  for (col in names(df)) {
    if (!is.null(VAR_LABELS[[col]])) attr(df[[col]], "label") <- VAR_LABELS[[col]]
  }
  df
}

# lags
add_lags <- function(df, vars, ks, id = "code", time = "period") {
  stopifnot(!anyDuplicated(df[c(id, time)]))
  key <- paste(df[[id]], df[[time]], sep = "\r")
  for (v in vars) for (k in ks) {
    src_key <- paste(df[[id]], df[[time]] + k, sep = "\r")
    df[[paste0(v, "_l", k)]] <- df[[v]][match(key, src_key)]
  }
  df
}

# panels
read_panel <- function(file) {
  read_parquet(file) |> arrange(code, year_numeric) |> mutate(period = year_numeric)
}

complete_on <- function(df, vars) df[stats::complete.cases(df[, vars, drop = FALSE]), ]

ce <- function(m, term) {
  if (!is.null(m$coefficients)) {
    c(est = unname(m$coefficients[term]), se = unname(m$std.error[term]))
  } else {
    c(est = unname(m$coef[term]), se = unname(m$se[term]))
  }
}

mod_nobs <- function(m) m$nobs
mod_nc   <- function(m) { a <- attr(m, "n_country"); if (!is.null(a)) a else m$n_country }
mod_r2   <- function(m) if (!is.null(m$r.squared)) unname(m$r.squared) else NA_real_

# wald
wald_stat <- function(m, terms, n_clusters) {
  b <- m$coefficients[terms]
  V <- m$vcov[terms, terms, drop = FALSE]
  q <- length(terms)
  Fstat <- as.numeric(t(b) %*% solve(V) %*% b) / q
  list(F = Fstat, p = pf(Fstat, q, n_clusters - 1, lower.tail = FALSE), q = q)
}

wald_p <- function(m, terms, n_clusters) wald_stat(m, terms, n_clusters)$p

# effects
FE_TERMS  <- c("factor(year)", "factor(code)")
FE_YEAR   <- FE_TERMS[1]

fit_ols <- function(df, lhs, rhs, country_fe, cluster = "code") {
  terms <- c(rhs, FE_YEAR, if (country_fe) FE_TERMS[2])
  f <- reformulate(terms, response = lhs)
  d <- complete_on(df, c(lhs, rhs, cluster))
  m <- lm_robust(f, data = d, clusters = d[[cluster]], se_type = "stata")
  attr(m, "n_country") <- n_distinct(d$code)
  m
}

fit_iv <- function(df, lhs, endog, inst, exog = character(), country_fe = TRUE,
                   cluster = "code") {
  dummies <- c(FE_YEAR, if (country_fe) FE_TERMS[2])
  rhs  <- paste(c(endog, exog, dummies), collapse = " + ")
  inst_rhs <- paste(c(inst, exog, dummies), collapse = " + ")
  f <- as.formula(paste(lhs, "~", rhs, "|", inst_rhs))
  d <- complete_on(df, c(lhs, endog, inst, exog, cluster))
  m <- iv_robust(f, data = d, clusters = d[[cluster]], se_type = "stata")
  attr(m, "n_country") <- n_distinct(d$code)
  m
}

fit_first_stage <- function(df, endog, inst, exog = character(), country_fe = TRUE,
                            cluster = "code") {
  fit_ols(df, lhs = endog, rhs = c(inst, exog), country_fe = country_fe, cluster = cluster)
}

# rank
drop_collinear <- function(M) {
  q <- qr(M)
  M[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
}

# schemes
AGG_SCHEMES <- c("none", "lag", "period", "period_mean", "full")
AGG_FAMILIES <- c("block", "period_geom")

SCHEME_LABEL <- c(
  none = "One per lag and period (uncollapsed)",
  lag = "One per lag distance (Roodman)",
  period = "One per period, summed",
  period_mean = "One per period, averaged",
  full = "One per variable, all lags collapsed")

# sweeps
LAG_WINDOW  <- 2:8
KEEP_LABEL  <- c(none = "uncollapsed", lag = "collapsed")
GMM_MODELS  <- c(onestep = "One-step", twostep = "Two-step")
OVERID_TEST <- c(onestep = "Sargan", twostep = "Hansen")
GEOM_RHO    <- seq(0.1, 1, by = 0.1)
SUBSET_DRAWS <- 200L
SUBSET_SEED  <- 20260811L

# inference
CI_LEVEL   <- 0.05
AR_GRID    <- seq(-6, 6, by = 0.002)
CLR_GRID   <- seq(-6, 6, by = 0.005)
CLR_SIM    <- 400000L
CLR_SEED   <- 20260728L

# thresholds
STOCK_YOGO_10 <- list(`1` = c(`1` = 16.38, `2` = 19.93, `3` = 22.30, `4` = 24.58),
                      `2` = c(`2` = 7.03))
sy10 <- function(n_endog, k)
  unname(STOCK_YOGO_10[[as.character(n_endog)]][as.character(k)])
MOP_CV_10 <- 23.109

# citations
CITE_ARELLANO_OPT <- "Arellano (2003b)"

TEST_LABEL <- c(
  ar = "Anderson-Rubin (clustered)",
  clr_robust = "Moreira CLR, cluster-robust (Kleibergen)",
  clr_hom = "Moreira CLR (homoskedastic)")

# instruments
gmm_instruments <- function(df_full, est, dep_level, group = "code", period = "period",
                            lag_start = 2L, lag_max = Inf, scheme = "none",
                            block_size = 1L, rho = 1) {
  stopifnot(scheme %in% c(AGG_SCHEMES, AGG_FAMILIES), block_size >= 1L, rho > 0)
  key <- paste(df_full[[group]], df_full[[period]], sep = "\r")
  periods <- sort(unique(df_full[[period]]))
  p0 <- min(periods)
  cols <- list(); hits <- list()
  for (vname in dep_level) {
    lev <- df_full[[vname]]
    for (p in periods) for (q in periods[periods <= p - lag_start & periods >= p - lag_max]) {
      src <- lev[match(paste(est[[group]], q, sep = "\r"), key)]
      avail <- est[[period]] == p & !is.na(src)
      v <- ifelse(avail, src, 0)
      if (scheme == "period_geom") v <- v * rho^(p - q - lag_start)
      if (!any(v != 0)) next
      k <- paste0(vname, "|", switch(scheme,
                  none = paste0("g", p, "_", q),
                  lag  = paste0("l", p - q),
                  period = paste0("p", p),
                  period_mean = paste0("p", p),
                  period_geom = paste0("p", p),
                  full = "all",
                  block = paste0("l", p - q, "b", ceiling((p - p0 + 1) / block_size))))
      if (is.null(cols[[k]])) { cols[[k]] <- v; hits[[k]] <- avail + 0 }
      else { cols[[k]] <- cols[[k]] + v; hits[[k]] <- hits[[k]] + avail }
    }
  }
  if (scheme == "period_mean")
    for (k in names(cols)) cols[[k]] <- ifelse(hits[[k]] > 0, cols[[k]] / pmax(hits[[k]], 1), 0)
  do.call(cbind, cols)
}

# gmm
fit_abgmm <- function(df_full, est, dep_level, endog, exog = character(),
                      inst_extra = character(), group = "code", period = "period",
                      yearvar = "year", prevyear = "year_l1", gmm_lag_start = 2L,
                      lag_max = Inf, scheme = "none", block_size = 1L, rho = 1,
                      keep_gmm = NULL, model = "onestep") {
  stopifnot(model %in% names(GMM_MODELS))
  est <- est[order(est[[group]], est[[period]]), ]

  yrs <- sort(unique(est[[yearvar]]))
  yd <- vapply(yrs, function(tt) as.integer(est[[yearvar]] == tt) -
                 as.integer(est[[prevyear]] == tt), integer(nrow(est)))
  yd <- matrix(yd, nrow = nrow(est)); colnames(yd) <- paste0("dyr", yrs)

  X <- drop_collinear(cbind(as.matrix(est[, c(endog, exog), drop = FALSE]), yd))

  Zgmm <- gmm_instruments(df_full, est, dep_level, group, period,
                          lag_start = gmm_lag_start, lag_max = lag_max, scheme = scheme,
                          block_size = block_size, rho = rho)
  if (!is.null(keep_gmm)) Zgmm <- Zgmm[, keep_gmm, drop = FALSE]
  Z <- drop_collinear(cbind(Zgmm, yd,
         if (length(exog))       as.matrix(est[, exog, drop = FALSE]),
         if (length(inst_extra)) as.matrix(est[, inst_extra, drop = FALSE])))

  grp <- est[[group]]; pe <- est[[period]]
  A <- matrix(0, ncol(Z), ncol(Z))
  for (g in unique(grp)) {
    idx <- which(grp == g); Ti <- length(idx); pp <- pe[idx]
    H <- matrix(0, Ti, Ti); diag(H) <- 2
    for (j in seq_len(Ti)) for (k in seq_len(Ti)) if (abs(pp[j] - pp[k]) == 1) H[j, k] <- -1
    Zi <- Z[idx, , drop = FALSE]; A <- A + crossprod(Zi, H %*% Zi)
  }
  ZtX <- crossprod(Z, X); Zty <- crossprod(Z, est$y)
  solve_gmm <- function(W) {
    bread <- solve(t(ZtX) %*% W %*% ZtX)
    b <- as.numeric(bread %*% (t(ZtX) %*% W %*% Zty)); names(b) <- colnames(X)
    list(bread = bread, beta = b, e = as.numeric(est$y - X %*% b))
  }
  meat <- function(e) {
    M <- matrix(0, ncol(Z), ncol(Z))
    for (g in unique(grp)) {
      idx <- which(grp == g); gi <- crossprod(Z[idx, , drop = FALSE], e[idx])
      M <- M + gi %*% t(gi)
    }
    M
  }
  W1 <- MASS::ginv(A)
  s1 <- solve_gmm(W1)
  M1 <- meat(s1$e)
  if (model == "onestep") {
    fit <- s1
    V <- s1$bread %*% (t(ZtX) %*% W1 %*% M1 %*% W1 %*% ZtX) %*% s1$bread
  } else {
    W2 <- MASS::ginv(M1)
    fit <- solve_gmm(W2)
    V <- fit$bread
  }
  W2 <- MASS::ginv(M1)
  gbar <- crossprod(Z, fit$e)
  df_j <- ncol(Z) - ncol(X)
  hansen <- as.numeric(t(gbar) %*% W2 %*% gbar)
  list(coef = fit$beta, se = sqrt(diag(V)), nobs = length(est$y),
       n_country = length(unique(grp)), n_inst = ncol(Z), n_gmm = ncol(Zgmm),
       n_par = ncol(X), model = model, hansen = hansen, hansen_df = df_j,
       hansen_p = if (df_j > 0L) pchisq(hansen, df_j, lower.tail = FALSE) else NA_real_)
}

# rows
LBL <- list(
  dem = "Democracy_t-1", inc = "Log GDP per capita_t-1",
  obs = "Observations", ctry = "Countries", r2 = "R-squared",
  fsr2 = "First-stage R-squared", logpop = "Log population_t-1",
  educ = "Education_t-1", labor = "Labor share_t-1",
  twdem_t = "Trade-weighted democracy_t", twdem = "Trade-weighted democracy",
  sav2 = "Savings rate_t-2", sav3 = "Savings rate_t-3",
  twgdp1 = "Trade-weighted log GDP_t-1", twgdp2 = "Trade-weighted log GDP_t-2")
fstg <- function(x) paste("First stage:", x)

# measures
MEASURES <- list(
  list(dep = "fhpolrigaug", inc = "lrgdpch", label = "Freedom House"),
  list(dep = "polity4",     inc = "lrgdpch", label = "Polity"))

ROWS_DYNAMIC <- c(LBL$dem, LBL$inc, LBL$obs, LBL$ctry, LBL$r2)

# format
num <- function(x, d = 3) ifelse(is.na(x), "", sprintf(paste0("%.", d, "f"), x))

NUMWORD <- c("one", "two", "three", "four", "five", "six", "seven", "eight",
             "nine", "ten", "eleven", "twelve")
spell <- function(n) if (n >= 1L && n <= length(NUMWORD)) NUMWORD[n] else as.character(n)
sentence_case <- function(x) sub("^(.)", "\\U\\1", x, perl = TRUE)
commas <- function(x, d = NULL) paste(
  if (!is.numeric(x)) x
  else if (is.null(d)) format(x, trim = TRUE, drop0trailing = TRUE)
  else num(x, d), collapse = ", ")

# wrapping
MD_WIDTH <- 78

md_doc <- function(...) {
  blocks <- unlist(list(...), use.names = FALSE)
  out <- character(0)
  for (b in blocks) {
    if (!nzchar(b)) { out <- c(out, ""); next }
    if (grepl("^#", b)) { out <- c(out, b, ""); next }
    if (grepl("^(\\||-{2,}|[ ]{2})", b)) { out <- c(out, b); next }
    out <- c(out, strwrap(b, width = MD_WIDTH), "")
  }
  out <- out[c(TRUE, head(out, -1) != "" | out[-1] != "")]
  while (length(out) && !nzchar(out[length(out)])) out <- head(out, -1)
  out
}

write_doc <- function(file, ...) writeLines(md_doc(...), file.path(PATH_DOCS, file))

# plm
gmm_panel <- function(s, dep, inc)
  pdata.frame(s[stats::complete.cases(s[, c(dep, inc)]), ], index = c("code", "year_numeric"))
gmm_formula <- function(dep, inc, lags)
  as.formula(sprintf("%s ~ lag(%s, 1) + lag(%s, 1) | lag(%s, %s) + lag(%s, %s)",
                     dep, dep, inc, dep, lags, inc, lags))
pgmm_co <- function(m) summary(m, robust = TRUE)$coefficients
co_row  <- function(co, var) co[grep(paste0("lag\\(", var), rownames(co))[1], 1:2]
pgmm_countries <- function(m)
  sum(vapply(m$residuals, function(r) any(r != 0), logical(1)))

# builder
table_builder <- function(dep, ...) {
  rows <- list()
  extra <- list(...)
  push <- function(column, row, value, se = NA_real_, type = "coef")
    rows[[length(rows) + 1]] <<- do.call(tibble, c(extra,
      list(column = column, row = row,
           value = as.numeric(value), se = as.numeric(se), type = type)))
  counts <- function(col, m, fs = NULL, r2 = FALSE) {
    push(col, LBL$obs, mod_nobs(m), type = "count")
    push(col, LBL$ctry, mod_nc(m), type = "count")
    if (r2) push(col, LBL$r2, mod_r2(m), type = "r2")
    if (!is.null(fs)) push(col, LBL$fsr2, mod_r2(fs), type = "r2")
  }
  sls <- function(col, dat, inst, exog = character(), endog = "Linc",
                  second = c(), fs = c()) {
    m <- fit_iv(dat, dep, endog = endog, inst = inst, exog = exog, country_fe = TRUE)
    f <- fit_first_stage(dat, endog = endog, inst = c(inst, exog), country_fe = TRUE)
    for (lab in names(second)) push(col, lab, ce(m, second[[lab]])["est"], ce(m, second[[lab]])["se"])
    for (lab in names(fs))     push(col, lab, ce(f, fs[[lab]])["est"], ce(f, fs[[lab]])["se"])
    counts(col, m, f); invisible(NULL)
  }
  list(push = push, counts = counts, sls = sls, collect = function() bind_rows(rows))
}

source(here::here("R", "_engine.R"))
source(here::here("R", "_weakiv.R"))

# weakiv
DEP_WEAKIV <- MEASURES[[1]]$dep

# specs
TSLS_SPECS <- list(
  list(tab = "5", col = 4, panel = "savings",     inst = "z2",            exog = character()),
  list(tab = "5", col = 5, panel = "savings",     inst = "z2",            exog = "Ldep"),
  list(tab = "5", col = 7, panel = "savings",     inst = "z2",            exog = "Llabor"),
  list(tab = "5", col = 8, panel = "savings",     inst = "z2",
       exog = c("Ldep", "L2dep", "L3dep")),
  list(tab = "5", col = 9, panel = "savings",     inst = c("z2", "z3"),   exog = character()),
  list(tab = "6", col = 4, panel = "worldincome", inst = "z1",            exog = character()),
  list(tab = "6", col = 5, panel = "worldincome", inst = "z1",            exog = "Ldep"),
  list(tab = "6", col = 7, panel = "worldincome", inst = "z1",            exog = "wdem"),
  list(tab = "6", col = 8, panel = "worldincome", inst = "z2wi",          exog = character()),
  list(tab = "6", col = 9, panel = "worldincome", inst = c("z1", "z2wi"), exog = character()))

AH_SPEC <- list(tab = c("2", "3"), col = 3, endog = c("dLdep", "dLinc"),
                inst = c("L2dep", "L2inc"))

spec_label <- function(tab, col) sprintf("Table %s col %d", tab, col)

tsls_panel <- local({
  builders <- list(savings = prep_savings_panel, worldincome = prep_worldincome_panel)
  cache <- list()
  function(name) {
    if (is.null(cache[[name]])) cache[[name]] <<- filter(builders[[name]](), sample == 1)
    cache[[name]]
  }
})

spec_data <- function(sp, dep = DEP_WEAKIV)
  complete_on(tsls_panel(sp$panel), c(dep, "Linc", sp$inst, sp$exog, "code"))

# published
published <- local({
  cache <- NULL
  function(...) {
    if (is.null(cache))
      cache <<- read_csv(file.path(PATH_DOCS, "published_values.csv"), show_col_types = FALSE)
    filter(cache, ...)
  }
})
