suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(readxl)
  library(estimatr)
  library(here)
})

here::i_am("R/00_setup.R")

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

COLS_5YR <- c("code", "country", "code_numeric", "year", "year_numeric",
              "sample", "samplebalancefe", "samplebalancegmm", "socialist",
              "fhpolrigaug", "polity4", "lrgdpch", "nsave", "worldincome",
              "worlddemocracy", "laborshare", "lpop", "medage", "education",
              "age_veryyoung", "age_young", "age_midage", "age_old")
COLS_FREQ <- c("code", "country", "code_numeric", "year", "year_numeric",
               "sample", "fhpolrigaug", "polity4", "lrgdpch")
COLS_MAD  <- c("code", "country", "code_numeric", "year", "year_numeric",
               "sample", "noextrapolation", "madid", "polity4", "lrgdpmad")

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

add_lags <- function(df, vars, ks, id = "code", time = "period") {
  stopifnot(!anyDuplicated(df[c(id, time)]))
  key <- paste(df[[id]], df[[time]], sep = "\r")
  for (v in vars) for (k in ks) {
    src_key <- paste(df[[id]], df[[time]] + k, sep = "\r")
    df[[paste0(v, "_l", k)]] <- df[[v]][match(key, src_key)]
  }
  df
}

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

wald_p <- function(m, terms, n_clusters) {
  b <- m$coefficients[terms]
  V <- m$vcov[terms, terms, drop = FALSE]
  q <- length(terms)
  Fstat <- as.numeric(t(b) %*% solve(V) %*% b) / q
  pf(Fstat, q, n_clusters - 1, lower.tail = FALSE)
}

fit_ols <- function(df, lhs, rhs, country_fe, cluster = "code") {
  terms <- c(rhs, "factor(year)", if (country_fe) "factor(code)")
  f <- reformulate(terms, response = lhs)
  d <- complete_on(df, c(lhs, rhs, cluster))
  m <- lm_robust(f, data = d, clusters = d[[cluster]], se_type = "stata")
  attr(m, "n_country") <- n_distinct(d$code)
  m
}

fit_iv <- function(df, lhs, endog, inst, exog = character(), country_fe = TRUE,
                   cluster = "code") {
  dummies <- c("factor(year)", if (country_fe) "factor(code)")
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

fit_abgmm <- function(df_full, est, dep_level, endog, exog = character(),
                      inst_extra = character(), group = "code", period = "period",
                      yearvar = "year", prevyear = "year_l1", gmm_lag_start = 2L) {
  est <- est[order(est[[group]], est[[period]]), ]

  yrs <- sort(unique(est[[yearvar]]))
  yd <- sapply(yrs, function(tt) as.integer(est[[yearvar]] == tt) -
                 as.integer(est[[prevyear]] == tt))
  yd <- matrix(yd, nrow = nrow(est)); colnames(yd) <- paste0("dyr", yrs)

  Xfull <- cbind(as.matrix(est[, c(endog, exog), drop = FALSE]), yd)
  X <- Xfull[, sort(qr(Xfull)$pivot[seq_len(qr(Xfull)$rank)]), drop = FALSE]

  key <- paste(df_full[[group]], df_full[[period]])
  lev <- df_full[[dep_level]]
  periods <- sort(unique(df_full[[period]]))
  gmm_cols <- list()
  for (p in periods) for (q in periods[periods <= p - gmm_lag_start]) {
    col <- ifelse(est[[period]] == p, lev[match(paste(est[[group]], q), key)], 0)
    col[is.na(col)] <- 0
    if (any(col != 0)) gmm_cols[[paste0("g", p, "_", q)]] <- col
  }
  Zgmm <- do.call(cbind, gmm_cols)
  Zfull <- cbind(Zgmm, yd,
                 if (length(exog))       as.matrix(est[, exog, drop = FALSE]),
                 if (length(inst_extra)) as.matrix(est[, inst_extra, drop = FALSE]))
  Z <- Zfull[, sort(qr(Zfull)$pivot[seq_len(qr(Zfull)$rank)]), drop = FALSE]

  grp <- est[[group]]; pe <- est[[period]]
  A <- matrix(0, ncol(Z), ncol(Z))
  for (g in unique(grp)) {
    idx <- which(grp == g); Ti <- length(idx); pp <- pe[idx]
    H <- matrix(0, Ti, Ti); diag(H) <- 2
    for (j in seq_len(Ti)) for (k in seq_len(Ti)) if (abs(pp[j] - pp[k]) == 1) H[j, k] <- -1
    Zi <- Z[idx, , drop = FALSE]; A <- A + crossprod(Zi, H %*% Zi)
  }
  W <- MASS::ginv(A)
  ZtX <- crossprod(Z, X); Zty <- crossprod(Z, est$y)
  bread <- solve(t(ZtX) %*% W %*% ZtX)
  beta <- as.numeric(bread %*% (t(ZtX) %*% W %*% Zty)); names(beta) <- colnames(X)
  e <- as.numeric(est$y - X %*% beta)
  M <- matrix(0, ncol(Z), ncol(Z))
  for (g in unique(grp)) {
    idx <- which(grp == g); gi <- crossprod(Z[idx, , drop = FALSE], e[idx])
    M <- M + gi %*% t(gi)
  }
  V <- bread %*% (t(ZtX) %*% W %*% M %*% W %*% ZtX) %*% bread
  list(coef = beta, se = sqrt(diag(V)), nobs = length(est$y),
       n_country = length(unique(grp)), n_inst = ncol(Z))
}

source(here::here("R", "_engine.R"))
