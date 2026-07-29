source(here::here("R", "00_setup.R"))

DEP <- "fhpolrigaug"
AR_GRID <- seq(-6, 6, by = 0.002)
CLR_GRID <- seq(-6, 6, by = 0.005)
CLR_SIM <- 400000L
LEVEL <- 0.05
STOCK_YOGO <- c(`1` = 16.38, `2` = 19.93, `3` = 22.30, `4` = 24.58)

partial_out <- function(dat, cols, exog) {
  W <- model.matrix(reformulate(c("factor(year)", "factor(code)", exog)), data = dat)
  q <- qr(W)
  list(R = qr.resid(q, as.matrix(dat[, cols, drop = FALSE])), rank = q$rank)
}

ar_pcurve <- function(dat, inst, exog, grid) {
  po <- partial_out(dat, c(DEP, "Linc", inst), exog)
  yt <- po$R[, DEP]; xt <- po$R[, "Linc"]; Zt <- po$R[, inst, drop = FALSE]
  g <- dat$code
  G <- n_distinct(g); N <- nrow(dat); k <- length(inst)
  cst <- (G / (G - 1)) * ((N - 1) / (N - po$rank - k))
  ZZi <- solve(crossprod(Zt))
  uy <- as.numeric(yt - Zt %*% (ZZi %*% crossprod(Zt, yt)))
  ux <- as.numeric(xt - Zt %*% (ZZi %*% crossprod(Zt, xt)))
  P  <- rowsum(Zt * yt, g);  Q  <- rowsum(Zt * xt, g)
  Ps <- rowsum(Zt * uy, g);  Qs <- rowsum(Zt * ux, g)
  fstat <- vapply(grid, function(b0) {
    s <- colSums(P - b0 * Q)
    as.numeric(crossprod(s, solve(crossprod(Ps - b0 * Qs), s))) / (k * cst)
  }, numeric(1))
  list(p = pf(fstat, k, G - 1, lower.tail = FALSE), stat = fstat, k = k, G = G,
       rank = po$rank, yt = yt, xt = xt, Zt = Zt, g = g)
}

clr_pcurve <- function(dat, inst, exog, grid, nsim = CLR_SIM, seed = 20260728L) {
  po <- partial_out(dat, c(DEP, "Linc", inst), exog)
  Y <- cbind(po$R[, DEP], po$R[, "Linc"]); Zt <- po$R[, inst, drop = FALSE]
  k <- length(inst); dfree <- nrow(Y) - po$rank - k
  ZZ <- crossprod(Zt); ZY <- crossprod(Zt, Y)
  Om <- crossprod(Y - Zt %*% solve(ZZ, ZY)) / dfree
  ei <- eigen(ZZ, symmetric = TRUE)
  ZZmh <- ei$vectors %*% diag(1 / sqrt(ei$values), k) %*% t(ei$vectors)
  Omi <- solve(Om)
  if (k > 1L) {
    old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv), add = TRUE)
    set.seed(seed)
    z1 <- rnorm(nsim); w <- rchisq(nsim, k - 1)
    QSs <- z1^2 + w
  }
  vapply(grid, function(b0) {
    b <- c(1, -b0); a <- c(b0, 1)
    S <- ZZmh %*% (ZY %*% b) / sqrt(as.numeric(t(b) %*% Om %*% b))
    Tt <- ZZmh %*% (ZY %*% (Omi %*% a)) / sqrt(as.numeric(t(a) %*% Omi %*% a))
    QS <- sum(S^2); QT <- sum(Tt^2); QST <- sum(S * Tt)
    LR <- 0.5 * (QS - QT + sqrt(max((QS + QT)^2 - 4 * (QS * QT - QST^2), 0)))
    if (k == 1L) return(pchisq(LR, 1, lower.tail = FALSE))
    LRs <- 0.5 * (QSs - QT + sqrt(pmax((QSs + QT)^2 - 4 * QT * w, 0)))
    mean(LRs >= LR)
  }, numeric(1))
}

set_info <- function(grid, p, level = LEVEL) {
  keep <- p >= level
  if (!any(keep)) return(list(text = "empty", lo = NA_real_, hi = NA_real_,
                              pieces = 0L, bounded = TRUE))
  r <- rle(keep)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  idx <- which(r$values)
  first <- starts[idx[1]]; last <- ends[idx[length(idx)]]
  list(text = paste(vapply(idx, function(i) {
         lo <- if (starts[i] == 1L) "-inf" else sprintf("%.2f", grid[starts[i]])
         hi <- if (ends[i] == length(grid)) "+inf" else sprintf("%.2f", grid[ends[i]])
         sprintf("[%s, %s]", lo, hi)
       }, character(1)), collapse = " U "),
       lo = grid[first], hi = grid[last], pieces = length(idx),
       bounded = first > 1L && last < length(grid))
}

s5 <- filter(prep_savings_panel(), sample == 1)
s6 <- filter(prep_worldincome_panel(), sample == 1)

SPECS <- list(
  list(tab = "5", col = 4, d = s5, inst = "z2",           exog = character()),
  list(tab = "5", col = 5, d = s5, inst = "z2",           exog = "Ldep"),
  list(tab = "5", col = 7, d = s5, inst = "z2",           exog = "Llabor"),
  list(tab = "5", col = 8, d = s5, inst = "z2",           exog = c("Ldep", "L2dep", "L3dep")),
  list(tab = "5", col = 9, d = s5, inst = c("z2", "z3"),  exog = character()),
  list(tab = "6", col = 4, d = s6, inst = "z1",           exog = character()),
  list(tab = "6", col = 5, d = s6, inst = "z1",           exog = "Ldep"),
  list(tab = "6", col = 7, d = s6, inst = "z1",           exog = "wdem"),
  list(tab = "6", col = 8, d = s6, inst = "z2wi",         exog = character()),
  list(tab = "6", col = 9, d = s6, inst = c("z1", "z2wi"), exog = character()))

pub <- read_csv(file.path(PATH_DOCS, "published_values.csv"), show_col_types = FALSE) |>
  filter(table %in% c("5", "6"), row == "inc", type == "coef")

rows <- list(); curves <- list()
for (sp in SPECS) {
  dat <- complete_on(sp$d, c(DEP, "Linc", sp$inst, sp$exog, "code"))
  m <- fit_iv(dat, DEP, endog = "Linc", inst = sp$inst, exog = sp$exog, country_fe = TRUE)
  fs <- fit_first_stage(dat, endog = "Linc", inst = c(sp$inst, sp$exog), country_fe = TRUE)
  G <- mod_nc(m)
  est <- ce(m, "Linc")
  tcrit <- qt(1 - LEVEL / 2, G - 1)

  want <- filter(pub, table == sp$tab, column == sp$col)$value
  stopifnot(length(want) == 1L, abs(round(est["est"], 3) - want) < 1e-9)

  fF <- wald_stat(fs, sp$inst, G)
  ar <- ar_pcurve(dat, sp$inst, sp$exog, AR_GRID)
  clr_p <- clr_pcurve(dat, sp$inst, sp$exog, CLR_GRID)

  if (length(sp$inst) == 1L) {
    at_hat <- which.min(abs(AR_GRID - est["est"]))
    stopifnot(ar$p[at_hat] > 0.99)
  }

  lab <- sprintf("Table %s col %d", sp$tab, sp$col)
  ai <- set_info(AR_GRID, ar$p); ci <- set_info(CLR_GRID, clr_p)
  rows[[length(rows) + 1]] <- tibble(
    spec = lab, table = sp$tab, column = sp$col,
    instruments = paste(sp$inst, collapse = " + "), k = length(sp$inst),
    countries = G, obs = mod_nobs(m),
    first_stage_F = fF$F,
    stock_yogo = unname(STOCK_YOGO[as.character(length(sp$inst))]),
    tsls = unname(est["est"]), tsls_se = unname(est["se"]),
    tsls_lo = unname(est["est"] - tcrit * est["se"]),
    tsls_hi = unname(est["est"] + tcrit * est["se"]),
    ar_set = ai$text, ar_lo = ai$lo, ar_hi = ai$hi, ar_bounded = ai$bounded,
    ar_p0 = ar$p[which.min(abs(AR_GRID))],
    clr_set = ci$text, clr_lo = ci$lo, clr_hi = ci$hi,
    clr_p0 = clr_p[which.min(abs(CLR_GRID))])
  curves[[length(curves) + 1]] <- bind_rows(
    tibble(spec = lab, beta = AR_GRID, p = ar$p, test = "Anderson-Rubin (clustered)"),
    tibble(spec = lab, beta = CLR_GRID, p = clr_p, test = "Moreira CLR (homoskedastic)"))
}

tab <- bind_rows(rows)
curve <- bind_rows(curves) |> mutate(spec = factor(spec, levels = tab$spec))
write_csv(tab, file.path(PATH_OUTPUT, "weakiv.csv"))

check_grid <- c(-0.5, 0, 0.25)
for (sp in SPECS[c(1, 5, 8, 10)]) {
  dat <- complete_on(sp$d, c(DEP, "Linc", sp$inst, sp$exog, "code"))
  a <- ar_pcurve(dat, sp$inst, sp$exog, check_grid)

  px <- a$Zt %*% solve(crossprod(a$Zt), crossprod(a$Zt, a$xt))
  b_hand <- as.numeric(crossprod(px, a$yt) / crossprod(px, a$xt))
  b_fit <- unname(ce(fit_iv(dat, DEP, endog = "Linc", inst = sp$inst,
                            exog = sp$exog, country_fe = TRUE), "Linc")["est"])
  stopifnot(abs(b_hand - b_fit) < 1e-9)

  slow <- vapply(check_grid, function(b0) {
    d2 <- dat; d2$e0 <- d2[[DEP]] - b0 * d2$Linc
    mm <- fit_ols(d2, "e0", c(sp$inst, sp$exog), country_fe = TRUE)
    wald_stat(mm, sp$inst, mod_nc(mm))$p
  }, numeric(1))
  stopifnot(max(abs(slow - a$p)) < 1e-8)

  if (length(sp$inst) == 1L) {
    cp <- clr_pcurve(dat, sp$inst, sp$exog, check_grid)
    dfree <- length(a$yt) - a$rank - 1L
    hom <- vapply(check_grid, function(b0) {
      e <- a$yt - b0 * a$xt
      fitz <- a$Zt %*% solve(crossprod(a$Zt), crossprod(a$Zt, e))
      pchisq(sum(fitz^2) / (sum((e - fitz)^2) / dfree), 1, lower.tail = FALSE)
    }, numeric(1))
    stopifnot(max(abs(cp - hom)) < 1e-12)
  }
}
cat("Partialling reproduces the 2SLS estimates; the fast Anderson-Rubin curve matches",
    "the clustered regression it stands in for;\nCLR collapses to Anderson-Rubin when",
    "the model is just identified.\n")

txt <- c(sprintf("%-16s %-14s %3s %9s %9s %-20s %-22s %-22s",
                 "Spec", "Instrument(s)", "k", "1st-st F", "SY 10%",
                 "2SLS (SE)", "Conventional 95% CI", "Anderson-Rubin 95%"))
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  txt <- c(txt, sprintf("%-16s %-14s %3d %9.2f %9.2f %-20s %-22s %-22s",
    r$spec, r$instruments, r$k, r$first_stage_F, r$stock_yogo,
    sprintf("%s (%s)", num(r$tsls), num(r$tsls_se)),
    sprintf("[%.2f, %.2f]", r$tsls_lo, r$tsls_hi), r$ar_set))
}
txt <- c(txt, "", "Moreira CLR 95% sets (homoskedastic):")
for (i in seq_len(nrow(tab))) txt <- c(txt, sprintf("  %-16s %s", tab$spec[i], tab$clr_set[i]))
writeLines(txt, file.path(PATH_OUTPUT, "weakiv.txt"))

fig <- ggplot(curve, aes(beta, p, color = test)) +
  geom_hline(yintercept = LEVEL, color = "grey40", linetype = "dashed") +
  geom_line(linewidth = 0.5) +
  geom_vline(data = tab, aes(xintercept = tsls), color = "grey30",
             linetype = "dotted", linewidth = 0.4) +
  facet_wrap(~spec, ncol = 5) +
  coord_cartesian(xlim = c(-1, 0.6), ylim = c(0, 1)) +
  labs(x = "Value of the income coefficient being tested", y = "p-value",
       color = NULL,
       title = "Weak-instrument-robust p-values for the effect of income",
       subtitle = "Dashed line is 5%; anything above it stays in the confidence set. Dotted line is the 2SLS estimate.") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")
ggsave(file.path(PATH_OUTPUT, "weakiv.png"), fig, width = 11, height = 5.5, dpi = 150)

n_below <- sum(tab$first_stage_F < tab$stock_yogo)
n_zero <- sum(tab$ar_p0 >= LEVEL)
widen <- tab |>
  mutate(conv = tsls_hi - tsls_lo, ar = ar_hi - ar_lo, ratio = ar / conv) |>
  arrange(desc(ratio))
ord <- arrange(widen, first_stage_F)
inv <- NULL
for (i in seq_len(nrow(ord) - 1L)) for (j in (i + 1L):nrow(ord))
  if (is.null(inv) && ord$ratio[i] < ord$ratio[j]) inv <- c(i, j)
strongest <- tab[which.max(tab$first_stage_F), ]
weakest <- tab[which.min(tab$first_stage_F), ]

md <- c(
"# Are the instruments strong enough?",
"",
"Professor Torgovitsky asked for weak instrument diagnostics, then for",
"Anderson-Rubin confidence sets and the Moreira conditional likelihood ratio",
"approach. All three apply to the two-stage least squares columns of Tables 5",
"and 6, which are the parts of the paper with one variable being instrumented,",
"income, and a small number of instruments for it: the savings rate in Table 5",
"and trade-weighted world income in Table 6. That gives ten columns to check.",
"The numbers are in output/weakiv.csv and output/weakiv.txt, and the p-value",
"curves are drawn in output/weakiv.png. The two tests are from Anderson and",
"Rubin (1949) and Moreira (2003), the first-stage cutoffs from Stock and Yogo",
"(2005), and full details are in docs/references.md.",
"",
"## Why this matters",
"",
"An instrument is weak when it barely moves the variable it is supposed to be",
"standing in for. When that happens the usual standard errors and confidence",
"intervals are not trustworthy, because they are built on the assumption that",
"the instrument has real pull. The tests below are built to stay honest whether",
"the instrument is strong or weak, so comparing them against the ordinary",
"interval shows how much the ordinary one was relying on that assumption.",
"",
"## How strong is the first stage",
"",
sprintf("The first-stage F statistics run from %.1f to %.1f, clustered by country the",
        min(tab$first_stage_F), max(tab$first_stage_F)),
"same way the paper clusters everything else. Every column clears the old rule",
"of thumb that an F above 10 is good enough. Against the stricter Stock and",
sprintf("Yogo threshold, %d of the %d columns fall short. So the instruments are not",
        n_below, nrow(tab)),
"badly weak, but a few are not comfortably strong either, which is exactly the",
"situation the next two tests are designed for.",
"",
"One caveat on that comparison. The Stock and Yogo cutoffs were worked out for a",
"setting with well behaved, unclustered errors, and the F statistics here are",
"clustered by country. Lining the two up is what most applied work does, and it",
"is what the standard software prints, but it is not strictly justified. I use",
"the threshold as a rough marker rather than a real test, which is another",
"reason to lean on the confidence sets below instead. The properly justified",
"cutoff for clustered errors would be the effective F statistic of Montiel Olea",
"and Pflueger (2013), which I have not implemented.",
"",
"## Anderson-Rubin",
"",
"The Anderson-Rubin test asks a simple question. Pick a candidate value for the",
"effect of income, subtract off what that value would imply, and see whether",
"anything is left over that the instruments can still explain. If something is,",
"that candidate value is rejected. Collecting every value that survives gives a",
"confidence set, and it stays valid no matter how weak the instruments are. I",
"clustered it by country to match the paper.",
"",
sprintf("%s %d sets come out bounded, and %s contain zero. Where the first stage",
        ifelse(all(tab$ar_bounded), "All", paste(sum(tab$ar_bounded), "of the")),
        nrow(tab),
        ifelse(n_zero == nrow(tab), "all of them", paste(n_zero, "of them"))),
sprintf("is strongest (%s, F of %.1f) the set is %s, which is",
        strongest$spec, strongest$first_stage_F, strongest$ar_set),
sprintf("nearly the same as the ordinary interval of [%.2f, %.2f]. Where it is",
        strongest$tsls_lo, strongest$tsls_hi),
sprintf("weakest (%s, F of %.1f) the set widens to %s against an",
        weakest$spec, weakest$first_stage_F, weakest$ar_set),
sprintf("ordinary interval of [%.2f, %.2f]. The Anderson-Rubin set is at least as wide",
        weakest$tsls_lo, weakest$tsls_hi),
"as the ordinary interval in every column, so the ordinary intervals are",
"optimistic throughout, not just at the extremes.",
"",
"I should not push that comparison further than it goes. Lining the ten columns",
"up by first-stage F does not give a tidy ordering. The widening ratio runs",
sprintf("from %.2f to %.2f, and it does not fall neatly as the first stage gets",
        min(widen$ratio), max(widen$ratio)),
sprintf("stronger. %s has the weaker first stage of the two,",  ord$spec[inv[1]]),
sprintf("F of %.1f against %.1f, yet widens less, %.2f against %.2f, than %s does.",
        ord$first_stage_F[inv[1]], ord$first_stage_F[inv[2]],
        ord$ratio[inv[1]], ord$ratio[inv[2]], ord$spec[inv[2]]),
"These specifications",
"differ in their controls and samples as well as in instrument strength, so the",
"two ends of the range are suggestive rather than a clean relationship.",
"",
"## Moreira's conditional likelihood ratio",
"",
"The conditional likelihood ratio test can be sharper than Anderson-Rubin.",
"Rather than comparing against a fixed cutoff, it adjusts the cutoff using a",
"quantity that carries the information about how strong the instruments are,",
"which buys back power that Anderson-Rubin gives away. I coded it directly,",
"since no package for it was installed.",
"",
sprintf("The catch is that %d of the %d columns here have exactly one instrument, and",
        sum(tab$k == 1), nrow(tab)),
"with one instrument there is nothing to buy back: the two tests are the same",
"test. So the extra power is only in play for the two columns with a pair of",
"instruments. For the single-instrument columns the cutoff is a chi-squared",
"quantile that can be written down exactly, and the code uses that rather than",
"simulating a distribution it already knows.",
"",
"There is a caveat worth stating plainly. The version of this test I implemented",
"assumes the errors are well behaved and independent, and it does not cluster by",
"country. The paper clusters, and the Anderson-Rubin sets above cluster. So the",
"conditional likelihood ratio sets come out narrower than the Anderson-Rubin",
"ones, but that gap is mostly the clustering, not the extra power of the test.",
"The code checks this: for the columns with a single instrument the conditional",
"likelihood ratio test and the unclustered Anderson-Rubin test agree, and the",
"run stops if they ever do not. The right comparison is therefore",
"Anderson-Rubin, and I treat the other as a benchmark.",
"",
"## What it adds up to",
"",
sprintf("Under inference that does not lean on the instruments being strong, %s",
        ifelse(n_zero == nrow(tab), "every one", sprintf("%d", n_zero))),
sprintf("of the %d two-stage least squares columns leaves zero inside the confidence",
        nrow(tab)),
"set. So none of them establishes an effect of income on democracy in either",
"direction.",
"",
"This does not overturn anything in the paper, whose argument is that the effect",
"is not there once the fixed differences between countries are taken out. It",
"does mean the instrumental-variables columns should be read as uninformative",
"about the sign of the effect rather than as evidence for a negative one.",
"",
"## Checks",
"",
"This is fiddly enough that I did not want to rely on code I wrote myself. Three",
"checks run every time the script does. The partialling step has to reproduce",
"the published two-stage least squares estimates, the fast Anderson-Rubin curve",
"has to match the slower clustered regression it stands in for, and the",
"conditional likelihood ratio test has to collapse onto Anderson-Rubin for the",
"columns with a single instrument. Any failure stops the run.",
"",
"On top of that, R/15_ivcrosscheck.R sends all ten columns through the ivmodel",
"package and stops unless its conditional likelihood ratio sets agree with mine",
"and its two-stage least squares estimates match exactly. They do, for all ten.",
"That script needs the ivmodel package, so it sits outside the main run; use",
"Rscript R/15_ivcrosscheck.R.")
writeLines(md, file.path(PATH_DOCS, "weak-instruments.md"))

cat("Weak-instrument diagnostics written.\n")
print(transmute(tab, spec, k, F = round(first_stage_F, 1),
                tsls = round(tsls, 3), ar_set, clr_set), n = 20)
