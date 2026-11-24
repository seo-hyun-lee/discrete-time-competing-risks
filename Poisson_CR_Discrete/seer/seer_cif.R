################################################################################
## seer_cif.R
## Cumulative incidence functions for one covariate profile, from the trend-filtering fits of
## seer_tf3.R, with pointwise 95% bootstrap confidence intervals (layout of Lee et al. 2018, Fig. 1)
##
##   BOOT=1:50 Rscript seer_cif.R boot     bootstrap replications 1-50 (run_seer_cif.sh: 4 x 50)
##   Rscript seer_cif.R figure             figure + table from the replications saved so far
##
## Bootstrap: patients resampled with replacement; for each link and cause, (alpha, beta) refitted
## by admm_bcd_exact at the (b, rho) selected on the full data, warm-started at the full-data
## estimate. (b, rho) are held fixed, so the intervals do not reflect the selection step.
## CIF (Section 2.3):    F_k(t_m; x) = sum_{l <= m} lambda_k(t_l; x) S(t_{l-1}; x),
##                       S(t_l; x) = prod_{r <= l} {1 - lambda_1(t_r; x) - lambda_2(t_r; x)}
##                       (CIF_FORM=poisson: S(t_l; x) = exp[-sum_{r <= l} {lambda_1 + lambda_2}]).
## Each replication stores (alpha, beta), so the CIF can be recomputed in either form.
## Requires seer_tf3.R (its data step is reused), seer_out3/seer_fits.rds (seer_tf3.R summary),
## admm_tf_long.R, admm_bcd_exact.R, seer_colon_wide.rds.
################################################################################
CIF_MODE <- commandArgs(TRUE)[1]; if (is.na(CIF_MODE)) CIF_MODE <- "figure"
## data step of seer_tf3.R (identical cohort and covariate coding)
src <- readLines("seer_tf3.R")
eval(parse(text = src[seq_len(grep("^## =+ 2\\. unpenalized fits", src) - 1)]))
FITS <- readRDS(file.path(OUT_DIR, "seer_fits.rds")); res <- FITS$res
BOOT_DIR <- file.path(OUT_DIR, "cif_boot2"); dir.create(BOOT_DIR, showWarnings = FALSE)
B_TOTAL <- as.integer(Sys.getenv("B", "200"))
CIF_FORM <- Sys.getenv("CIF_FORM", "exact")

## covariate profile, as in Lee et al. (2018, Fig. 1): married non-Hispanic White woman aged 66,
## grade III/IV tumour of the distal colon, diagnosed in the last year of the cohort (2003)
x0 <- c(distal = 1, grade34 = 1, single = 0, hispanic = 0, nhblack = 0, nhaian = 0, nhapi = 0,
        female = 1, age = 66 - 66, year = 2003 - 1994)
stopifnot(identical(names(x0), colnames(X)))

hazard <- function(alpha, beta, link) {
  eta <- alpha + sum(x0 * beta)
  if (link == "log") exp(eta) else plogis(eta)
}
cif <- function(l1, l2, form = CIF_FORM) {
  l1[!is.finite(l1)] <- 0; l2[!is.finite(l2)] <- 0          # unpenalized -Inf intervals
  ## CIF_FORM = "exact" (default, eq:cif of Section 2.3): S(t_{l-1}) = prod_{r<l} (1 - l1 - l2);
  ## "poisson": S(t_{l-1}) = exp{-sum_{r<l} (l1 + l2)} (kept for the comparison below)
  S_prev <- if (form == "exact") c(1, cumprod(1 - l1 - l2))[seq_along(l1)] else
            exp(-c(0, cumsum(l1 + l2))[seq_along(l1)])
  cbind(cumsum(l1 * S_prev), cumsum(l2 * S_prev))
}
cif_fit <- function(F, link, form = CIF_FORM) cif(hazard(F[[1]]$alpha, F[[1]]$beta, link),
                                 hazard(F[[2]]$alpha, F[[2]]$beta, link), form)
est <- list(log   = cif_fit(lapply(res, `[[`, "log_tf"), "log"),
            logit = cif_fit(lapply(res, `[[`, "logit_tf"), "logit"),
            lee   = cif_fit(lapply(res, `[[`, "lee"), "logit", "exact"))   # Lee: exact binomial model

## ============================================================ bootstrap
if (CIF_MODE == "boot") {
  rr <- eval(parse(text = Sys.getenv("BOOT", "1:50")))
  for (r in rr) {
    out <- file.path(BOOT_DIR, sprintf("rep_%03d.rds", r))
    if (file.exists(out)) next
    t0 <- proc.time()[3]
    set.seed(20261014 + r)
    idx <- sample(nrow(w), replace = TRUE)
    dat_b <- make_long(w$j[idx], w$event[idx], X[idx, , drop = FALSE])
    R <- list(); conv <- c()
    for (link in c("log", "logit")) {
      F <- lapply(1:2, function(k) {
        T <- res[[k]][[paste0(link, "_tf")]]
        admm_bcd_exact(dat_b, k, T$b, T$rho, link = link, cvec = GAMMA_C,
                       alpha0 = T$alpha, beta0 = T$beta, u0 = "kkt")
      })
      R[[link]] <- list(fits = lapply(F, `[`, c("alpha", "beta")),   # kept so the CIF can be
                        beta = sapply(F, `[[`, "beta"))               # recomputed in any form
      conv[link] <- all(sapply(F, `[[`, "converged"))
    }
    R$converged <- conv
    saveRDS(R, out)
    stamp("rep %3d done (%.0f s)%s", r, proc.time()[3] - t0,
          if (!all(conv)) "  NOT CONVERGED" else "")
  }
  quit(save = "no")
}

## ============================================================ figure and table
f <- sort(list.files(BOOT_DIR, "^rep_\\d+\\.rds$", full.names = TRUE))
if (!length(f)) stop("no bootstrap replications in ", BOOT_DIR)
BT <- lapply(f, readRDS)
ok <- sapply(BT, function(R) all(R$converged))
stamp("%d bootstrap replications (%d of %d planned), %d with a non-converged fit (dropped)",
      length(BT), length(BT), B_TOTAL, sum(!ok))
BT <- BT[ok]
band <- function(link, k) {
  M <- sapply(BT, function(R) cif_fit(R[[link]]$fits, link)[, k])  # J x B
  t(apply(M, 1, quantile, c(.025, .975)))
}
bands <- lapply(c(log = "log", logit = "logit"), function(link) lapply(1:2, function(k) band(link, k)))

## bootstrap SE of the trend-filtering beta (for Table tab:seer-beta, if wanted)
se_tf <- do.call(rbind, lapply(c("log", "logit"), function(link) do.call(rbind, lapply(1:2, function(k)
  data.frame(link = link, cause = k, covariate = colnames(X),
             estimate = res[[k]][[paste0(link, "_tf")]]$beta,
             se_boot = apply(sapply(BT, function(R) R[[link]]$beta[, k]), 1, sd))))))
write.csv(se_tf, file.path(OUT_DIR, "seer_tf_beta_bootse.csv"), row.names = FALSE)

## CIF at 1, 3, 5, 10 years
yrs <- c(12, 36, 60, 120)
tab <- do.call(rbind, lapply(1:2, function(k) do.call(rbind, lapply(c("log", "logit"), function(link)
  data.frame(cause = k, method = paste(link, "TF"), month = yrs, cif = est[[link]][yrs, k],
             lo = bands[[link]][[k]][yrs, 1], hi = bands[[link]][[k]][yrs, 2])))))
tab <- rbind(tab, do.call(rbind, lapply(1:2, function(k)
  data.frame(cause = k, method = "Lee", month = yrs, cif = est$lee[yrs, k], lo = NA, hi = NA))))
write.csv(tab, file.path(OUT_DIR, "seer_cif.csv"), row.names = FALSE)
print(within(tab, { cif <- round(cif, 3); lo <- round(lo, 3); hi <- round(hi, 3) }), row.names = FALSE)

pdf(file.path(OUT_DIR, "seer_cif.pdf"), width = 9, height = 4.4)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.2, 0.8), oma = c(1.6, 0, 0, 0))
tt <- 1:J
for (k in 1:2) {
  plot(NA, xlim = c(0, 120), ylim = c(0, 1), xaxt = "n", las = 1,
       xlab = "Time (Months)", ylab = "Cumulative incidence estimates",
       main = c("Death from colon cancer", "Death from other causes")[k])
  axis(1, at = seq(0, 120, 24))
  for (link in c("log", "logit")) {
    col <- if (link == "log") "black" else "red"
    lines(c(0, tt), c(0, est[[link]][, k]), col = col, lwd = 1.5)
    matlines(c(0, tt), rbind(0, bands[[link]][[k]]), col = col, lty = 2, lwd = 1)
  }
  lines(c(0, tt), c(0, est$lee[, k]), col = "grey45", lty = 3, lwd = 1.5)
  if (k == 1)
    legend("topleft", bty = "n", cex = .8, col = c("black", "black", "red", "red", "grey45"),
           lty = c(1, 2, 1, 2, 3), lwd = c(1.5, 1, 1.5, 1, 1.5),
           legend = c("Log link, Trend filtering", "95% CI", "Logit link, Trend filtering", "95% CI",
                      "Lee et al. (2018)"))
  mtext(c("(a)", "(b)")[k], side = 1, line = 4.6, cex = .95)
}
invisible(dev.off())
cat("\noutputs:", file.path(OUT_DIR, c("seer_cif.pdf", "seer_cif.csv", "seer_tf_beta_bootse.csv")), "\n")

## ============================================================ exact vs Poisson-approximated form
## Both forms from the same bootstrap replications: point estimate, 95% CI, bootstrap SE at
## 1, 3, 5, 10 years, and the largest difference over months 1-120.
forms <- c("exact", "poisson")
cmp <- do.call(rbind, lapply(forms, function(fm) do.call(rbind, lapply(c("log", "logit"), function(link) {
  F0 <- cif_fit(lapply(res, `[[`, paste0(link, "_tf")), link, fm)
  do.call(rbind, lapply(1:2, function(k) {
    M <- sapply(BT, function(R) cif_fit(R[[link]]$fits, link, fm)[, k])
    data.frame(form = fm, link = link, cause = k, month = 1:J, cif = F0[, k],
               lo = apply(M, 1, quantile, .025), hi = apply(M, 1, quantile, .975), se = apply(M, 1, sd))
  }))
}))))
write.csv(cmp[cmp$month %in% yrs, ], file.path(OUT_DIR, "seer_cif_compare.csv"), row.names = FALSE)
cat("\n== exact vs Poisson form (trend filtering), CIF in % ==\n")
for (link in c("log", "logit")) for (k in 1:2) {
  e <- cmp[cmp$form == "exact" & cmp$link == link & cmp$cause == k, ]
  p <- cmp[cmp$form == "poisson" & cmp$link == link & cmp$cause == k, ]
  cat(sprintf("%-5s cause %d | max |diff| over 120 months: estimate %.4f, lower %.4f, upper %.4f, SE %.5f\n",
              link, k, max(abs(e$cif - p$cif)), max(abs(e$lo - p$lo)), max(abs(e$hi - p$hi)),
              max(abs(e$se - p$se))))
  for (m in yrs) cat(sprintf("        month %3d: exact %5.2f (%5.2f-%5.2f)  poisson %5.2f (%5.2f-%5.2f)\n", m,
                             100 * e$cif[m], 100 * e$lo[m], 100 * e$hi[m], 100 * p$cif[m], 100 * p$lo[m], 100 * p$hi[m]))
}
pdf(file.path(OUT_DIR, "seer_cif_compare.pdf"), width = 9, height = 4.4)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.2, 0.8), oma = c(1.6, 0, 0, 0))
for (k in 1:2) {
  e <- cmp[cmp$form == "exact" & cmp$link == "log" & cmp$cause == k, ]
  p <- cmp[cmp$form == "poisson" & cmp$link == "log" & cmp$cause == k, ]
  plot(NA, xlim = c(0, 120), ylim = c(0, if (k == 1) .6 else .2), xaxt = "n", las = 1,
       xlab = "Time (Months)", ylab = "Cumulative incidence estimates",
       main = c("Death from colon cancer", "Death from other causes")[k])
  axis(1, at = seq(0, 120, 24))
  lines(c(0, e$month), c(0, e$cif), lwd = 1.5); matlines(c(0, e$month), rbind(0, cbind(e$lo, e$hi)), col = 1, lty = 2)
  lines(c(0, p$month), c(0, p$cif), lwd = 1.5, col = "red"); matlines(c(0, p$month), rbind(0, cbind(p$lo, p$hi)), col = "red", lty = 3)
  if (k == 1) legend("topleft", bty = "n", cex = .8, col = c(1, 1, "red", "red"), lty = c(1, 2, 1, 3),
                     legend = c("Exact form", "95% CI", "Poisson-approximated form", "95% CI"))
  mtext(c("(a)", "(b)")[k], side = 1, line = 4.6, cex = .95)
}
invisible(dev.off())
cat("outputs:", file.path(OUT_DIR, c("seer_cif_compare.csv", "seer_cif_compare.pdf")), "\n")
