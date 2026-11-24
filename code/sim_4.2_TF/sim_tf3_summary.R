################################################################################
## sim_tf3_summary.R -- tables for the revised trend-filtering simulation (sim_tf3.R output)
##   Rscript sim_tf3_summary.R sim_tf3_smooth_n1000_J30.rds [sim_tf3_step_...rds ...]
## Writes sim_tf3_beta.csv, sim_tf3_baseline.csv, sim_tf3_conv.csv and prints them.
##   beta     : bias, SD, RMSE of each coefficient (truth = b1 for cause 1, b2 for cause 2)
##   baseline : on the hazard scale lambda_0 = g^{-1}(alpha) (x = 0)
##              RMSE = sqrt(mean over l and reps of (lambda_hat - lambda_0)^2)
##              TV   = mean over reps of sum_l |lambda_hat_l - lambda_hat_{l-1}| (roughness);
##              TV_true is the same quantity for the true baseline
##   estimators: unpen (rho = 0), lee (Lee et al. 2018), tf_b1, tf_b2, tf_sel (b chosen by BIC)
################################################################################
files <- commandArgs(TRUE)
if (!length(files)) files <- Sys.glob("sim_tf3_*_n*_J*.rds")

inv  <- function(dg) if (dg == "ph") exp else plogis
tv   <- function(lam) sum(abs(diff(lam)))
BT <- BL <- BC <- list()
for (fn in files) {
  S <- readRDS(fn)
  for (dg in names(S$res)) {
    R <- S$res[[dg]]; g <- inv(dg)
    ok <- apply(R$unpen, 1, function(a) !all(is.na(a)))          # reps with complete data
    nb <- length(S$BS)
    for (k in 1:2) {
      ix <- (k - 1) * 2 + 1:2; tb <- if (k == 1) S$b1 else S$b2
      sel <- apply(R$bic[ok, k, , drop = FALSE], 1, which.min)   # b chosen by BIC, per rep
      est_beta <- list(unpen = R$beta_u[ok, ix, drop = FALSE],
                       lee   = R$beta_lee[ok, ix, drop = FALSE])
      M <- function(a) matrix(a, nrow = sum(ok))
      est_lam  <- list(unpen = g(M(R$unpen[ok, , k])), lee = plogis(M(R$lee[ok, , k])))
      for (bo in seq_len(nb)) {
        est_beta[[paste0("tf_b", S$BS[bo])]] <- M(R$beta_t[ok, ix, bo])
        est_lam[[paste0("tf_b", S$BS[bo])]]  <- g(M(R$tf[ok, , k, bo]))
      }
      est_beta$tf_sel <- M(t(sapply(seq_along(sel), function(r) R$beta_t[which(ok)[r], ix, sel[r]])))
      est_lam$tf_sel  <- M(t(sapply(seq_along(sel), function(r) g(R$tf[which(ok)[r], , k, sel[r]]))))
      lam0 <- g(S$alp0[[dg]][[k]])
      for (e in names(est_beta)) {
        bt <- est_beta[[e]]
        for (j in 1:2) BT[[length(BT) + 1]] <- data.frame(
          shape = S$shape, dgp = dg, cause = k, coef = paste0("beta", j), estimator = e,
          true = tb[j], bias = mean(bt[, j]) - tb[j], sd = sd(bt[, j]),
          rmse = sqrt(mean((bt[, j] - tb[j])^2)))
        lh <- est_lam[[e]]
        BL[[length(BL) + 1]] <- data.frame(
          shape = S$shape, dgp = dg, cause = k, estimator = e,
          rmse = sqrt(mean(sweep(lh, 2, lam0)^2, na.rm = TRUE)), tv = mean(apply(lh, 1, tv)),
          tv_true = tv(lam0))
      }
      for (bo in seq_len(nb)) BC[[length(BC) + 1]] <- data.frame(
        shape = S$shape, dgp = dg, cause = k, b = S$BS[bo],
        reps = sum(ok), pct_selected = mean(sel == bo) * 100,
        median_rho = median(R$rho[ok, k, bo]),
        pct_at_min = mean(R$rho[ok, k, bo] == min(S$rho_grid)) * 100,    # rho at grid boundary:
        pct_at_max = mean(R$rho[ok, k, bo] == max(S$rho_grid)) * 100,    #   extend the grid if large
        median_df = median(R$df[ok, k, bo]),
        nonconv_fits = sum(R$nonconv[ok, k, bo]), max_iter = max(R$iter[ok, k, bo]),
        max_mu = max(R$max_mu[ok, k, bo]))
    }
  }
}
fmt <- function(d) { num <- sapply(d, is.numeric); d[num] <- lapply(d[num], signif, 4); d }
beta <- fmt(do.call(rbind, BT)); base <- fmt(do.call(rbind, BL)); conv <- fmt(do.call(rbind, BC))
write.csv(beta, "sim_tf3_beta.csv", row.names = FALSE)
write.csv(base, "sim_tf3_baseline.csv", row.names = FALSE)
write.csv(conv, "sim_tf3_conv.csv", row.names = FALSE)
cat("\n== beta ==\n");     print(beta, row.names = FALSE)
cat("\n== baseline (hazard scale) ==\n"); print(base, row.names = FALSE)
cat("\n== convergence / b selection ==\n"); print(conv, row.names = FALSE)
