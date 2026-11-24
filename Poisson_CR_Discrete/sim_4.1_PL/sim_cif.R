################################################################################
## sim_cif.R -- simulation for cumulative incidence prediction (Section 4.3)
##   SETTING=s1n300 Rscript sim_cif.R      (one setting per run; see CFG below)
##   Rscript sim_cif.R summary              (collects the finished settings -> sim_cif_summary.csv)
##
## Data: the two mechanisms of Section 4.1 (gen_data in admm_tf_long.R; "ph" = grouped
##   continuous-time PH, "logit" = discrete logit hazards), X = (N(0,1), Bernoulli(0.5)),
##   beta1 = (0.5, -0.3), beta2 = (-0.2, 0.4), administrative censoring at J.
##   Baseline hazards increase linearly, h1: 0.05 -> 0.15, h2: 0.03 -> 0.12 (scale 1), or one
##   fifth of these (scale 0.2).
##     s1n300 : scale 1,   J = 10, n = 300      s1n600 : scale 1,   J = 10, n = 600
##     s02J10 : scale 0.2, J = 10, n = 2000     s02J30 : scale 0.2, J = 30, n = 2000
## Predictors of F_k(t_m; x), x = (0,0) and (1,1), m = J/2 and J, all from the product form
##   F_k(t_m; x) = sum_{l<=m} lam_k(t_l; x) prod_{r<l} {1 - lam_1(t_r; x) - lam_2(t_r; x)}:
##   Proposed : beta by partial likelihood, Breslow-type baseline (fit_pl); SE by the
##              sandwich-based delta method of Section 2.3 (both causes stacked, so the
##              between-cause covariance is included)
##   Lee      : collapsed binomial likelihood (fit_lee, logit link); SE by the delta method
##              with the cluster-robust sandwich covariance of the stacked logistic scores
## Truth: product form at the true discrete hazards of the generating mechanism.
## Env: SETTING, B (default 1000), SEED (default 2026), MC_CORES (default 4; 1 on Windows).
## Needs admm_tf_long.R. Base R + parallel.
################################################################################
source("admm_tf_long.R")
MODE <- commandArgs(TRUE)[1]; if (is.na(MODE)) MODE <- "run"
CFG <- list(s1n300 = list(scale = 1,   J = 10, n = 300),
            s1n600 = list(scale = 1,   J = 10, n = 600),
            s02J10 = list(scale = 0.2, J = 10, n = 2000),
            s02J30 = list(scale = 0.2, J = 30, n = 2000))
b1 <- c(.5, -.3); b2 <- c(-.2, .4)
PROF <- list(x00 = c(0, 0), x11 = c(1, 1))

## ------------------------------------------------------------------ CIF and gradient
## eta: 2 x J linear predictors at x; link "log" or "logit"; returns, for m in ms,
## F[k] and dF_k/deta (2 x J; zero beyond m)
cif_grad <- function(eta, link, ms) {
  lam <- if (link == "log") exp(eta) else plogis(eta)
  dl  <- if (link == "log") lam else lam * (1 - lam)
  lam[!is.finite(eta)] <- 0; dl[!is.finite(eta)] <- 0
  J <- ncol(lam); L <- colSums(lam)
  Sprev <- c(1, cumprod(1 - L))[1:J]                     # S(t_{l-1})
  out <- list()
  for (m in ms) for (k in 1:2) {
    f <- lam[k, 1:m] * Sprev[1:m]                        # sub-densities
    tail <- rev(cumsum(rev(f)))                          # sum_{l >= q} f_l
    after <- c(tail[-1], 0)                              # sum_{l > q} f_l
    G <- matrix(0, 2, J)
    for (kk in 1:2) G[kk, 1:m] <- dl[kk, 1:m] * ((kk == k) * Sprev[1:m] - after / (1 - L[1:m]))
    out[[paste(m, k)]] <- list(F = sum(f), G = G)
  }
  out
}

## ---------------------------------------------- stacked scores and information
## For each cause: residual r = dk - w * mu (w = dc for the proposed working likelihood,
## w = 1 for the binomial), derivative weight for A; columns = active intervals, then x.
score_info <- function(dat, alpha, beta, k, type) {
  dk  <- if (k == 1L) dat$d1 else dat$d2
  act <- which(is.finite(alpha)); Ja <- length(act)
  eta <- alpha[dat$l] + as.vector(dat$X %*% beta)
  if (type == "prop") { mu <- ifelse(is.finite(eta), exp(eta), 0); r <- dk - dat$dc * mu; h <- dat$dc * mu
  } else              { mu <- plogis(eta); r <- dk - mu; h <- mu * (1 - mu) }
  Pa <- matrix(0, dat$n, dat$J); Pa[cbind(dat$id, dat$l)] <- r
  Pb <- matrix(0, dat$n, dat$p); rb <- rowsum(r * dat$X, dat$id, reorder = TRUE)
  Pb[as.integer(rownames(rb)), ] <- rb
  Haa <- gsum(h, dat$l, dat$J)[act]; Hab <- gsum(dat$X * h, dat$l, dat$J)[act, , drop = FALSE]
  A <- rbind(cbind(diag(Haa, Ja), Hab), cbind(t(Hab), crossprod(dat$X * h, dat$X)))
  list(Psi = cbind(Pa[, act, drop = FALSE], Pb), A = A, act = act)
}
## delta-method SE for every (m, k) at profile x, from the two causes' fits
delta_se <- function(S1, S2, alpha, beta, x, link, ms) {
  Ainv <- matrix(0, ncol(S1$Psi) + ncol(S2$Psi), ncol(S1$Psi) + ncol(S2$Psi))
  i1 <- seq_len(ncol(S1$Psi)); i2 <- ncol(S1$Psi) + seq_len(ncol(S2$Psi))
  Ainv[i1, i1] <- solve(S1$A); Ainv[i2, i2] <- solve(S2$A)
  Psi <- cbind(S1$Psi, S2$Psi)
  V <- Ainv %*% crossprod(Psi) %*% Ainv
  eta <- rbind(alpha[[1]] + sum(x * beta[[1]]), alpha[[2]] + sum(x * beta[[2]]))
  cg <- cif_grad(eta, link, ms)
  sapply(cg, function(o) {
    g <- c(o$G[1, S1$act], x * sum(o$G[1, ]), o$G[2, S2$act], x * sum(o$G[2, ]))
    c(F = o$F, se = sqrt(drop(t(g) %*% V %*% g)))
  })
}

true_cif <- function(h1, h2, x, dgp, ms) {
  if (dgp == "ph") {
    a1 <- h1 * exp(sum(x * b1)); a2 <- h2 * exp(sum(x * b2)); pev <- 1 - exp(-(a1 + a2))
    lam <- rbind(pev * a1 / (a1 + a2), pev * a2 / (a1 + a2))
  } else lam <- rbind(plogis(qlogis(h1) + sum(x * b1)), plogis(qlogis(h2) + sum(x * b2)))
  sapply(cif_grad(log(lam), "log", ms), `[[`, "F")
}

one_rep <- function(r, cf, h1, h2, ms) {
  set.seed(SEED * 100000 + r)
  out <- list()
  for (dgp in c("ph", "logit")) {
    d <- gen_data(cf$n, cf$J, h1, h2, b1, b2, dgp); dat <- make_long(d$Y, d$D, d$X)
    if (dat$J < cf$J) return(NULL)                     # no one reached interval J (does not occur here)
    pl <- lapply(1:2, function(k) fit_pl(dat, k))
    ap <- lapply(pl, function(f) log(f$lam0)); bp <- lapply(pl, `[[`, "beta")
    Sp <- lapply(1:2, function(k) score_info(dat, ap[[k]], bp[[k]], k, "prop"))
    le <- lapply(1:2, function(k) fit_lee(dat, k))
    al <- lapply(le, function(f) unname(f$alpha)); bl <- lapply(le, function(f) unname(f$beta))
    Sl <- lapply(1:2, function(k) score_info(dat, al[[k]], bl[[k]], k, "lee"))
    for (pn in names(PROF)) {
      x <- PROF[[pn]]
      out[[paste(dgp, pn, "Proposed")]] <- delta_se(Sp[[1]], Sp[[2]], ap, bp, x, "log", ms)
      out[[paste(dgp, pn, "Lee")]]      <- delta_se(Sl[[1]], Sl[[2]], al, bl, x, "logit", ms)
    }
  }
  out
}

SEED <- as.integer(Sys.getenv("SEED", "2026"))
if (MODE == "run") {
  SETTING <- Sys.getenv("SETTING", "s1n300"); cf <- CFG[[SETTING]]
  if (is.null(cf)) stop("SETTING must be one of ", paste(names(CFG), collapse = ", "))
  B <- as.integer(Sys.getenv("B", "1000")); CORES <- as.integer(Sys.getenv("MC_CORES", "4"))
  h1 <- cf$scale * seq(.05, .15, length = cf$J); h2 <- cf$scale * seq(.03, .12, length = cf$J)
  ms <- c(cf$J / 2, cf$J)
  t0 <- proc.time()[3]
  run <- function(r) tryCatch(one_rep(r, cf, h1, h2, ms), error = function(e) NULL)
  R <- if (CORES > 1) parallel::mclapply(1:B, run, mc.cores = CORES) else lapply(1:B, run)
  truth <- list()
  for (dgp in c("ph", "logit")) for (pn in names(PROF)) truth[[paste(dgp, pn)]] <- true_cif(h1, h2, PROF[[pn]], dgp, ms)
  saveRDS(list(R = R, truth = truth, cf = cf, ms = ms, B = B, seed = SEED, h1 = h1, h2 = h2),
          sprintf("sim_cif_%s.rds", SETTING))
  cat(sprintf("%s: %d of %d replications usable, %.1f min -> sim_cif_%s.rds\n", SETTING,
              sum(!sapply(R, is.null)), B, (proc.time()[3] - t0) / 60, SETTING))
  quit(save = "no")
}

## ------------------------------------------------------------------ summary
rows <- list()
for (s in names(CFG)) {
  f <- sprintf("sim_cif_%s.rds", s); if (!file.exists(f)) next
  S <- readRDS(f); R <- Filter(Negate(is.null), S$R)
  for (key in names(R[[1]])) {
    parts <- strsplit(key, " ")[[1]]; tr <- S$truth[[paste(parts[1], parts[2])]]
    M <- simplify2array(lapply(R, `[[`, key))           # 2 x (m,k) x reps
    for (j in seq_len(ncol(M))) {
      mk <- as.integer(strsplit(colnames(M)[j], " ")[[1]])
      est <- M["F", j, ]; se <- M["se", j, ]
      rows[[length(rows) + 1]] <- data.frame(setting = s, n = S$cf$n, J = S$cf$J, scale = S$cf$scale,
        dgp = parts[1], profile = parts[2], method = parts[3], cause = mk[2], m = mk[1],
        truth = tr[j], bias = mean(est) - tr[j], rel_bias = 100 * (mean(est) - tr[j]) / tr[j],
        ASE = mean(se), ESE = sd(est), CP = mean(abs(est - tr[j]) <= 1.96 * se), B = length(est))
    }
  }
}
if (!length(rows)) stop("no sim_cif_<setting>.rds found")
T <- do.call(rbind, rows)
write.csv(T, "sim_cif_summary.csv", row.names = FALSE)
cat("== Bias, ASE, ESE, CP (x = (0,0)); full table in sim_cif_summary.csv ==\n")
P <- T[T$profile == "x00", ]
P[, c("truth", "bias", "ASE", "ESE")] <- round(P[, c("truth", "bias", "ASE", "ESE")], 4)
P$rel_bias <- round(P$rel_bias, 1); P$CP <- round(P$CP, 3)
print(P[order(P$setting, P$dgp, P$cause, P$m, P$method),
        c("setting", "dgp", "cause", "m", "method", "truth", "bias", "rel_bias", "ASE", "ESE", "CP")],
      row.names = FALSE)
