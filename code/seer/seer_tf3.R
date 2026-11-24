################################################################################
## seer_tf3.R
## SEER stage III colon cancer (1994-2003, age >= 66, surgery), J = 120 months.
## Same estimator and algorithm as the simulation (sim_tf3.R):
##   (alpha, beta) estimated together from the penalized objective, log and logit links,
##   ADMM whose primal step is the BCD sweep repeated to convergence (admm_bcd_exact),
##   b in {1, 2}, rho by BIC along a decreasing grid with warm starts.
##
## Three modes
##   Rscript seer_tf3.R fast              unpenalized fits (log = PL/Breslow, logit at rho = 0),
##                                        Lee et al. (2018), Meir-Gorfine (2025)    [minutes]
##   JOBS=1,5 Rscript seer_tf3.R fit      trend-filtering paths; job ids:
##                                          1-4 log   (k, b) = (1,1) (1,2) (2,1) (2,2)
##                                          5-8 logit (k, b) = (1,1) (1,2) (2,1) (2,2)
##                                        each rho saved as it finishes; a stopped job resumes
##   Rscript seer_tf3.R summary           beta table, baseline figure, path summary
##   (run_seer_tf3.sh runs "fast" and the 8 jobs as 4 independent processes)
##
## Settings (environment variables)
##   RHO_MIN  grid points below this are dropped (default 5; selected rho were >= 34 before)
##   GAMMA_C  "c1,c2": gamma = gamma0 * max(1, c_b * rho), fixed along the run (default "0.3,1";
##            only affects speed; chosen by timing the cause-1 log-link paths on the full data)
##   QUICK=1  smoke test: 1500 subjects, output to seer_out3_quick (a few minutes per job)
## Requires admm_tf_long.R, admm_bcd_exact.R, seer_colon_wide.rds (seer_prep_base.R).
## The data never leave this machine; outputs contain estimates only.
################################################################################
suppressPackageStartupMessages(library(survival))
source("admm_tf_long.R"); source("admm_bcd_exact.R")
if (!file.exists("seer_colon_wide.rds")) stop("seer_colon_wide.rds not found -> run seer_prep_base.R first")

MODE    <- commandArgs(TRUE)[1]; if (is.na(MODE)) MODE <- "summary"
OUT_DIR <- Sys.getenv("SEER_OUT", "seer_out3")
NRHO    <- 15
RHO_MIN <- as.numeric(Sys.getenv("RHO_MIN", "5"))
GAMMA_C <- as.numeric(strsplit(Sys.getenv("GAMMA_C", "0.3,1"), ",")[[1]])
BS      <- 1:2
RUN_MG  <- Sys.getenv("RUN_MG", "1") == "1"
QUICK   <- Sys.getenv("QUICK") == "1"
if (QUICK) OUT_DIR <- paste0(OUT_DIR, "_quick")
dir.create(OUT_DIR, showWarnings = FALSE)
stamp <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")

## ============================================================ 1. data
## Covariates coded as in Lee et al. (2018, Table 3); substage and comorbidity unavailable.
w <- readRDS("seer_colon_wide.rds"); chr <- function(v) as.character(v)
proximal <- c("Cecum", "Ascending Colon", "Hepatic Flexure", "Transverse Colon")
distal   <- c("Splenic Flexure", "Descending Colon", "Sigmoid Colon")
race_lv  <- c("Non-Hispanic White", "Hispanic (All Races)", "Non-Hispanic Black",
              "Non-Hispanic American Indian/Alaska Native", "Non-Hispanic Asian or Pacific Islander")
cv <- data.frame(
  distal  = ifelse(chr(w$site) %in% distal, 1, ifelse(chr(w$site) %in% proximal, 0, NA)),
  grade34 = ifelse(grepl("Grade III|Grade IV", chr(w$grade)), 1,
                   ifelse(grepl("Grade I|Grade II", chr(w$grade)), 0, NA)),
  single  = ifelse(chr(w$marital) == "Unknown", NA,
                   as.numeric(chr(w$marital) != "Married (including common law)")),
  race    = factor(ifelse(chr(w$race) %in% race_lv, chr(w$race), NA), levels = race_lv),
  female  = as.numeric(chr(w$sex) == "Female"),
  age     = pmin(w$age, 90), year = w$year_dx)
keep <- complete.cases(cv); cv <- cv[keep, ]; w <- w[keep, ]
if (QUICK) { set.seed(1); i <- sort(sample(nrow(cv), 1500)); cv <- cv[i, ]; w <- w[i, ] }
X <- cbind(distal = cv$distal, grade34 = cv$grade34, single = cv$single,
           hispanic = as.numeric(cv$race == race_lv[2]), nhblack = as.numeric(cv$race == race_lv[3]),
           nhaian = as.numeric(cv$race == race_lv[4]), nhapi = as.numeric(cv$race == race_lv[5]),
           female = cv$female, age = cv$age - 66, year = cv$year - 1994)
dat <- make_long(w$j, w$event, X)
J <- dat$J; p <- dat$p; n <- dat$n

rho_grid <- sort(exp(seq(log(.05), log(n / 2), length = NRHO)), decreasing = TRUE)
rho_grid <- rho_grid[rho_grid >= RHO_MIN]
JOBS_ALL <- data.frame(id = 1:8, link = rep(c("log", "logit"), each = 4),
                       k = rep(rep(1:2, each = 2), 2), b = rep(BS, 4), stringsAsFactors = FALSE)

## ============================================================ 2. unpenalized fits
## sandwich SE for the partial-likelihood beta (= log-link estimator at rho = 0)
pl_sandwich <- function(dat, k, fit) {
  dk <- if (k == 1L) dat$d1 else dat$d2
  r  <- dat$dc == 1L; Xr <- dat$X[r, , drop = FALSE]; lr <- dat$l[r]
  wr <- exp(as.vector(Xr %*% fit$beta))
  S0 <- gsum(wr, lr, dat$J); S1 <- gsum(Xr * wr, lr, dat$J)
  xbar <- S1 / ifelse(S0 > 0, S0, 1); ev <- dk == 1L
  U_ev <- dat$X[ev, , drop = FALSE] - xbar[dat$l[ev], , drop = FALSE]
  cr   <- wr * (fit$events / ifelse(S0 > 0, S0, 1))[lr]
  U_rs <- -cr * (Xr - xbar[lr, , drop = FALSE])
  U <- rbind(rowsum(U_ev, dat$id[ev]), rowsum(U_rs, dat$id[r]))
  U <- rowsum(U, as.integer(rownames(U)))
  Iinv <- solve(fit$I); V <- Iinv %*% crossprod(U) %*% Iinv
  list(se_model = sqrt(diag(Iinv)), se_robust = sqrt(diag(V)), score_sum = colSums(U))
}

## Newton for (alpha, beta) under a logit link with sandwich SE
##   loss = "pois":  proposed objective at rho = 0;  loss = "binom": Lee et al. (2018)
## Intervals with no cause-k event have alpha = -Inf and are left out of the fit.
newton_logit <- function(dat, k, loss = c("pois", "binom"), maxit = 100, tol = 1e-9) {
  loss <- match.arg(loss)
  dk <- if (k == 1L) dat$d1 else dat$d2
  ek <- tabulate(dat$l[dk == 1L], nbins = dat$J); act <- which(ek > 0); Ja <- length(act)
  rows <- dat$l %in% act; li <- match(dat$l[rows], act); X <- dat$X[rows, , drop = FALSE]
  dkr <- dk[rows]; dcr <- dat$dc[rows]; idr <- dat$id[rows]
  alpha <- qlogis(ek[act] / tabulate(li, nbins = Ja)); beta <- numeric(dat$p)
  terms <- function(eta) {
    mu <- plogis(eta); v <- mu * (1 - mu)
    if (loss == "pois") list(f = sum(dcr * mu) - sum(log(mu[dkr == 1L])),
                             g1 = dcr * v - dkr * (1 - mu), g2 = dcr * v * (1 - 2 * mu) + dkr * v)
    else list(f = -sum(dkr * eta - log1p(exp(eta))), g1 = mu - dkr, g2 = v)
  }
  hess <- function(t) {
    Haa <- diag(gsum(t$g2, li, Ja), Ja); Hab <- gsum(X * t$g2, li, Ja)
    rbind(cbind(Haa, Hab), cbind(t(Hab), crossprod(X * t$g2, X)))
  }
  tt <- terms(alpha[li] + as.vector(X %*% beta))
  for (it in 1:maxit) {
    g <- c(gsum(tt$g1, li, Ja), as.vector(crossprod(X, tt$g1)))
    step <- solve(hess(tt), g); th <- c(alpha, beta); sc <- 1
    repeat {
      thn <- th - sc * step
      tn <- terms(thn[li] + as.vector(X %*% thn[Ja + 1:dat$p]))
      if (is.finite(tn$f) && tn$f <= tt$f + 1e-12) break
      sc <- sc / 2; if (sc < 1e-8) break
    }
    alpha <- thn[1:Ja]; beta <- thn[Ja + 1:dat$p]; tt <- tn
    if (max(abs(sc * step)) < tol) break
  }
  Hinv <- solve(hess(tt))
  Sa <- matrix(0, dat$n, Ja); Sa[cbind(idr, li)] <- tt$g1
  S  <- cbind(Sa, rowsum(X * tt$g1, idr, reorder = TRUE)[as.character(1:dat$n), , drop = FALSE])
  S[is.na(S)] <- 0
  V <- Hinv %*% crossprod(S) %*% Hinv; ib <- Ja + 1:dat$p
  a_full <- rep(-Inf, dat$J); a_full[act] <- alpha
  list(alpha = a_full, beta = beta, f = tt$f, iter = it,
       se_model = sqrt(diag(Hinv)[ib]), se_robust = sqrt(diag(V)[ib]))
}

## Meir & Gorfine (2025): clogit stratified by interval + uniroot baseline
fit_mg <- function(dat, k) {
  dk <- if (k == 1L) dat$d1 else dat$d2
  df <- data.frame(dk = dk, l = dat$l, dat$X)
  fml <- as.formula(paste("dk ~", paste(colnames(dat$X), collapse = " + "), "+ strata(l)"))
  fit <- survival::clogit(fml, data = df, method = "breslow")
  beta <- coef(fit); eta <- as.vector(dat$X %*% beta); alpha <- rep(-Inf, dat$J)
  for (l in 1:dat$J) {
    at <- dat$l == l; N <- sum(dk[at]); if (N == 0 || N >= sum(at)) next
    alpha[l] <- uniroot(function(a) sum(plogis(a + eta[at])) - N, c(-50, 50), tol = 1e-10)$root
  }
  list(alpha = alpha, beta = unname(beta), se_model = unname(sqrt(diag(vcov(fit)))))
}

fast_file <- file.path(OUT_DIR, "fits_fast.rds")
run_fast <- function() {
  if (file.exists(fast_file)) return(readRDS(fast_file))
  res <- list()
  for (k in 1:2) {
    stamp("cause %d: unpenalized fits", k)
    pl <- fit_pl(dat, k); sw <- pl_sandwich(dat, k, pl)
    lg  <- newton_logit(dat, k, "pois")
    lee <- newton_logit(dat, k, "binom")
    mg  <- if (RUN_MG) fit_mg(dat, k) else NULL
    res[[k]] <- list(log = list(beta = pl$beta, se = sw$se_robust, alpha = log(pl$lam0)),
                     logit = list(beta = lg$beta, se = lg$se_robust, alpha = lg$alpha),
                     lee = list(beta = lee$beta, se = lee$se_model, alpha = lee$alpha),
                     mg = if (!is.null(mg)) list(beta = mg$beta, se = mg$se_model, alpha = mg$alpha),
                     events = pl$events)
  }
  saveRDS(res, fast_file); res
}

## ============================================================ 3. trend-filtering paths
## one job = (link, k, b): decreasing rho grid with warm starts; after each rho the path so far
## is saved, so a stopped job resumes at the next rho. Progress: <OUT_DIR>/log_<link>_k?_b?.txt
job_file <- function(j) file.path(OUT_DIR, sprintf("tf_%s_k%d_b%d.rds", j$link, j$k, j$b))
run_job <- function(id) {
  j <- JOBS_ALL[id, ]; out <- job_file(j)
  lg <- file.path(OUT_DIR, sprintf("log_%s_k%d_b%d.txt", j$link, j$k, j$b))
  P <- if (file.exists(out)) readRDS(out) else list(fits = list())
  done <- length(P$fits)
  if (done == length(rho_grid)) return(invisible(P))
  a0 <- if (done) P$fits[[done]]$alpha; b0 <- if (done) P$fits[[done]]$beta
  for (i in (done + 1):length(rho_grid)) {
    t0 <- proc.time()[3]
    f <- admm_bcd_exact(dat, j$k, j$b, rho_grid[i], link = j$link, cvec = GAMMA_C,
                        alpha0 = a0, beta0 = b0)
    a0 <- f$alpha; b0 <- f$beta
    P$fits[[i]] <- c(f[c("alpha", "beta", "z", "gamma", "iter", "n_sweep", "converged",
                         "max_mu", "f", "df", "bic")],
                     rho = rho_grid[i], time = unname(proc.time()[3] - t0))
    saveRDS(P, out)
    cat(sprintf("%s %-5s k%d b%d rho %2d/%d = %8.3g  iter %5d  sweeps %6d  conv %-5s  df %3d  BIC %.2f  (%.0f s)\n",
                format(Sys.time(), "%H:%M:%S"), j$link, j$k, j$b, i, length(rho_grid), rho_grid[i],
                f$iter, f$n_sweep, f$converged, f$df, f$bic, proc.time()[3] - t0),
        file = lg, append = TRUE)
  }
  invisible(P)
}

if (MODE == "fast") { run_fast(); stamp("wrote %s", fast_file); quit(save = "no") }
if (MODE == "fit") {
  ids <- as.integer(strsplit(Sys.getenv("JOBS", "1,2,3,4,5,6,7,8"), ",")[[1]])
  for (id in ids) {
    stamp("job %d: %s link, cause %d, b = %d, %d rho values (>= %g), gamma c = (%s)", id,
          JOBS_ALL$link[id], JOBS_ALL$k[id], JOBS_ALL$b[id], length(rho_grid), RHO_MIN,
          paste(GAMMA_C, collapse = ","))
    run_job(id)
  }
  stamp("done: jobs %s", paste(ids, collapse = ","))
  quit(save = "no")
}

## ============================================================ 4. summary
res <- run_fast()
miss <- !sapply(seq_len(nrow(JOBS_ALL)), function(id) {
  f <- job_file(JOBS_ALL[id, ]); file.exists(f) && length(readRDS(f)$fits) == length(rho_grid) })
if (any(miss)) stop("unfinished trend-filtering jobs: ", paste(which(miss), collapse = ","))
paths <- lapply(seq_len(nrow(JOBS_ALL)), function(id) readRDS(job_file(JOBS_ALL[id, ]))$fits)

## per (link, k): rho by BIC within each b, then b by BIC
path_tab <- do.call(rbind, lapply(seq_along(paths), function(id) {
  j <- JOBS_ALL[id, ]
  do.call(rbind, lapply(paths[[id]], function(f)
    data.frame(link = j$link, cause = j$k, b = j$b, rho = f$rho, df = f$df, bic = f$bic,
               iter = f$iter, sweeps = f$n_sweep, converged = f$converged,
               max_mu = f$max_mu, gamma = f$gamma, time_s = f$time)))
}))
write.csv(path_tab, file.path(OUT_DIR, "seer_tf_paths.csv"), row.names = FALSE)
for (link in c("log", "logit")) for (k in 1:2) {
  ids <- which(JOBS_ALL$link == link & JOBS_ALL$k == k)
  best <- lapply(ids, function(id) { b <- sapply(paths[[id]], `[[`, "bic")
                                     paths[[id]][[which.min(b)]] })
  bb <- sapply(best, `[[`, "bic"); i <- which.min(bb); f <- best[[i]]
  res[[k]][[paste0(link, "_tf")]] <- list(alpha = f$alpha, beta = f$beta, b = JOBS_ALL$b[ids[i]],
    rho = f$rho, df = f$df, converged = f$converged, bic_by_b = setNames(bb, paste0("b", BS)),
    edge = f$rho %in% range(rho_grid))
  stamp("cause %d %-5s TF: b = %d, rho = %.3g, df = %d, converged = %s, BIC by b = %s%s",
        k, link, JOBS_ALL$b[ids[i]], f$rho, f$df, f$converged,
        paste(sprintf("%.1f", bb), collapse = " / "),
        if (f$rho %in% range(rho_grid)) "  <- at the edge of the grid" else "")
}
stamp("all paths: %d of %d fits converged; logit max hazard %.3f (must be < 0.5)",
      sum(path_tab$converged), nrow(path_tab), max(path_tab$max_mu[path_tab$link == "logit"]))
saveRDS(list(res = res, X_names = colnames(X), J = J, n = n, rho_grid = rho_grid),
        file.path(OUT_DIR, "seer_fits.rds"))

## ---- beta table (TF beta reported without SE)
row <- function(method, k, R) data.frame(cause = k, method = method, covariate = colnames(X),
  estimate = R$beta, se = if (is.null(R$se)) NA else R$se,
  p = if (is.null(R$se)) NA else 2 * pnorm(-abs(R$beta / R$se)))
tab <- do.call(rbind, lapply(1:2, function(k) { R <- res[[k]]
  rbind(row("Log, unpenalized (PL)", k, R$log),  row("Log, TF", k, R$log_tf),
        row("Logit, unpenalized", k, R$logit),   row("Logit, TF", k, R$logit_tf),
        row("Lee et al. (2018)", k, R$lee),
        if (!is.null(R$mg)) row("Meir and Gorfine (2025)", k, R$mg)) }))
write.csv(tab, file.path(OUT_DIR, "seer_beta.csv"), row.names = FALSE)
for (k in 1:2) {
  T <- tab[tab$cause == k, ]
  f <- ifelse(is.na(T$se), sprintf("%.3f", T$estimate), sprintf("%.3f (%.3f)", T$estimate, T$se))
  cat(sprintf("\n== cause %d: estimate (robust SE) ==\n", k))
  print(noquote(tapply(f, list(factor(T$covariate, colnames(X)), factor(T$method, unique(T$method))),
                       identity)))
}

## ---- baseline figure: unpenalized (points) and TF (line), log and logit links
pdf(file.path(OUT_DIR, "seer_baseline.pdf"), width = 9, height = 7.5)
par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 0.8))
fin <- function(a) { a[!is.finite(a)] <- NA; a }
for (k in 1:2) {
  R <- res[[k]]; lab <- c("colon cancer", "other causes")[k]
  pos <- c("bottomleft", "bottomright")[k]; zero <- which(R$events == 0)
  yl <- range(fin(c(R$log$alpha, R$log_tf$alpha)), na.rm = TRUE)
  plot(1:J, fin(R$log$alpha), pch = 16, cex = .45, col = "grey55", ylim = yl, las = 1,
       xlab = "Months since diagnosis", ylab = "Baseline log-hazard",
       main = sprintf("Death from %s: log link", lab))
  lines(1:J, R$log_tf$alpha, col = "#2166ac", lwd = 2)
  if (length(zero)) rug(zero, col = "grey40")
  legend(pos, bty = "n", cex = .75, pch = c(16, NA), lwd = c(NA, 2), col = c("grey55", "#2166ac"),
         legend = c("Unpenalized", sprintf("Trend filtering (b = %d)", R$log_tf$b)))
  yl <- range(fin(c(R$logit$alpha, R$logit_tf$alpha, R$lee$alpha)), na.rm = TRUE)
  plot(1:J, fin(R$logit$alpha), pch = 16, cex = .45, col = "grey55", ylim = yl, las = 1,
       xlab = "Months since diagnosis", ylab = "Baseline log-odds",
       main = sprintf("Death from %s: logit link", lab))
  lines(1:J, fin(R$lee$alpha), col = "#d7301f", lwd = 1)
  lines(1:J, R$logit_tf$alpha, col = "#2166ac", lwd = 2)
  if (length(zero)) rug(zero, col = "grey40")
  legend(pos, bty = "n", cex = .75, pch = c(16, NA, NA), lwd = c(NA, 1, 2),
         col = c("grey55", "#d7301f", "#2166ac"),
         legend = c("Unpenalized", "Lee et al. (2018)", sprintf("Trend filtering (b = %d)", R$logit_tf$b)))
}
invisible(dev.off())
cat("\noutputs:", file.path(OUT_DIR, c("seer_fits.rds", "seer_beta.csv", "seer_tf_paths.csv",
                                       "seer_baseline.pdf")), "\n")
