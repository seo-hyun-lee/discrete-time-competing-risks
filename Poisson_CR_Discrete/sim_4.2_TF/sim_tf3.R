################################################################################
## sim_tf3.R -- trend-filtering simulation (Section 4.2), revised
##   Both DGPs: (alpha, beta) estimated together from the penalized Poisson-approximated
##   objective, unpenalized vs trend filtering; beta is never fixed in advance.
##     DGP "ph"    (grouped PH)  -> log link
##     DGP "logit" (logit)       -> logit link
##   Estimators per cause:
##     unpen : rho = 0 (under the log link its beta equals the partial-likelihood beta;
##             intervals without cause-k events get alpha = -Inf)
##     lee   : Lee et al. (2018) collapsed binomial GLM (comparator)
##     tf    : ADMM, BCD primal step solved to convergence, fixed gamma (admm_bcd_exact.R),
##             b = 1, 2; rho by BIC within b; BIC across b stored for selecting b
##   Env: SHAPE (smooth | step | abrupt), N, J, B, SEED, MC_CORES (keep 1 on Windows), OUT,
##        DGPS ("ph,logit"), NRHO (rho grid length, default 12),
##        GAMMA_C ("c1,c2": gamma = gamma0 * max(1, c_b * rho); default "1,1", fastest at J = 30),
##        REPS (e.g. "1:25": run only these replications; for splitting one shape over several
##              independent Rscript processes. The final .rds is written by whichever process
##              finds all B replications saved),
##        RHO_MIN (drop grid points below this; default 0 = full grid. The smallest rho are
##                 the slowest to fit and are rarely chosen by BIC; check pct_at_min in the summary)
##   Needs admm_tf_long.R (DGP, make_long, fit_lee, helpers) and admm_bcd_exact.R.
##   Same DGP, true baselines, rho grid and seed handling as sim_tf2.R.
################################################################################
source("admm_tf_long.R")
source("admm_bcd_exact.R")

B      <- as.integer(Sys.getenv("B", "100"))
n      <- as.integer(Sys.getenv("N", "1000"))
J      <- as.integer(Sys.getenv("J", "30"))
SHAPE  <- Sys.getenv("SHAPE", "smooth")
SEED   <- as.integer(Sys.getenv("SEED", "1"))
CORES  <- as.integer(Sys.getenv("MC_CORES", "1"))
NRHO   <- as.integer(Sys.getenv("NRHO", "12"))
DGPS   <- strsplit(Sys.getenv("DGPS", "ph,logit"), ",")[[1]]
OUT    <- Sys.getenv("OUT", sprintf("sim_tf3_%s_n%d_J%d.rds", SHAPE, n, J))
BS     <- 1:2
b1 <- c(0.5, -0.3); b2 <- c(-0.2, 0.4)
GAMMA_C <- as.numeric(strsplit(Sys.getenv("GAMMA_C", "1,1"), ",")[[1]])
RHO_MIN <- as.numeric(Sys.getenv("RHO_MIN", "0"))
rg <- sort(exp(seq(log(.05), log(n / 2), length = NRHO)), decreasing = TRUE)
rg <- rg[rg >= RHO_MIN]

## ----------------------------------------------------------- true baselines
make_shape <- function(shape, J) {
  t <- 1:J; s <- J / 30                       # shapes defined on J = 30, rescaled
  tt <- ceiling(t / s)
  if (shape == "smooth") {
    h1 <- 0.015 + 0.065 * exp(-((tt - 12)/6)^2)
    h2 <- 0.008 + 0.042 * tt / 30
  } else if (shape == "step") {
    h1 <- ifelse(tt <= 10, .03, ifelse(tt <= 20, .07, .03))
    h2 <- ifelse(tt <= 15, .02, .04)
  } else if (shape == "abrupt") {
    h1 <- rep(.02, J); h1[tt %in% 9:12] <- c(.05, .09, .09, .05)[match(tt[tt %in% 9:12], 9:12)]
    h2 <- 0.01 + 0.07 * exp(-(tt - 1)/3)
  } else stop("unknown SHAPE")
  list(h1 = h1, h2 = h2)
}
sh <- make_shape(SHAPE, J); h1 <- sh$h1; h2 <- sh$h2
alp0 <- list(ph    = list(true_alpha(h1, h2, 1, "ph"),    true_alpha(h1, h2, 2, "ph")),
             logit = list(true_alpha(h1, h2, 1, "logit"), true_alpha(h1, h2, 2, "logit")))

## ------------------------------------------------ unpenalized fit (rho = 0)
## Newton on (alpha, beta) jointly. Intervals without cause-k events have alpha_l = -Inf at the
## minimum; their records then contribute nothing, so they are dropped and alpha_l set to -Inf.
fit_unpen <- function(dat, k, link, maxit = 100, tol = 1e-9) {
  J <- dat$J; p <- dat$p
  dk <- if (k == 1L) dat$d1 else dat$d2
  dl <- gsum(dk, dat$l, J); keep <- which(dl > 0)
  r  <- dat$l %in% keep
  l  <- match(dat$l[r], keep); X <- dat$X[r, , drop = FALSE]; dc <- dat$dc[r]; dk <- dk[r]
  Jk <- length(keep)
  q  <- pmin(dl[keep] / tabulate(l, nbins = Jk), .45)
  th <- c(if (link == "logit") qlogis(q) else log(q), numeric(p))
  obj <- function(th) {
    eta <- th[l] + as.vector(X %*% th[Jk + 1:p])
    if (link == "logit" && max(eta) > 0) return(Inf)
    obj_f(eta, dc, dk, link)
  }
  F0 <- obj(th)
  for (it in 1:maxit) {
    eta <- th[l] + as.vector(X %*% th[Jk + 1:p])
    tt  <- link_terms(eta, dc, dk, link)
    g   <- c(gsum(tt$g1, l, Jk), as.vector(crossprod(X, tt$g1)))
    Hab <- gsum(X * tt$g2, l, Jk)
    H   <- rbind(cbind(diag(gsum(tt$g2, l, Jk) + 1e-12, Jk), Hab),
                 cbind(t(Hab), crossprod(X * tt$g2, X)))
    st  <- solve(H, g); t <- 1
    repeat { F1 <- obj(th - t * st); if (F1 <= F0 - 1e-4 * t * sum(g * st) || t < 1e-10) break; t <- t / 2 }
    th <- th - t * st; done <- abs(F0 - F1) < tol * (1 + abs(F0)); F0 <- F1
    if (done || max(abs(g)) < 1e-8) break
  }
  alpha <- rep(-Inf, J); alpha[keep] <- th[1:Jk]
  list(alpha = alpha, beta = th[Jk + 1:p], converged = it < maxit, n_zero = J - Jk)
}

## ------------------------------------------------------------------ one rep
one_rep <- function(rep) {
  out <- list()
  for (dg in DGPS) {
    link <- if (dg == "ph") "log" else "logit"
    d <- gen_data(n, J, h1, h2, b1, b2, dg); dat <- make_long(d$Y, d$D, d$X)
    R <- list(unpen = matrix(NA, J, 2), lee = matrix(NA, J, 2),
              tf = array(NA_real_, c(J, 2, length(BS))), rho = matrix(NA, 2, length(BS)),
              df = matrix(NA, 2, length(BS)), bic = matrix(NA, 2, length(BS)),
              nonconv = matrix(NA, 2, length(BS)), iter = matrix(NA, 2, length(BS)),
              max_mu = matrix(NA, 2, length(BS)), n_zero = rep(NA, 2),
              beta_u = rep(NA, 4), beta_lee = rep(NA, 4), beta_t = matrix(NA, 4, length(BS)))
    if (dat$J < J) { out[[dg]] <- R; next }
    for (k in 1:2) {
      ix <- (k - 1) * 2 + 1:2
      u0 <- fit_unpen(dat, k, link)
      R$unpen[, k] <- u0$alpha; R$beta_u[ix] <- u0$beta; R$n_zero[k] <- u0$n_zero
      le <- fit_lee(dat, k); R$lee[, k] <- le$alpha; R$beta_lee[ix] <- le$beta
      for (bo in BS) {
        pa <- tf_path_bcd(dat, k, bo, rg, link = link, cvec = GAMMA_C)
        j  <- pa$best; fj <- pa$fits[[j]]
        R$tf[, k, bo] <- fj$alpha; R$beta_t[ix, bo] <- fj$beta
        R$rho[k, bo] <- pa$rho; R$df[k, bo] <- fj$df; R$bic[k, bo] <- min(pa$bic)
        R$nonconv[k, bo] <- pa$n_nonconv; R$iter[k, bo] <- fj$iter; R$max_mu[k, bo] <- fj$max_mu
      }
    }
    out[[dg]] <- R
  }
  out
}

## Each replication has its own seed (SEED, rep), so results do not depend on MC_CORES, and is
## saved to sim_tf3_parts/<shape>_n<n>_J<J>/rep_###.rds as soon as it finishes: a stopped run
## resumes where it left off. Progress: one line per finished replication in
## sim_tf3_progress_<shape>_n<n>_J<J>.log (time, rep, seconds, number done).
PART <- file.path("sim_tf3_parts", sprintf("%s_n%d_J%d", SHAPE, n, J))
PROG <- sprintf("sim_tf3_progress_%s_n%d_J%d.log", SHAPE, n, J)
dir.create(PART, recursive = TRUE, showWarnings = FALSE)
part_file <- function(r) file.path(PART, sprintf("rep_%03d.rds", r))
run_rep <- function(r) {
  f <- part_file(r)
  if (file.exists(f)) return(invisible(NULL))
  RNGkind("L'Ecuyer-CMRG"); set.seed(SEED * 100000 + r)
  t1 <- proc.time()[3]
  res <- one_rep(r)
  saveRDS(res, f)
  cat(sprintf("%s  rep %3d done  %5.0f s  (%d / %d finished)\n", format(Sys.time(), "%m-%d %H:%M"), r,
              proc.time()[3] - t1, length(list.files(PART, "^rep_")), B), file = PROG, append = TRUE)
  invisible(NULL)
}
t0 <- proc.time()
REPS <- eval(parse(text = Sys.getenv("REPS", "1:B")))
todo <- REPS[!file.exists(part_file(REPS))]
cat(sprintf("%s  start %s: %d of %d replications to run on %d core(s)\n",
            format(Sys.time(), "%m-%d %H:%M"), SHAPE, length(todo), B, CORES), file = PROG, append = TRUE)
if (CORES > 1) {
  out <- parallel::mclapply(todo, function(r) try(run_rep(r)), mc.cores = CORES, mc.preschedule = FALSE)
  bad <- sapply(out, inherits, "try-error"); if (any(bad)) stop(out[bad][[1]])
} else for (r in todo) run_rep(r)
if (!all(file.exists(part_file(1:B)))) {
  cat(sprintf("reps %s finished; %d of %d saved so far, final .rds not written yet\n",
              Sys.getenv("REPS", "1:B"), sum(file.exists(part_file(1:B))), B))
  quit(save = "no")
}
reps <- lapply(1:B, function(r) readRDS(part_file(r)))

## ------------------------------------------------------------------ collect
pick <- function(dg, nm) lapply(reps, function(r) r[[dg]][[nm]])
arr  <- function(L) { d <- dim(L[[1]]); if (is.null(d)) d <- length(L[[1]])
                      a <- array(unlist(L), c(d, length(L)))
                      aperm(a, c(length(d) + 1, seq_along(d))) }
nms <- c("unpen", "lee", "tf", "rho", "df", "bic", "nonconv", "iter", "max_mu", "n_zero",
         "beta_u", "beta_lee", "beta_t")
res <- list()
for (dg in DGPS) res[[dg]] <- lapply(setNames(nm = nms), function(nm) arr(pick(dg, nm)))
saveRDS(list(res = res, alp0 = alp0, h1 = h1, h2 = h2, n = n, J = J, B = B, shape = SHAPE,
             b1 = b1, b2 = b2, BS = BS, rho_grid = rg, gamma_c = GAMMA_C), OUT)
cat("done", OUT, "in", round((proc.time() - t0)[3]), "s\n")
