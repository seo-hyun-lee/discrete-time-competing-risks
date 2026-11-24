################################################################################
## admm_bcd_exact.R
## Trend-filtered estimation of (alpha, beta), log or logit link, as in the original algorithm
## (ADMM with a block-coordinate primal step), except that the primal step is solved EXACTLY:
##   within each ADMM iteration, the BCD sweep of the original algorithm (one damped Newton step
##   on alpha, then one on beta) is repeated until the full gradient of the primal subproblem
##   satisfies ||g||_inf < newton_tol.
##   (alpha, beta) = argmin  f(alpha, beta) + rho * ||D_b alpha||_1,
##   f(alpha, beta) = sum_{person-period records} { dc * mu(eta) - dk * log mu(eta) },
##   eta = alpha_l + x' beta,  mu = exp(eta) (log) or expit(eta) (logit).
## logit: minimised over Theta = { eta <= 0 for every record } (hazard <= 1/2), on which f is
##        convex; the line search rejects points outside Theta. Returns max_mu to check that the
##        constraint is inactive (max_mu < 1/2).
## gamma FIXED for the run: gamma0 * max(1, c_b * rho), b in {1, 2}, c = (0, 1).
## Stopping rule of Boyd et al. (2011, Sec. 3.3.1) with unscaled u.
## Needs admm_tf_long.R (gsum, getD, link_terms, obj_f).
################################################################################

admm_bcd_exact <- function(dat, k, b, rho, link = c("log", "logit"), cvec = c(0, 1),
                           gamma = NULL, alpha0 = NULL, beta0 = NULL, u0 = NULL, maxit = 20000,
                           tol_abs = 1e-8, tol_rel = 1e-5, newton_tol = 1e-6,
                           max_sweep = 200, center = TRUE, trace = 0) {
  link <- match.arg(link)
  J <- dat$J; p <- dat$p; l <- dat$l; dc <- dat$dc; id <- dat$id
  Xs <- dat$X[match(seq_len(dat$n), id), , drop = FALSE]   # subject-level covariates
  ## centring x is a reparametrisation alpha_l -> alpha_l + xbar' beta that leaves D_b alpha,
  ## eta and the objective unchanged (D_b annihilates constants); it reduces the alpha-beta
  ## coupling and hence the number of BCD sweeps. alpha is returned on the original scale.
  xbar <- if (center) colMeans(dat$X) else numeric(p)
  Xs <- sweep(Xs, 2, xbar)
  dk <- if (k == 1L) dat$d1 else dat$d2
  D <- getD(J, b); DtD <- crossprod(D); m <- nrow(D)

  if (is.null(alpha0)) {                                   # crude interval rates
    ek <- tabulate(l[dk == 1L], nbins = J); nr <- tabulate(l, nbins = J)
    q  <- pmin(pmax((ek + .5) / (nr + 1), 1e-6), .45)
    alpha0 <- if (link == "logit") qlogis(q) else log(q)
  }
  beta  <- if (is.null(beta0)) numeric(p) else beta0
  alpha <- alpha0 + sum(xbar * beta)                       # to the centred scale
  if (is.null(gamma)) {
    a <- rep(log(max(sum(dk), 1) / length(dk)), J)
    if (link == "logit") a <- a - log1p(-exp(a[1]))
    g0 <- max(stats::median(gsum(link_terms(a[l], dc, dk, link)$g2, l, J)), 1e-3)
    gamma <- g0 * max(1, cvec[b] * rho)
  }
  z <- as.vector(D %*% alpha); u <- numeric(m)
  ## optional dual warm start. u0 = "kkt": the u that satisfies the alpha-stationarity condition
  ## grad_alpha f + D'u = 0 at the starting point (least squares), clipped to [-rho, rho].
  ## Only the starting point of the iteration changes, not the algorithm or its limit.
  if (identical(u0, "kkt")) {
    g1 <- link_terms(alpha[l] + as.vector(Xs %*% beta)[id], dc, dk, link)$g1
    u  <- -as.vector(solve(tcrossprod(D), D %*% gsum(g1, l, J)))
    u  <- pmin(pmax(u, -rho), rho)
  } else if (is.numeric(u0)) u <- u0
  ## log-link aggregates: events per interval and summed (centred) covariates of the events
  dkl <- gsum(dk, l, J); evx <- colSums(Xs[id[dk == 1L], , drop = FALSE])
  subj_sum <- function(v) { o <- numeric(dat$n); r <- rowsum(v, id, reorder = TRUE)
                            o[as.integer(rownames(r))] <- r; o }
  n_sweep <- 0L
  soft <- function(v, kappa) sign(v) * pmax(abs(v) - kappa, 0)
  eta_of <- function(al, be) al[l] + as.vector(Xs %*% be)[id]
  Lx <- function(al, be) {                       # primal objective; logit: +Inf outside Theta
    eta <- eta_of(al, be)
    if (link == "logit" && max(eta) > 0) return(Inf)
    r <- as.vector(D %*% al) - z
    obj_f(eta, dc, dk, link) + sum(u * r) + gamma / 2 * sum(r^2)
  }
  grads <- function(al, be) {                   # gradient pieces and objective at (al, be)
    eta <- eta_of(al, be); tt <- link_terms(eta, dc, dk, link)
    r <- as.vector(D %*% al) - z
    list(tt = tt, L = obj_f(eta, dc, dk, link) + sum(u * r) + gamma / 2 * sum(r^2),
         ga = gsum(tt$g1, l, J) + as.vector(crossprod(D, u + gamma * (D %*% al - z))),
         gb = as.vector(crossprod(Xs, rowsum(tt$g1, id, reorder = TRUE))))
  }
  ## Armijo backtracking. Near the solution the predicted decrease falls below the rounding level
  ## of the objective, and the comparison is then decided by floating-point noise, which differs
  ## across platforms (e.g. x86 vs Apple arm64) and can reject every step, stalling the sweeps.
  ## Two safeguards: (i) if the predicted decrease is itself below the rounding level, the full
  ## Newton step is taken (if feasible); (ii) the sufficient-decrease test allows a slack of the
  ## size of the rounding error, tolL = 1e-10 * (1 + |L0|). Both only affect steps whose effect
  ## on the objective is below numerical resolution; the stopping rule (||g||_inf) is unchanged.
  ls_step <- function(Lt, L0, slope) {
    tolL <- 1e-10 * (1 + abs(L0))
    if (slope <= tolL && is.finite(Lt(1))) return(1)
    t <- 1
    while (t > 1e-10 && Lt(t) > L0 - 1e-4 * t * slope + tolL) t <- t / 2
    t
  }
  backtrack <- function(al, be, da, db, g, L0)
    ls_step(function(t) Lx(al - t * da, be - t * db), L0, sum(g * c(da, db)))

  for (s in seq_len(maxit)) {
    ## ---- primal step: BCD sweeps (one Newton step on alpha, then one on beta), repeated
    ##      until the full gradient of the primal subproblem satisfies ||g||_inf < newton_tol
    if (link == "log") {
      ## log link: f = sum_l A_l(beta) e^{alpha_l} - sum_l d_l alpha_l - sum_events x' beta
      ##                = sum_i e^{x_i' beta} C_i(alpha) - ...,  so each block needs only interval-
      ##           or subject-level aggregates (one cheap pass over the records per block).
      for (sw in seq_len(max_sweep)) {
        eb <- exp(as.vector(Xs %*% beta)); A <- gsum(dc * eb[id], l, J)
        C  <- subj_sum(dc * exp(alpha)[l])
        ga <- A * exp(alpha) - dkl + as.vector(crossprod(D, u + gamma * (D %*% alpha - z)))
        gb <- as.vector(crossprod(Xs, eb * C)) - evx
        if (max(abs(ga), abs(gb)) < newton_tol) break
        ## alpha-step (A fixed)
        La <- function(al) { r <- as.vector(D %*% al) - z
                             sum(A * exp(al) - dkl * al) + sum(u * r) + gamma / 2 * sum(r^2) }
        da <- as.vector(solve(diag(A * exp(alpha) + 1e-12, J) + gamma * DtD, ga))
        t <- ls_step(function(t) La(alpha - t * da), La(alpha), sum(ga * da))
        alpha <- alpha - t * da
        ## beta-step (C fixed)
        C  <- subj_sum(dc * exp(alpha)[l])
        Lb <- function(be) sum(exp(as.vector(Xs %*% be)) * C) - sum(evx * be)
        w  <- eb * C
        gb <- as.vector(crossprod(Xs, w)) - evx
        db <- as.vector(solve(crossprod(Xs * w, Xs), gb))
        t <- ls_step(function(t) Lb(beta - t * db), Lb(beta), sum(gb * db))
        beta <- beta - t * db
        n_sweep <- n_sweep + 1L
      }
    } else for (sw in seq_len(max_sweep)) {
      G <- grads(alpha, beta)
      if (max(abs(G$ga), abs(G$gb)) < newton_tol) break
      da <- as.vector(solve(diag(gsum(G$tt$g2, l, J) + 1e-12, J) + gamma * DtD, G$ga))
      t  <- backtrack(alpha, beta, da, numeric(p), c(G$ga, numeric(p)), G$L)
      alpha <- alpha - t * da
      G  <- grads(alpha, beta)
      w  <- as.vector(rowsum(G$tt$g2, id, reorder = TRUE))
      db <- as.vector(solve(crossprod(Xs * w, Xs), G$gb))
      t  <- backtrack(alpha, beta, numeric(J), db, c(numeric(J), G$gb), G$L)
      beta <- beta - t * db
      n_sweep <- n_sweep + 1L
    }
    ## ---- z- and u-updates
    Da <- as.vector(D %*% alpha); z_old <- z
    z  <- soft(Da + u / gamma, rho / gamma)
    u  <- u + gamma * (Da - z)
    r_pri <- sqrt(sum((Da - z)^2))
    r_dua <- gamma * sqrt(sum(crossprod(D, z - z_old)^2))
    e_pri <- sqrt(m) * tol_abs + tol_rel * max(sqrt(sum(Da^2)), sqrt(sum(z^2)))
    e_dua <- sqrt(J) * tol_abs + tol_rel * sqrt(sum(crossprod(D, u)^2))
    if (trace > 0 && s %% trace == 0)
      message(sprintf("iter %6d  sweeps %7d  r_pri/e_pri %.3g  r_dua/e_dua %.3g",
                      s, n_sweep, r_pri / e_pri, r_dua / e_dua))
    if (r_pri <= e_pri && r_dua <= e_dua) break
  }
  eta  <- eta_of(alpha, beta)
  fval <- obj_f(eta, dc, dk, link)
  df   <- sum(abs(z) > 1e-8) + b + p
  alpha <- alpha - sum(xbar * beta)                        # back to the original scale
  list(alpha = alpha, beta = beta, z = z, u = u, gamma = gamma, iter = s, n_sweep = n_sweep,
       converged = (r_pri <= e_pri && r_dua <= e_dua), r_pri = r_pri, r_dua = r_dua,
       max_mu = if (link == "logit") plogis(max(eta)) else exp(max(eta)),
       f = fval, df = df, bic = 2 * fval + log(dat$n) * df)
}

## path over a decreasing rho grid with warm starts
tf_path_bcd <- function(dat, k, b, rho_grid, link = "log", ...) {
  a0 <- NULL; b0 <- NULL; out <- vector("list", length(rho_grid))
  for (i in seq_along(rho_grid)) {
    out[[i]] <- admm_bcd_exact(dat, k, b, rho_grid[i], link = link, alpha0 = a0, beta0 = b0, ...)
    a0 <- out[[i]]$alpha; b0 <- out[[i]]$beta
  }
  bic <- sapply(out, `[[`, "bic"); j <- which.min(bic)
  list(rho_grid = rho_grid, fits = out, bic = bic, best = j, rho = rho_grid[j],
       alpha = out[[j]]$alpha, beta = out[[j]]$beta, df = out[[j]]$df,
       n_nonconv = sum(!sapply(out, `[[`, "converged")))
}
