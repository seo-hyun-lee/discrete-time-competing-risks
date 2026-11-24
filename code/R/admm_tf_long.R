################################################################################
## admm_tf_long.R
## Discrete-time competing risks: Poisson partial likelihood + trend filtering
## Long (person-period) format, base R only.
##   - log  link : two-step  (beta from partial likelihood, then alpha + TF)
##   - logit link: joint ADMM (beta, alpha) + TF
################################################################################

## ------------------------------------------------- safe per-interval accumulator
gsum <- function(v, l, J) {
  out <- if (is.matrix(v)) matrix(0, J, ncol(v)) else numeric(J)
  o <- rowsum(v, l, reorder = TRUE)
  idx <- as.integer(rownames(o))
  if (is.matrix(v)) out[idx, ] <- o else out[idx] <- o
  out
}

## ---------------------------------------------------------------- difference
getD <- function(J, b) diff(diag(J), differences = b)   # (J-b) x J

## ----------------------------------------------------------------- long data
## Y in 1..J, D in {0,1,2} (0 = censored), X: n x p
make_long <- function(Y, D, X) {
  reps <- as.integer(Y)
  id   <- rep.int(seq_along(Y), reps)
  l    <- sequence(reps)
  last <- (l == Y[id])
  d1   <- as.integer(last & D[id] == 1L)
  d2   <- as.integer(last & D[id] == 2L)
  list(id = id, l = l, d1 = d1, d2 = d2, dc = 1L - d1 - d2,
       X = X[id, , drop = FALSE], n = length(Y), J = max(Y), p = ncol(X))
}

## ------------------------------------------------- partial likelihood (log link)
## Route A: reduced risk set {dc == 1}; Breslow-type baseline
fit_pl <- function(dat, k, maxit = 50, tol = 1e-10) {
  dk <- if (k == 1L) dat$d1 else dat$d2
  J  <- dat$J; p <- dat$p
  r  <- dat$dc == 1L                       # reduced risk set rows
  Xr <- dat$X[r, , drop = FALSE]; lr <- dat$l[r]
  ev <- dat$l[dk == 1L]                    # event times (cause k)
  e  <- tabulate(ev, nbins = J)            # events per interval
  sx <- colSums(dat$X[dk == 1L, , drop = FALSE])
  beta <- numeric(p)
  for (it in 1:maxit) {
    w  <- exp(as.vector(Xr %*% beta))
    S0 <- gsum(w, lr, J)
    S1 <- gsum(Xr * w, lr, J)
    xb <- S1 / ifelse(S0 > 0, S0, NA)                        # J x p
    xb[is.na(xb)] <- 0
    U  <- sx - colSums(e * xb)
    I  <- matrix(0, p, p)
    for (a in 1:p) for (bq in a:p) {
      S2 <- gsum(Xr[, a] * Xr[, bq] * w, lr, J)
      q  <- ifelse(S0 > 0, S2 / ifelse(S0 > 0, S0, 1) - xb[, a] * xb[, bq], 0)
      I[a, bq] <- I[bq, a] <- sum(e * q)
    }
    step <- solve(I, U); beta <- beta + step
    if (max(abs(step)) < tol) break
  }
  S0 <- gsum(exp(as.vector(Xr %*% beta)), lr, J)
  list(beta = beta, I = I, lam0 = e / S0, S0 = S0, events = e)
}

## ------------------------------------------------------------------- link bits
## returns mu, dmu/deta, and the pieces of  f = sum[ dc*mu - dk*log(mu) ]
link_terms <- function(eta, dc, dk, link) {
  if (link == "log") {
    mu <- exp(eta)
    g1 <- dc * mu - dk                  # df/deta
    g2 <- dc * mu                       # d2f/deta2   (>= 0)
  } else {                              # logit
    mu <- 1 / (1 + exp(-eta)); v <- mu * (1 - mu)
    g1 <- dc * v - dk * (1 - mu)
    g2 <- dc * v * (1 - 2 * mu) + dk * v
    g2 <- pmax(g2, 1e-8)                # safeguard (mu < .5 => already > 0)
  }
  list(mu = mu, g1 = g1, g2 = g2)
}

obj_f <- function(eta, dc, dk, link) {
  mu <- if (link == "log") exp(eta) else 1 / (1 + exp(-eta))
  sum(dc * mu) - sum(dk[dk == 1] * log(mu[dk == 1]))
}

## ------------------------------------------------------------------ ADMM core
## beta_fix = NULL -> joint update (logit route);  numeric -> two-step (log route)
admm_tf <- function(dat, k, link = "log", b = 2, rho = 0, gamma = NULL,
                    beta_fix = NULL, alpha0 = NULL, beta0 = NULL,
                    maxit = 2000, tol_abs = 1e-9, tol_rel = 1e-7,
                    adapt = TRUE) {
  J <- dat$J; p <- dat$p
  dk <- if (k == 1L) dat$d1 else dat$d2
  dc <- dat$dc; lv <- dat$l; X <- dat$X
  D  <- getD(J, b); DtD <- crossprod(D); m <- nrow(D)

  alpha <- if (is.null(alpha0)) rep(log(max(sum(dk), 1) / length(dk)), J) else alpha0
  if (link == "logit" && is.null(alpha0)) alpha <- alpha - log1p(-exp(alpha[1]))
  beta  <- if (!is.null(beta_fix)) beta_fix else if (is.null(beta0)) numeric(p) else beta0
  z <- as.vector(D %*% alpha); u <- numeric(m)
  joint <- is.null(beta_fix)
  if (is.null(gamma)) {            # scale the augmented parameter to the Hessian
    h0 <- gsum(link_terms(alpha[dat$l] + as.vector(X %*% beta), dc, dk, link)$g2,
               lv, J)
    gamma <- max(stats::median(h0), 1e-3)
  }

  ## cached person-period quantities (recomputed only when beta changes)
  xb <- as.vector(X %*% beta)
  ## log link with fixed beta: the interval-l aggregates of the person-period
  ## likelihood do not change across iterations, so cache them once.
  fast <- (link == "log") && !joint
  if (fast) { A <- gsum(dc * exp(xb), lv, J); ek <- gsum(dk, lv, J) }

  for (s in 1:maxit) {
    ## --- alpha block: one Newton step -------------------------------------
    if (fast) {
      ea <- A * exp(alpha); ga <- ea - ek; ha <- ea
    } else {
      tt  <- link_terms(alpha[dat$l] + xb, dc, dk, link)
      ga  <- gsum(tt$g1, lv, J)
      ha  <- gsum(tt$g2, lv, J)
    }
    H   <- diag(ha + 1e-10, J) + gamma * DtD
    g   <- ga + as.vector(crossprod(D, u)) + gamma * as.vector(crossprod(D, D %*% alpha - z))
    alpha <- alpha - as.vector(solve(H, g))

    ## --- beta block: one Newton step --------------------------------------
    if (joint) {
      tt <- link_terms(alpha[dat$l] + xb, dc, dk, link)
      Ub <- as.vector(crossprod(X, tt$g1))
      Hb <- crossprod(X * tt$g2, X)
      beta <- beta - as.vector(solve(Hb + diag(1e-10, p), Ub))
      xb <- as.vector(X %*% beta)
    }

    ## --- z and u ----------------------------------------------------------
    Da <- as.vector(D %*% alpha)
    zo <- z
    v  <- Da + u / gamma
    z  <- sign(v) * pmax(abs(v) - rho / gamma, 0)
    u  <- u + gamma * (Da - z)

    if (!all(is.finite(alpha)) || !all(is.finite(beta)) || max(abs(alpha)) > 1e6)
      return(list(alpha = alpha, beta = beta, z = z, u = u, iter = s, df = NA,
                  f = NA, converged = FALSE, diverged = TRUE))
    r_pri <- sqrt(sum((Da - z)^2))
    r_dua <- gamma * sqrt(sum((crossprod(D, z - zo))^2))
    e_pri <- sqrt(m) * tol_abs + tol_rel * max(sqrt(sum(Da^2)), sqrt(sum(z^2)))
    e_dua <- sqrt(J) * tol_abs + tol_rel * sqrt(sum((crossprod(D, u))^2))
    if (r_pri <= e_pri && r_dua <= e_dua) break
    if (adapt && s %% 5 == 0) {                      # residual balancing
      if (r_pri > 10 * r_dua) gamma <- gamma * 2
      else if (r_dua > 10 * r_pri) gamma <- gamma / 2
    }
  }
  eta <- alpha[dat$l] + as.vector(X %*% beta)
  df  <- sum(abs(z) > 1e-8) + b + (if (joint) p else 0)
  list(alpha = alpha, beta = beta, z = z, u = u, iter = s, df = df,
       f = obj_f(eta, dc, dk, link), r_pri = r_pri, r_dua = r_dua,
       converged = (s < maxit), diverged = FALSE)
}

## ---------------------------------------------------------------- BIC on a path
## nbic = number of independent subjects (clustered rows)
tf_path <- function(dat, k, link = "log", b = 2, beta_fix = NULL,
                    rho_grid = NULL, gamma = NULL, nrho = 30, ...) {
  if (is.null(rho_grid)) rho_grid <- exp(seq(log(1e-3), log(5e2), length = nrho))
  rho_grid <- sort(rho_grid, decreasing = TRUE)      # warm start from smoothest
  fits <- vector("list", length(rho_grid)); bic <- numeric(length(rho_grid))
  a0 <- NULL; b0 <- NULL
  for (i in seq_along(rho_grid)) {
    fit <- admm_tf(dat, k, link, b, rho = rho_grid[i], gamma = gamma,
                   beta_fix = beta_fix, alpha0 = a0, beta0 = b0, ...)
    if (isTRUE(fit$diverged)) { fits[[i]] <- fit; bic[i] <- Inf; a0 <- NULL; b0 <- NULL; next }
    a0 <- fit$alpha; b0 <- if (is.null(beta_fix)) fit$beta else NULL
    fits[[i]] <- fit
    bic[i] <- 2 * fit$f + log(dat$n) * fit$df
  }
  i <- which.min(bic)
  list(fit = fits[[i]], rho = rho_grid[i], bic = bic, rho_grid = rho_grid,
       df = fits[[i]]$df, fits = fits)
}

## ------------------------------------------------------------------------ DGPs
## h1, h2 : length-J baseline hazard vectors (continuous-time rates for "ph")
gen_data <- function(n, J, h1, h2, beta1, beta2, dgp = c("ph", "logit"),
                     exhaust = FALSE) {
  dgp <- match.arg(dgp)
  X <- cbind(rnorm(n), rbinom(n, 1, .5))
  e1 <- exp(as.vector(X %*% beta1)); e2 <- exp(as.vector(X %*% beta2))
  Y <- rep(J, n); D <- rep(0L, n); alive <- rep(TRUE, n)
  for (t in 1:J) {
    idx <- which(alive); if (!length(idx)) break
    if (dgp == "ph") {
      a1 <- h1[t] * e1[idx]; a2 <- h2[t] * e2[idx]
      pev <- 1 - exp(-(a1 + a2)); pk1 <- a1 / (a1 + a2)
      ev <- runif(length(idx)) < pev
      k  <- ifelse(runif(length(idx)) < pk1, 1L, 2L)
    } else {
      l1 <- plogis(qlogis(h1[t]) + as.vector(X[idx, ] %*% beta1))
      l2 <- plogis(qlogis(h2[t]) + as.vector(X[idx, ] %*% beta2))
      uu <- runif(length(idx))
      ev <- uu < (l1 + l2); k <- ifelse(uu < l1, 1L, 2L)
    }
    if (exhaust && t == J) { ev <- rep(TRUE, length(idx)) }   # terminal exhaustion
    hit <- idx[ev]
    Y[hit] <- t; D[hit] <- k[ev]; alive[hit] <- FALSE
  }
  list(Y = Y, D = D, X = X)
}

## true alpha on the model scale used by the fit
true_alpha <- function(h1, h2, k, dgp) {
  if (dgp == "logit") return(qlogis(if (k == 1) h1 else h2))
  tot <- h1 + h2; pev <- 1 - exp(-tot)
  lam <- (if (k == 1) h1 else h2) / tot * pev          # discrete hazard at x = 0
  log(lam)
}

## ---------------------------------------------- Lee et al. (2018) comparator
## collapsed binomial GLM with unstructured baseline (treat competing event as censored)
fit_lee <- function(dat, k) {
  dk <- if (k == 1L) dat$d1 else dat$d2
  fl <- factor(dat$l, levels = 1:dat$J)
  g <- suppressWarnings(glm.fit(cbind(model.matrix(~ fl - 1), dat$X), dk,
                                family = binomial()))
  cf <- g$coefficients
  list(alpha = cf[1:dat$J], beta = cf[-(1:dat$J)])
}
