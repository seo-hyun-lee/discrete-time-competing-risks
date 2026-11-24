################################################################################
## sim_tf3_tables.R -- LaTeX tables and the numbers quoted in Section 4.2 (revised simulation)
##   Rscript sim_tf3_tables.R
## Reads sim_tf3_{smooth,step,abrupt}_n1000_J30.rds and, if present,
##       sim_tf3_smooth_n500_J30.rds, sim_tf3_smooth_n2000_J30.rds.
## Writes
##   sim_tf3_tables.tex   tab:ise (n = 1000), tab:ise_n (smooth, by n; if the files exist),
##                        tab:ise_split (interior / tail), tab:beta_tf (regression coefficients)
##   sim_tf3_numbers.txt  every number quoted in the text, with how it was computed
## Conventions (as in the text): hazard scale lambda = g^{-1}(alpha) at x = 0, an interval
## without cause-k events contributes lambda_hat = 0; ISE_k = sum_l (lambda_hat - lambda)^2,
## reported x 10^3; interior = intervals l <= 2J/3 (t <= 20 for J = 30), tail = the rest.
## Unpenalized fit "degenerate": some alpha_hat non-finite or below -10.
################################################################################
SHAPES <- c(smooth = "Smooth", step = "Step", abrupt = "Abrupt")
DGPN   <- c(ph = "1 (PH)", logit = "2 (logit)")
rd <- function(n, s) { f <- sprintf("sim_tf3_%s_n%d_J30.rds", s, n); if (file.exists(f)) readRDS(f) }
S1000 <- lapply(setNames(nm = names(SHAPES)), function(s) rd(1000, s))
if (any(sapply(S1000, is.null))) stop("need sim_tf3_{smooth,step,abrupt}_n1000_J30.rds")

inv <- function(dg) if (dg == "ph") exp else plogis
## hazard-scale estimates (reps x J) for every estimator, one (shape, dgp, cause)
est <- function(S, dg, k) {
  R <- S$res[[dg]]; g <- inv(dg)
  ok <- apply(R$unpen, 1, function(a) !all(is.na(a)))
  M <- function(a) matrix(a, nrow = sum(ok))
  sel <- apply(R$bic[ok, k, , drop = FALSE], 1, which.min)
  un <- M(R$unpen[ok, , k])
  out <- list(unpen = g(un), lee = plogis(M(R$lee[ok, , k])))
  for (bo in seq_along(S$BS)) out[[paste0("b", S$BS[bo])]] <- g(M(R$tf[ok, , k, bo]))
  out$sel <- M(t(sapply(seq_along(sel), function(r) g(R$tf[which(ok)[r], , k, sel[r]]))))
  out <- lapply(out, function(m) { m[!is.finite(m)] <- 0; m })
  list(lam = out, lam0 = g(S$alp0[[dg]][[k]]), sel = sel, nb = length(S$BS),
       degen = mean(apply(un, 1, function(a) any(!is.finite(a) | a < -10))) * 100,
       beta = list(unpen = M(R$beta_u[ok, (k - 1) * 2 + 1:2]), lee = M(R$beta_lee[ok, (k - 1) * 2 + 1:2]),
                   sel = M(t(sapply(seq_along(sel), function(r) R$beta_t[which(ok)[r], (k - 1) * 2 + 1:2, sel[r]])))),
       B = sum(ok))
}
ise  <- function(m, l0, idx = seq_along(l0)) rowSums(sweep(m[, idx, drop = FALSE], 2, l0[idx])^2) * 1e3
bv   <- function(m, l0, idx) c(bias2 = sum((colMeans(m[, idx, drop = FALSE]) - l0[idx])^2) * 1e3,
                               var   = sum(apply(m[, idx, drop = FALSE], 2, var)) * 1e3)
f3 <- function(x) formatC(x, format = "f", digits = 3)
NUM <- character()
note <- function(...) NUM <<- c(NUM, sprintf(...))

## ---------------------------------------------------------------- tab:ise (n = 1000)
rows <- character(); RAT <- list()
for (s in names(SHAPES)) {
  first_s <- TRUE
  for (dg in c("ph", "logit")) for (k in 1:2) {
    E <- est(S1000[[s]], dg, k); L <- E$lam
    m <- sapply(L, function(x) mean(ise(x, E$lam0)))
    fixed <- grep("^b", names(m)); best <- fixed[which.min(m[fixed])]
    cells <- f3(m); cells[best] <- paste0("\\textbf{", cells[best], "}")
    pct <- paste(sprintf("%.0f", tabulate(E$sel, E$nb) / length(E$sel) * 100), collapse = " / ")
    rows <- c(rows, sprintf("%s & %s & %d & %s & %s & %s & %s & %s & %s & %.0f \\\\",
      if (first_s) SHAPES[[s]] else "", if (k == 1) DGPN[[dg]] else "", k,
      cells["unpen"], cells["lee"], cells["b1"], cells["b2"], cells["sel"], pct, E$degen))
    first_s <- FALSE
    RAT[[length(RAT) + 1]] <- data.frame(shape = s, dgp = dg, cause = k,
      unpen = m["unpen"], lee = m["lee"], b1 = m["b1"], b2 = m["b2"], sel = m["sel"],
      sel_b1 = mean(E$sel == 1) * 100, degen = E$degen,
      best_fixed = names(m)[best], sel_over_best = m["sel"] / min(m[fixed]),
      pct_best = mean(E$sel == as.integer(sub("b", "", names(m)[best]))) * 100)
    if (k == 2 && dg == "ph") rows <- c(rows, "\\addlinespace")
  }
  rows <- c(rows, "\\midrule")
}
rows[length(rows)] <- "\\bottomrule"
RAT <- do.call(rbind, RAT); rownames(RAT) <- NULL
tab_ise <- c(
"\\begin{table}[htbp]", "\\centering", "\\footnotesize",
"\\caption{Mean integrated squared error ($\\times 10^{3}$) of the estimated baseline hazards, $J=30$, $n=1000$, $B=100$. Bold: smaller of the two fixed orders. ``Selected'': percentage of replicates in which BIC selected $b=1/2$. ``Degenerate'': percentage of replicates in which the unpenalized estimate had a non-finite value or a value below $-10$.}",
"\\label{tab:ise}", "\\setlength{\\tabcolsep}{3.5pt}",
"\\begin{tabular}{lll rr rr r c r}", "\\toprule",
" &  &  &  &  & \\multicolumn{2}{c}{Trend filtering} &  &  &  \\\\", "\\cmidrule(lr){6-7}",
"Baseline & DGP & Cause & \\begin{tabular}[c]{@{}c@{}}Unpen-\\\\alized\\end{tabular} & \\begin{tabular}[c]{@{}c@{}}Lee et al.\\\\(2018)\\end{tabular} & $b=1$ & $b=2$ & \\begin{tabular}[c]{@{}c@{}}BIC-\\\\selected $b$\\end{tabular} & \\begin{tabular}[c]{@{}c@{}}Selected\\\\(\\%)\\end{tabular} & \\begin{tabular}[c]{@{}c@{}}Degenerate\\\\(\\%)\\end{tabular} \\\\",
"\\midrule", rows, "\\end{tabular}", "\\end{table}")

note("== tab:ise (n = 1000) ==")
note("sel < unpen and sel < lee in all 12: %s / %s", all(RAT$sel < RAT$unpen), all(RAT$sel < RAT$lee))
note("both fixed orders < unpen and < lee in all 12: %s", all(pmax(RAT$b1, RAT$b2) < pmin(RAT$unpen, RAT$lee)))
note("unpen / sel: %.1f to %.1f;  lee / sel: %.1f to %.1f", min(RAT$unpen / RAT$sel), max(RAT$unpen / RAT$sel),
     min(RAT$lee / RAT$sel), max(RAT$lee / RAT$sel))
note("degenerate unpenalized: max %.0f%% (%s)", max(RAT$degen),
     paste(with(RAT[which.max(RAT$degen), ], paste(shape, dgp, "cause", cause)), collapse = ""))
note("best fixed order per configuration:"); for (i in seq_len(nrow(RAT)))
  note("  %-6s %-5s cause %d: best %s, BIC picks b=1 in %3.0f%%, picks best in %3.0f%%, sel/best = %.2f",
       RAT$shape[i], RAT$dgp[i], RAT$cause[i], RAT$best_fixed[i], RAT$sel_b1[i], RAT$pct_best[i], RAT$sel_over_best[i])
note("sel within 10%% of best fixed order in %d of 12", sum(RAT$sel_over_best <= 1.10))

## ---------------------------------------------------------------- tab:ise_split + decomposition
J <- length(S1000$smooth$alp0$ph[[1]]); int <- seq_len(round(2 * J / 3)); tl <- setdiff(seq_len(J), int)
rows <- character(); SP <- list()
for (s in names(SHAPES)) {
  first_s <- TRUE
  for (dg in c("ph", "logit")) for (k in 1:2) {
    E <- est(S1000[[s]], dg, k); L <- E$lam
    v <- c(); for (reg in list(int, tl)) for (e in c("unpen", "lee", "sel")) v <- c(v, mean(ise(L[[e]], E$lam0, reg)))
    rows <- c(rows, sprintf("%s & %s & %d & %s \\\\", if (first_s) SHAPES[[s]] else "", if (k == 1) DGPN[[dg]] else "",
                            k, paste(f3(v), collapse = " & ")))
    first_s <- FALSE
    for (rg in c("interior", "tail")) {
      idx <- if (rg == "interior") int else tl
      d_un <- ise(L$unpen, E$lam0, idx) - ise(L$sel, E$lam0, idx)
      d_le <- ise(L$lee, E$lam0, idx) - ise(L$sel, E$lam0, idx)
      SP[[length(SP) + 1]] <- data.frame(shape = s, dgp = dg, cause = k, region = rg,
        d_unpen = mean(d_un), se_unpen = sd(d_un) / sqrt(E$B), d_lee = mean(d_le), se_lee = sd(d_le) / sqrt(E$B),
        bias2_sel = bv(L$sel, E$lam0, idx)[1], var_sel = bv(L$sel, E$lam0, idx)[2],
        bias2_lee = bv(L$lee, E$lam0, idx)[1], var_lee = bv(L$lee, E$lam0, idx)[2],
        bias2_unpen = bv(L$unpen, E$lam0, idx)[1], var_unpen = bv(L$unpen, E$lam0, idx)[2])
    }
    if (k == 2 && dg == "ph") rows <- c(rows, "\\addlinespace")
  }
  rows <- c(rows, "\\midrule")
}
rows[length(rows)] <- "\\bottomrule"
SP <- do.call(rbind, SP); rownames(SP) <- NULL
tab_split <- c(
"\\begin{table}[htbp]", "\\centering", "\\footnotesize",
sprintf("\\caption{Mean integrated squared error ($\\times 10^{3}$) of the estimated baseline hazards over the interior intervals ($t\\le %d$) and the tail intervals ($t\\ge %d$), $J=30$, $n=1000$, $B=100$. ``BIC-$b$'': trend filtering with $b$ and $\\rho$ selected by BIC.}", max(int), min(tl)),
"\\label{tab:ise_split}", "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{lll rrr rrr}", "\\toprule",
" &  &  & \\multicolumn{3}{c}{Interior} & \\multicolumn{3}{c}{Tail} \\\\", "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}",
"Baseline & DGP & Cause & Unpenalized & \\begin{tabular}[c]{@{}c@{}}Lee et al.\\\\(2018)\\end{tabular} & BIC-$b$ & Unpenalized & \\begin{tabular}[c]{@{}c@{}}Lee et al.\\\\(2018)\\end{tabular} & BIC-$b$ \\\\",
"\\midrule", rows, "\\end{tabular}", "\\end{table}")

note("\n== tab:ise_split: paired differences (x 10^3), MC SE = sd/sqrt(B) ==")
note("TF better than unpen by > 2 SE in all 24 region cells: %s", all(SP$d_unpen > 2 * SP$se_unpen))
note("TF better than Lee by > 2 SE in %d of 24 cells; the others:", sum(SP$d_lee > 2 * SP$se_lee))
for (i in which(!(SP$d_lee > 2 * SP$se_lee)))
  note("  %-6s %-5s cause %d %-8s: Lee - TF = %.3f (SE %.3f)", SP$shape[i], SP$dgp[i], SP$cause[i], SP$region[i], SP$d_lee[i], SP$se_lee[i])
note("\n== bias^2 / variance decomposition (x 10^3) ==")
for (rg in c("interior", "tail")) { x <- SP[SP$region == rg, ]
  note("%s: var(Lee)/var(TF) %.1f to %.1f; bias2 TF %.2f to %.2f; bias2 Lee <= %.2f; bias2 unpen %.2f to %.2f",
       rg, min(x$var_lee / x$var_sel), max(x$var_lee / x$var_sel), min(x$bias2_sel), max(x$bias2_sel),
       max(x$bias2_lee), min(x$bias2_unpen), max(x$bias2_unpen)) }
x <- SP[SP$region == "tail", ]; i <- which.max(x$bias2_sel / (x$bias2_sel + x$var_sel))
note("largest tail bias share for TF: %s %s cause %d, %.0f%% of tail ISE", x$shape[i], x$dgp[i], x$cause[i],
     100 * x$bias2_sel[i] / (x$bias2_sel[i] + x$var_sel[i]))

## ---------------------------------------------------------------- tab:ise_n (smooth, by n)
Sn <- list(`500` = rd(500, "smooth"), `1000` = S1000$smooth, `2000` = rd(2000, "smooth"))
Sn <- Sn[!sapply(Sn, is.null)]
tab_n <- character()
if (length(Sn) > 1) {
  rows <- character()
  for (dg in c("ph", "logit")) { for (k in 1:2) for (nn in names(Sn)) {
      E <- est(Sn[[nn]], dg, k); m <- sapply(E$lam, function(x) mean(ise(x, E$lam0)))
      cells <- f3(m); b <- grep("^b", names(m)); bb <- b[which.min(m[b])]; cells[bb] <- paste0("\\textbf{", cells[bb], "}")
      rows <- c(rows, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %.0f \\\\", if (k == 1 && nn == names(Sn)[1]) DGPN[[dg]] else "",
        if (nn == names(Sn)[1]) k else "", nn, cells["unpen"], cells["lee"], cells["b1"], cells["b2"], cells["sel"], E$degen))
      if (k == 1) note("smooth %s cause 1 n=%s: lee/sel = %.2f, degenerate %.0f%%", dg, nn, m["lee"] / m["sel"], E$degen)
    }
    rows <- c(rows, "\\midrule") }
  rows[length(rows)] <- "\\bottomrule"
  tab_n <- c("\\begin{table}[htbp]", "\\centering", "\\footnotesize",
    "\\caption{Mean integrated squared error ($\\times 10^{3}$) of the estimated baseline hazards for the smooth configuration by sample size, $J=30$, $B=100$. Bold: smaller of the two fixed orders.}",
    "\\label{tab:ise_n}", "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{lcr rr rr r r}", "\\toprule",
    " &  &  &  &  & \\multicolumn{2}{c}{Trend filtering} &  &  \\\\", "\\cmidrule(lr){6-7}",
    "DGP & Cause & $n$ & \\begin{tabular}[c]{@{}c@{}}Unpen-\\\\alized\\end{tabular} & \\begin{tabular}[c]{@{}c@{}}Lee et al.\\\\(2018)\\end{tabular} & $b=1$ & $b=2$ & \\begin{tabular}[c]{@{}c@{}}BIC-\\\\selected $b$\\end{tabular} & \\begin{tabular}[c]{@{}c@{}}Degenerate\\\\(\\%)\\end{tabular} \\\\",
    "\\midrule", rows, "\\end{tabular}", "\\end{table}")
} else note("\n(tab:ise_n skipped: sim_tf3_smooth_n500_J30.rds / n2000 not found)")

## ---------------------------------------------------------------- regression coefficients
rows <- character(); BT <- list(); bl <- c("$\\beta_{%d1}$", "$\\beta_{%d2}$")
for (s in names(SHAPES)) { first_s <- TRUE
  for (dg in c("ph", "logit")) for (k in 1:2) {
    E <- est(S1000[[s]], dg, k); tb <- if (k == 1) S1000[[s]]$b1 else S1000[[s]]$b2
    for (j in 1:2) {
      st <- sapply(E$beta, function(m) c(bias = mean(m[, j]) - tb[j], ese = sd(m[, j]), rmse = sqrt(mean((m[, j] - tb[j])^2))))
      rows <- c(rows, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", if (first_s) SHAPES[[s]] else "",
        if (k == 1 && j == 1) DGPN[[dg]] else "", sprintf(bl[j], k),
        f3(st["bias", "unpen"]), f3(st["ese", "unpen"]), f3(st["bias", "sel"]), f3(st["ese", "sel"]),
        f3(st["bias", "lee"]), f3(st["ese", "lee"])))
      first_s <- FALSE
      BT[[length(BT) + 1]] <- data.frame(shape = s, dgp = dg, coef = sprintf("b%d%d", k, j),
        bias_un = st["bias", "unpen"], bias_tf = st["bias", "sel"], ese_un = st["ese", "unpen"], ese_tf = st["ese", "sel"],
        rmse_un = st["rmse", "unpen"], rmse_tf = st["rmse", "sel"], bias_lee = st["bias", "lee"])
    }
    if (k == 2 && dg == "ph") rows <- c(rows, "\\addlinespace")
  }
  rows <- c(rows, "\\midrule") }
rows[length(rows)] <- "\\bottomrule"
BT <- do.call(rbind, BT)
tab_beta <- c("\\begin{table}[htbp]", "\\centering", "\\footnotesize",
  "\\caption{Bias and empirical standard error (ESE) of the regression coefficients, unpenalized versus trend filtering with $b$ and $\\rho$ selected by BIC, $J=30$, $n=1000$, $B=100$. Under the first scenario the unpenalized estimator coincides with the partial-likelihood estimator.}",
  "\\label{tab:beta_tf}", "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{lll rr rr rr}", "\\toprule",
  " &  &  & \\multicolumn{2}{c}{Unpenalized} & \\multicolumn{2}{c}{Trend filtering} & \\multicolumn{2}{c}{Lee et al.\\ (2018)} \\\\",
  "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  "Baseline & DGP & & Bias & ESE & Bias & ESE & Bias & ESE \\\\", "\\midrule", rows, "\\end{tabular}", "\\end{table}")
note("\n== regression coefficients (24 = 3 shapes x 2 DGPs x 4 coefficients) ==")
for (dg in c("ph", "logit")) { x <- BT[BT$dgp == dg, ]
  note("%s: RMSE(TF)/RMSE(unpen) %.3f to %.3f (mean %.3f); |bias| unpen <= %.4f, TF <= %.4f; ESE ratio TF/unpen %.3f to %.3f",
       dg, min(x$rmse_tf / x$rmse_un), max(x$rmse_tf / x$rmse_un), mean(x$rmse_tf / x$rmse_un),
       max(abs(x$bias_un)), max(abs(x$bias_tf)), min(x$ese_tf / x$ese_un), max(x$ese_tf / x$ese_un)) }
x <- BT[BT$dgp == "logit" & BT$coef == "b11", ]
note("logit beta11 bias, unpen -> TF: %s", paste(sprintf("%s %.3f -> %.3f", x$shape, x$bias_un, x$bias_tf), collapse = "; "))
note("logit beta11 bias, Lee: %s", paste(sprintf("%s %.4f", x$shape, x$bias_lee), collapse = "; "))

## ---------------------------------------------------------------- CIF at x = 0
cif <- function(l1, l2) { surv <- t(apply(1 - l1 - l2, 1, cumprod)); prev <- cbind(1, surv[, -ncol(surv), drop = FALSE])
  list(t(apply(l1 * prev, 1, cumsum)), t(apply(l2 * prev, 1, cumsum))) }
note("\n== cumulative incidence at x = 0: mean absolute error over intervals and replicates ==")
for (s in names(SHAPES)) for (dg in c("ph", "logit")) {
  E1 <- est(S1000[[s]], dg, 1); E2 <- est(S1000[[s]], dg, 2)
  tr <- cif(matrix(E1$lam0, 1), matrix(E2$lam0, 1))
  mae <- sapply(c("unpen", "lee", "sel"), function(e) { cc <- cif(E1$lam[[e]], E2$lam[[e]])
    c(sapply(1:2, function(k) mean(abs(sweep(cc[[k]], 2, tr[[k]][1, ]))))) })
  note("%-6s %-5s cause 1: unpen %.4f lee %.4f TF %.4f | cause 2: unpen %.4f lee %.4f TF %.4f  (TF/unpen %.2f, %.2f)",
       s, dg, mae[1, 1], mae[1, 2], mae[1, 3], mae[2, 1], mae[2, 2], mae[2, 3], mae[1, 3] / mae[1, 1], mae[2, 3] / mae[2, 1])
}

writeLines(c("%% generated by sim_tf3_tables.R", tab_ise, "", tab_n, "", tab_split, "", tab_beta), "sim_tf3_tables.tex")
writeLines(NUM, "sim_tf3_numbers.txt")
cat(NUM, sep = "\n"); cat("\nwrote sim_tf3_tables.tex and sim_tf3_numbers.txt\n")
