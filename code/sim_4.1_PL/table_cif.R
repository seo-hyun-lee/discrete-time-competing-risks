################################################################################
## table_cif.R -- summarise sim_cif_*.rds (output of sim_cif.R)
##   Usage: source("table_cif.R")   or   Rscript table_cif.R [files...]
## Outputs: cif_summary.csv and one LaTeX table per (hazard scale, J):
## table_cif_scale<SCALE>_J<J>.tex, label tab:cif_s<SCALE>_J<J>. Both covariate
## profiles are shown side by side; rows = n x (cause, t_j) x {Bias, ASE, ESE, CP};
## Scenario 1 = grouped PH, Scenario 2 = logit; columns Lee et al. (2018) and Proposed.
## Bias = mean(F-hat) - F;  ASE = mean SE;  ESE = SD of F-hat;
## CP = coverage of the TRUE F by the 95% Wald interval.
################################################################################
files <- commandArgs(TRUE)
if (!length(files)) files <- list.files(pattern = "^sim_cif_s.*_n.*\\.rds$")
if (!length(files)) stop("no sim_cif_*.rds found in ", getwd())
PROFILE <- Sys.getenv("PROFILE", "x=(0,0)")

summ <- function(r) {
  one <- function(F, se) c(bias = mean(F) - r$truth[1], ASE = mean(se), ESE = sd(F),
                           CP = mean(abs(F - r$truth[1]) <= 1.96 * se))
  rbind(Lee = one(r$F_lee, r$se_lee), Proposed = one(r$F_prop, r$se_prop)) }
rows <- list()
for (f in files) {
  S <- readRDS(f); R <- S$R
  for (key in unique(paste(R$dgp, R$profile, R$cause, R$j))) {
    r <- R[paste(R$dgp, R$profile, R$cause, R$j) == key, ]; m <- summ(r)
    for (me in rownames(m)) rows[[length(rows) + 1]] <- data.frame(
      scale = S$scale, n = S$n, J = S$J, B = S$B, dgp = r$dgp[1], profile = r$profile[1],
      cause = r$cause[1], j = r$j[1], truth = r$truth[1], method = me,
      bias = m[me, "bias"], ASE = m[me, "ASE"], ESE = m[me, "ESE"], CP = m[me, "CP"],
      row.names = NULL) } }
L <- do.call(rbind, rows); write.csv(L, "cif_summary.csv", row.names = FALSE)
cat("wrote cif_summary.csv (", nrow(L), "rows )\n")

f4 <- function(v) ifelse(v < 0, sprintf("$-$%.4f", -v), sprintf("%.4f", v))
combos <- unique(L[, c("scale", "J")]); combos <- combos[order(combos$scale, combos$J), ]
for (ci in seq_len(nrow(combos))) {
  sc <- combos$scale[ci]; Jc <- combos$J[ci]
  T1 <- L[L$scale == sc & L$J == Jc, ]
  profs <- sort(unique(T1$profile)); np <- length(profs)
  ns <- sort(unique(T1$n)); Bn <- T1$B[1]
  sct <- gsub("\\.", "", format(sc)); lab0 <- sprintf("cif_s%s_J%d", sct, Jc)
  ptex <- sapply(profs, function(pr) sprintf("$x=%s^{\\top}$", sub("^x=", "", pr)))
  hdr1 <- paste0(" &  & ", paste(sprintf("& \\multicolumn{4}{c}{%s}", ptex), collapse = " "), " \\\\")
  cm1 <- paste(sprintf("\\cmidrule(lr){%d-%d}", 4 + 4 * (seq_len(np) - 1), 7 + 4 * (seq_len(np) - 1)), collapse = "")
  hdr2 <- paste0(" &  & ", paste(rep("& \\multicolumn{2}{c}{Scenario 1} & \\multicolumn{2}{c}{Scenario 2}", np), collapse = " "), " \\\\")
  cm2 <- paste(sprintf("\\cmidrule(lr){%d-%d}", 4 + 2 * (seq_len(2 * np) - 1), 5 + 2 * (seq_len(2 * np) - 1)), collapse = "")
  hdr3 <- paste0("$n$ & Cumulative incidence & ", paste(rep("& Lee & Proposed & Lee & Proposed", np), collapse = " "), " \\\\")
  tex <- c("\\begin{table}[htbp]", "\\centering", "\\footnotesize",
    sprintf("\\caption{Prediction of the cumulative incidence with baseline hazards multiplied by %s, based on $B=%d$ replications with $J=%d$ discrete time points. Lee: predictor of \\citet{lee2018analysis}; Proposed: predictor of Section~2.3. Bias: mean of the estimate minus the true value; ASE: average estimated standard error; ESE: empirical standard error; CP: empirical coverage probability of the true value by the nominal $95\\%%$ Wald-type confidence interval.}",
            format(sc), Bn, Jc),
    sprintf("\\label{tab:%s}", lab0), "\\setlength{\\tabcolsep}{3pt}",
    sprintf("\\begin{tabular}{ll l %s}", paste(rep("rr rr", np), collapse = " ")), "\\toprule",
    hdr1, cm1, hdr2, cm2, hdr3, "\\midrule")
  for (ni in seq_along(ns)) {
    if (ni > 1) tex <- c(tex, "\\midrule")
    first <- TRUE
    for (k in 1:2) for (jj in sort(unique(T1$j))) {
      if (!first) tex <- c(tex, "\\addlinespace")
      for (st in c("bias", "ASE", "ESE", "CP")) {
        fm <- if (st == "bias") f4 else if (st == "CP") function(z) sprintf("%.3f", z) else function(z) sprintf("%.4f", z)
        cells <- c()
        for (pr in profs) for (dg in c("ph", "logit")) for (me in c("Lee", "Proposed"))
          cells <- c(cells, fm(T1[[st]][T1$n == ns[ni] & T1$profile == pr & T1$dgp == dg & T1$cause == k & T1$j == jj & T1$method == me]))
        lab <- if (st == "bias") sprintf("$F_%d(t_{%d})$", k, jj) else ""
        nn <- if (first && st == "bias") format(ns[ni]) else ""
        tex <- c(tex, paste(nn, "&", lab, "&", c(bias = "Bias", ASE = "ASE", ESE = "ESE", CP = "CP")[st], "&",
                            paste(cells, collapse = " & "), "\\\\"))
        first <- FALSE } } }
  tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  fn <- sprintf("table_cif_scale%s_J%d.tex", sct, Jc); writeLines(tex, fn)
  cat("wrote ", fn, " (label tab:", lab0, ")\n", sep = "") }
