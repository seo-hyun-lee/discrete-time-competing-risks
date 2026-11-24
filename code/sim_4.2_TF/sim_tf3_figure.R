################################################################################
## sim_tf3_figure.R -- baseline figure for the trend-filtering simulation (Section 4.2)
##   Rscript sim_tf3_figure.R        (uses sim_tf3_{smooth,step}_n1000_J30.rds)
## One 2 x 2 figure, cause 1 only -> fig_baseline.pdf
##   rows: smooth, step baseline; columns: first scenario (log link, log-hazard),
##         second scenario (logit link, log-odds)
## Lines: true baseline (black), unpenalized (red), trend filtering with b and rho chosen by BIC
## (blue); solid = pointwise median over replications, shaded = pointwise 2.5% and 97.5%
## quantiles. The true step baseline is drawn as a step function, so its jumps are vertical. Unpenalized estimates of -Inf (intervals without cause-1 events) are drawn at
## the bottom of the panel.
################################################################################
SHAPES <- c(smooth = "Smooth", step = "Step")
N <- as.integer(Sys.getenv("N", "1000")); J <- as.integer(Sys.getenv("J", "30")); K <- 1
files <- sprintf("sim_tf3_%s_n%d_J%d.rds", names(SHAPES), N, J)
if (!all(file.exists(files))) stop("need ", paste(files, collapse = ", "))
S <- setNames(lapply(files, readRDS), names(SHAPES))

q3 <- function(m) apply(m, 2, quantile, c(.025, .5, .975), na.rm = TRUE, names = FALSE)
pdf("fig_baseline.pdf", width = 9, height = 7.6)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))
panel <- 0
for (s in names(SHAPES)) for (dg in c("ph", "logit")) {
  panel <- panel + 1
  R <- S[[s]]$res[[dg]]; a0 <- S[[s]]$alp0[[dg]][[K]]
  ok  <- apply(R$unpen, 1, function(a) !all(is.na(a)))
  un  <- matrix(R$unpen[ok, , K], nrow = sum(ok))
  sel <- apply(R$bic[ok, K, , drop = FALSE], 1, which.min)
  tf  <- t(sapply(seq_along(sel), function(r) R$tf[which(ok)[r], , K, sel[r]]))
  qu <- q3(un); qt <- q3(tf)
  fin <- c(a0, qt, qu[is.finite(qu)])
  yl  <- c(floor(min(fin)) - .5, ceiling(max(fin)) + .5)
  qu[!is.finite(qu) | qu < yl[1]] <- yl[1] + .05                # -Inf drawn at the bottom
  ## estimates: interval values joined by lines, percentile band filled; the true step baseline
  ## is drawn flat over each interval (centred at the interval index), so its jumps are vertical
  Jx <- ncol(un); x <- seq_len(Jx)
  band <- function(lo, hi, col) polygon(c(x, rev(x)), c(lo, rev(hi)),
                                        col = adjustcolor(col, .22), border = NA)
  plot(NA, xlim = c(.5, Jx + .5), ylim = yl, xlab = "Interval", las = 1, xaxs = "i",
       ylab = if (dg == "ph") "Baseline log-hazard" else "Baseline log-odds",
       main = sprintf("(%s) %s baseline, %s link", letters[panel], tolower(SHAPES[s]),
                      if (dg == "ph") "log" else "logit"))
  band(qu[1, ], qu[3, ], "#D7301F")
  band(qt[1, ], qt[3, ], "#2B6CB0")
  lines(x, qu[2, ], col = "#D7301F", lwd = 1.5)
  lines(x, qt[2, ], col = "#2B6CB0", lwd = 2)
  if (s == "step") lines(c(rbind(x - .5, x + .5)), rep(a0, each = 2), col = "black", lwd = 2)
  else lines(x, a0, col = "black", lwd = 2)
  if (panel == 1) legend("topleft", bty = "n", cex = .8,
                         legend = c("True", "Unpenalized", "Trend filtering"),
                         col = c("black", "#D7301F", "#2B6CB0"), lwd = c(2, 1.5, 2))
}
invisible(dev.off())
cat("wrote fig_baseline.pdf\n")
