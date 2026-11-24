################################################################################
## sim_cif_figure.R -- figure for Section 4.3 from sim_cif_summary.csv
##   Rscript sim_cif_figure.R            -> fig_sim_cif.pdf
## x = (0,0). (a) relative bias (%) and (b) coverage of the 95% interval, for the proposed
## predictor (black) and Lee et al. (2018) (red). Each column of points is one setting and
## time point; the four points in it are the two scenarios x two causes (circle = scenario 1,
## triangle = scenario 2; filled = cause 1, open = cause 2). In (b) the dotted lines are the
## Monte Carlo limits 0.95 +- 1.96 sqrt(0.95 * 0.05 / 1000).
################################################################################
T <- read.csv("sim_cif_summary.csv")
T <- T[T$profile == "x00", ]
SET <- c(s1n300 = "Large hazards\nn = 300", s1n600 = "Large hazards\nn = 600",
         s02J10 = "Small hazards\nJ = 10", s02J30 = "Small hazards\nJ = 30")
T$pos  <- 2 * (match(T$setting, names(SET)) - 1) + ifelse(T$m == T$J, 2, 1)   # 1..8
T$xoff <- ifelse(T$method == "Proposed", -.15, .15) +
          ifelse(T$dgp == "ph", -.05, .05)
T$pch  <- ifelse(T$dgp == "ph", ifelse(T$cause == 1, 16, 1), ifelse(T$cause == 1, 17, 2))
T$col  <- ifelse(T$method == "Proposed", "black", "#D7301F")

pdf("fig_sim_cif.pdf", width = 9, height = 4.4)
par(mfrow = c(1, 2), mar = c(6.2, 4.4, 2.2, 0.8), oma = c(0, 0, 0, 0))
panel <- function(y, ylab, ylim, ref, lab) {
  plot(T$pos + T$xoff, y, xlim = c(.5, 8.5), ylim = ylim, xaxt = "n", las = 1,
       pch = T$pch, col = T$col, xlab = "", ylab = ylab, main = lab, cex = .9)
  abline(v = c(2.5, 4.5, 6.5), col = "grey80")
  abline(h = ref, lty = 3)
  axis(1, at = 1:8, labels = rep(c("mid", "end"), 4), cex.axis = .8)
  mtext(SET, side = 1, at = c(1.5, 3.5, 5.5, 7.5), line = 3.6, cex = .62)
}
panel(T$rel_bias, "Relative bias (%)", range(c(0, T$rel_bias)) + c(-1, 1), 0, "(a) Relative bias")
legend("topright", bty = "n", cex = .7, pch = c(15, 15, 16, 17, 1),
       col = c("black", "#D7301F", "grey30", "grey30", "grey30"),
       legend = c("Proposed", "Lee et al. (2018)", "Scenario 1, cause 1", "Scenario 2, cause 1",
                  "Open symbol: cause 2"))
mc <- 1.96 * sqrt(.95 * .05 / 1000)
panel(T$CP, "Coverage probability", c(min(.7, min(T$CP)) - .01, 1), c(.95 - mc, .95, .95 + mc),
      "(b) Coverage of 95% intervals")
invisible(dev.off())
cat("wrote fig_sim_cif.pdf\n")
