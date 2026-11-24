################################################################################
## sim_tf3_collect.R -- assemble the replications saved so far into one .rds (for a first look
## while sim_tf3.R is still running; does not touch the running processes)
##   SHAPE=smooth Rscript sim_tf3_collect.R      -> sim_tf3_smooth_n1000_J30_partial.rds
##   Rscript sim_tf3_summary.R sim_tf3_smooth_n1000_J30_partial.rds
## Settings must match the run (defaults: N=1000, J=30, B=100, RHO_MIN=1, GAMMA_C="1,1").
################################################################################
source("admm_tf_long.R")
n <- as.integer(Sys.getenv("N", "1000")); J <- as.integer(Sys.getenv("J", "30"))
B <- as.integer(Sys.getenv("B", "100")); SHAPE <- Sys.getenv("SHAPE", "smooth")
NRHO <- as.integer(Sys.getenv("NRHO", "12")); RHO_MIN <- as.numeric(Sys.getenv("RHO_MIN", "1"))
GAMMA_C <- as.numeric(strsplit(Sys.getenv("GAMMA_C", "1,1"), ",")[[1]])
BS <- 1:2; b1 <- c(0.5, -0.3); b2 <- c(-0.2, 0.4)
rg <- sort(exp(seq(log(.05), log(n / 2), length = NRHO)), decreasing = TRUE); rg <- rg[rg >= RHO_MIN]

make_shape <- function(shape, J) {                       # as in sim_tf3.R
  t <- 1:J; s <- J / 30; tt <- ceiling(t / s)
  if (shape == "smooth") { h1 <- 0.015 + 0.065 * exp(-((tt - 12)/6)^2); h2 <- 0.008 + 0.042 * tt / 30
  } else if (shape == "step") { h1 <- ifelse(tt <= 10, .03, ifelse(tt <= 20, .07, .03)); h2 <- ifelse(tt <= 15, .02, .04)
  } else if (shape == "abrupt") {
    h1 <- rep(.02, J); h1[tt %in% 9:12] <- c(.05, .09, .09, .05)[match(tt[tt %in% 9:12], 9:12)]
    h2 <- 0.01 + 0.07 * exp(-(tt - 1)/3)
  } else stop("unknown SHAPE")
  list(h1 = h1, h2 = h2)
}
sh <- make_shape(SHAPE, J); h1 <- sh$h1; h2 <- sh$h2
alp0 <- list(ph    = list(true_alpha(h1, h2, 1, "ph"),    true_alpha(h1, h2, 2, "ph")),
             logit = list(true_alpha(h1, h2, 1, "logit"), true_alpha(h1, h2, 2, "logit")))

PART <- file.path("sim_tf3_parts", sprintf("%s_n%d_J%d", SHAPE, n, J))
f <- sort(list.files(PART, "^rep_\\d+\\.rds$", full.names = TRUE))
if (!length(f)) stop("no saved replications in ", PART)
reps <- lapply(f, readRDS); DGPS <- names(reps[[1]])
pick <- function(dg, nm) lapply(reps, function(r) r[[dg]][[nm]])
arr  <- function(L) { d <- dim(L[[1]]); if (is.null(d)) d <- length(L[[1]])
                      a <- array(unlist(L), c(d, length(L)))
                      aperm(a, c(length(d) + 1, seq_along(d))) }
nms <- c("unpen", "lee", "tf", "rho", "df", "bic", "nonconv", "iter", "max_mu", "n_zero",
         "beta_u", "beta_lee", "beta_t")
res <- list()
for (dg in DGPS) res[[dg]] <- lapply(setNames(nm = nms), function(nm) arr(pick(dg, nm)))
OUT <- sprintf("sim_tf3_%s_n%d_J%d_partial.rds", SHAPE, n, J)
saveRDS(list(res = res, alp0 = alp0, h1 = h1, h2 = h2, n = n, J = J, B = length(f), shape = SHAPE,
             b1 = b1, b2 = b2, BS = BS, rho_grid = rg, gamma_c = GAMMA_C,
             reps = as.integer(sub("rep_(\\d+)\\.rds", "\\1", basename(f)))), OUT)
cat(sprintf("wrote %s with %d of %d replications\n", OUT, length(f), B))
