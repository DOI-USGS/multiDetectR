# 15_unseen_profile.R -- peer (paper chat) 2026-10-11: replace the stated prior on the number of unseen species U with a profile-likelihood
# interval, and run the design at its ends. a = 0.03; capsule and run multipliers held at the main fit's values.
#   lower end: maximise ll - lambda * U over the mixture weights for a grid of lambda > 0; U(lambda) falls as lambda grows; the end is
#              where the profile drops 1.92 log-likelihood units (95%, chi-square 1 df).
#   upper end: extra unseen mass is cheapest on the cell with the largest P0 (never-detected probability): each extra species there
#              costs (1 - P0max) / P0max units, so the 1.92 drop is reached at dU = 1.92 * P0max / (1 - P0max).
# Writes t_unseen_profile (ll by lambda) and fits_unseen_lo.rds / fits_unseen_hi.rds (in FITDIR) for 03 (env FITS).
source("common.R"); source("model.R")
fits <- readRDS(paste0(FITDIR, "fits.rds")); B <- readRDS(paste0(FITDIR, "build.rds")); a <- 0.03; IT <- 3000
MK <- if (nzchar(Sys.getenv("MARKER"))) Sys.getenv("MARKER") else MARKERS
LAM <- c(0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3)
prof <- list(); lo <- list(); hi <- list()
for (m in MK) { f0 <- fits[[paste(m, a)]]; L <- B$L[B$L$marker == m, ]; pr <- B$pres[[m]]
  fit <- function(...) fit_model(L, pr, a, run_mult = f0$run_mult, kap_capsule = f0$kap_capsule, iters = IT, ...)
  base <- fit(keep_post = FALSE); lg(m, "base ll", round(base$ll, 2), "U", round(base$unseen, 1), " (saved fit:", round(f0$ll, 2), round(f0$unseen, 1), ")")
  P <- bind_rows(data.frame(lambda = 0, ll = base$ll, U = base$unseen), bind_rows(lapply(LAM, function(l) { g <- fit(keep_post = FALSE, lambda = l)
    lg(m, "lambda", l, "ll", round(g$ll, 2), "U", round(g$unseen, 1)); data.frame(lambda = l, ll = g$ll, U = g$unseen) })))
  P$marker <- m; P$drop <- max(P$ll) - P$ll; prof[[m]] <- P
  o <- P[order(P$U), ]; j <- which(o$drop <= 1.92)[1]                     # smallest U still inside the interval, then interpolate in lambda
  lam_lo <- if (j > 1) approx(o$drop[(j - 1):j], o$lambda[(j - 1):j], xout = 1.92)$y else max(LAM)
  flo <- fit(lambda = lam_lo); lo[[paste(m, a)]] <- flo
  dU <- 1.92 * base$P0max / (1 - base$P0max); fhi <- fit(add_unseen = dU); hi[[paste(m, a)]] <- fhi
  lg(m, "upper end: P0max", signif(base$P0max, 6), "dU", round(dU, 1), "U", round(fhi$unseen, 1), "ll drop", round(max(P$ll) - fhi$ll, 2))
  lg(m, "lower end: lambda", round(lam_lo, 3), "U", round(flo$unseen, 1), "ll drop", round(max(P$ll) - flo$ll, 2)) }
SFX <- if (length(MK) == 1) paste0("_", MK) else ""
wcsv(bind_rows(prof), paste0("t_unseen_profile", SFX))
for (nm in c("lo", "hi")) { fp <- paste0(FITDIR, "fits_unseen_", nm, SFX, ".rds"); if (file.exists(fp)) stop("exists: ", fp); saveRDS(get(nm), fp) }
print(bind_rows(prof))
