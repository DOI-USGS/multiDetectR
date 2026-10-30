# 24_em_convergence.R -- peer (paper chat) 2026-10-11: what is the number of unseen species U after the EM has converged?
# The saved fits stop at 3,000 iterations, and U was still drifting up (15_unseen_profile.R). Here: the a = 0.03 fit (capsule and run
# multipliers as saved) run to 50,000 iterations, recording ll, U and the total at checkpoints, plus the observed-species design objective
# (expected observed species at the $5,000 12S+COI design, V4 3/2/0 f1 k4|12) so the effect on designs is visible. Env FITS = fits_template.rds
# repeats it for the template fit (written with a _template suffix).
source("common.R"); source("model.R")
OUT <- FITDIR; B <- readRDS(paste0(FITDIR, "build.rds")); FN <- Sys.getenv("FITS", "fits.rds"); fits <- readRDS(paste0(FITDIR, FN)); a <- 0.03
SFX <- if (FN == "fits.rds") "" else paste0("_", sub("^fits_(.*)\\.rds$", "\\1", FN))
TR <- c(1, 100, 300, 1000, 1500, 3000, 5000, 10000, 20000, 30000, 50000); IT <- max(TR) + 1
res <- list(); des <- list()
for (m in MARKERS) { f0 <- fits[[paste(m, a)]]; L <- B$L[B$L$marker == m, ]; pr <- B$pres[[m]]; t0 <- Sys.time()
  ff <- fit_model(L, pr, a, run_mult = f0$run_mult, kap_capsule = f0$kap_capsule, iters = IT, tmpl = f0$tmpl, trace_at = TR + 1)
  lg(m, "done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min; final U", round(ff$unseen, 1))
  tr <- ff$trace; tr$marker <- m; tr$Nobs <- ff$Nobs; tr$dll_from_previous <- c(NA, diff(tr$ll)); res[[m]] <- tr
  for (fit in list(list("saved (3,000 iterations)", f0), list("converged run", ff))) { ce <- fit[[2]]$cells
    P <- det_pc(ce, a, V = 4, b = c(GP = 3, CJ = 2, JB = 0), f = 1, k = if (m == "12S") 4 else 12, kap_capsule = fit[[2]]$kap_capsule)
    des[[paste(m, fit[[1]])]] <- data.frame(marker = m, fit = fit[[1]], iterations = if (fit[[1]] == "converged run") IT - 1 else 3000, U = fit[[2]]$unseen,
      S_observed_at_5000_design = sum(colSums(fit[[2]]$post) * P), S_all_at_5000_design = sum(ce$nu * P)) } }
R <- bind_rows(res); D <- bind_rows(des)
wcsv(R, paste0("t_em_convergence", SFX)); wcsv(D, paste0("t_em_convergence_design", SFX))
print(R %>% mutate(across(c(ll, U, Ntot, dll_from_previous), ~ signif(.x, 7))) %>% as.data.frame()); print(D)
