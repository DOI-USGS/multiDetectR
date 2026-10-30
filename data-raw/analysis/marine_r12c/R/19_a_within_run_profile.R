# 19_a_within_run_profile.R -- peer (paper chat) 2026-10-11, companion to 18: where does the full-data preference for a = 0.01 come from?
# In the full data every run has one PCR count except 12S JVB6097 (k 6 and 12), 12S JVB3735 (6 and 18) and COI JVB3735 (6 and 12); the
# other runs differ in k only BETWEEN runs, and each run has its own multiplier c_b. 02 compared a = 0.01 / 0.03 / 0.05 with c_b and the
# capsule multiplier re-chosen on doubling grids, so part of the 8-33 unit differences can be grid snapping of the overall scale.
# Here a is profiled with capsule and run multipliers held at the a = 0.03 fit, in two versions:
#   "within-run": M = c_b * kref_b * m_eff(k, a) / m_eff(kref_b, a), kref_b = the run's smallest k. A run's scale is the same at every a,
#                 so a is informed only by the k contrasts inside the three mixed runs.
#   "raw":        M = c_b * m_eff(k, a), as in 02 but with c_b held (so a also moves every run's scale against the p grid).
source("common.R"); source("model.R")
OUT <- FITDIR; B <- readRDS(paste0(FITDIR, "build.rds")); fits <- readRDS(paste0(FITDIR, "fits.rds")); IT <- 1500
AG <- c(0.001, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2)
P <- list()
for (m in MARKERS) { f0 <- fits[[paste(m, 0.03)]]; L <- B$L[B$L$marker == m, ]; pr <- B$pres[[m]]
  kref <- ave(L$k, L$batch, FUN = min)
  lg(m, "runs with more than one k:", paste(names(which(tapply(L$k, L$batch, function(z) length(unique(z))) > 1)), collapse = " "))
  for (v in c("within-run", "raw")) for (a in AG) {
    w <- if (v == "within-run") kref * m_eff(L$k, a) / m_eff(kref, a) else m_eff(L$k, a)
    f <- fit_model(L, pr, 0.03, run_mult = f0$run_mult, kap_capsule = f0$kap_capsule, iters = IT, keep_post = FALSE, M_lib = w)
    lg(m, v, "a", a, "ll", round(f$ll, 3)); P[[length(P) + 1]] <- data.frame(marker = m, version = v, a = a, ll = f$ll, U = f$unseen) } }
P <- bind_rows(P) %>% group_by(marker, version) %>% mutate(drop = max(ll) - ll) %>% ungroup()
S <- P %>% group_by(marker, version) %>% group_modify(function(d, k) { d <- d[order(d$a), ]; i <- which.max(d$ll)
  end <- function(ix) { z <- d[ix, ]; j <- which(z$drop > 1.92)[1]; if (is.na(j)) return(NA_real_); approx(z$drop[(j - 1):j], z$a[(j - 1):j], xout = 1.92)$y }
  data.frame(a_best_on_grid = d$a[i], lower = if (i > 1) end(i:1) else NA_real_, upper = end(i:nrow(d)), drop_at_0.001 = d$drop[1], drop_at_0.03 = d$drop[d$a == 0.03]) }) %>% ungroup()
wcsv(P, "t_a_within_run_profile"); wcsv(S, "t_a_within_run_summary")
print(as.data.frame(P %>% mutate(across(where(is.numeric), ~ round(.x, 3))))); print(as.data.frame(S))
