# 02_fit.R -- per marker and extract fraction a: capsule multiplier, run multipliers (vs REF_BATCH), final fit, level tests, calibration.
# One R process at a time (peak ~1.5 GB for COI).
# PORTED OPTIMIZATION (2026-10-26): the capsule (KC) and per-batch run-multiplier (CG) grid
# searches below are the hot loop -- ~75 fit_model() calls per sweep (10 non-ref batches x 7 CG
# values + 5 KC values), x2 sweeps x2 markers, which is the bulk of the ~3 h single-core runtime.
# Each grid is a set of INDEPENDENT llf() evaluations followed by which.max(), so we evaluate the
# grid with gmap() = parallel::mclapply over NCORE fork workers instead of sapply. This changes
# NOTHING numerically: same grid, same values, same which.max selection, same coordinate-ascent
# order (the batch loop stays sequential so cv[b] still feeds the next batch). NCORE defaults to 1
# => identical to the original serial behaviour unless the driver sets it. Reuses the exact
# NCORE/mclapply convention already in 29_rule_visit_bands.R.
source("common.R"); source("model.R")
B <- readRDS(paste0(OUT, "build.rds")); PI <- 300; IT <- 3000
.NCORE <- as.integer(Sys.getenv("NCORE", "1")); if (is.na(.NCORE) || .NCORE < 1L) .NCORE <- 1L
# gmap(): map fn over a grid, serially when NCORE==1 (byte-identical to sapply) else over forks.
# mc.preschedule=FALSE so uneven per-value times (EM convergence varies) balance across workers;
# the result order matches the input grid, so which.max picks the same element either way.
gmap <- if (.NCORE > 1L && requireNamespace("parallel", quietly = TRUE)) {
  function(grid, fn) unlist(parallel::mclapply(grid, fn, mc.cores = .NCORE, mc.preschedule = FALSE))
} else function(grid, fn) sapply(grid, fn)
# env A: fit one extract fraction only (outputs suffixed _a<A>; 2026-10-10 Jonah Ventures: an extract supports ~12 tests -> a = 1/12 if a test is one PCR)
AS <- if (nzchar(Sys.getenv("A"))) as.numeric(Sys.getenv("A")) else A_GRID; SFX <- if (nzchar(Sys.getenv("A"))) paste0("_a", Sys.getenv("A")) else ""
KC <- c(.125, .25, .5, 1, 2); CG <- c(.125, .25, .5, 1, 2, 4, 8)
pred_unit <- function(f, u) { tb <- f$tubes[f$tubes$visit_id == u$visit_id & f$tubes$site == u$site, ]; ce <- f$cells; lm <- 0
  for (bg in unique(tb$bagkey)) { lf <- 0
    for (j in which(tb$bagkey == bg)) { M <- sum(f$L$M[f$L$Barcode == tb$Barcode[j]]); kk <- if (tb$ftype[j] == "capsule") f$kap_capsule else 1
      lf <- lf + log(1 - odds_mult(ce$r, kk) * (1 - (1 - ce$p)^M)) }
    lm <- lm + log(1 - ce$q * (1 - exp(lf))) }
  sum(ce$nu * ce[[paste0("phi_", u$site)]] * (1 - exp(lm))) }
fits <- list(); summ <- list(); runm <- list(); levl <- list(); cal <- list()
for (m in MARKERS) {
  L <- B$L[B$L$marker == m, ]; pr <- B$pres[[m]]; n <- n_distinct(pr$taxon)
  bt <- setdiff(sort(unique(L$batch)), REF_BATCH); cv <- setNames(rep(1, length(bt)), bt); kc <- 1; nu <- NULL
  # cold starts only: the EM update is multiplicative, so a warm start locks out cells that reached ~0 (2026-10-10: warm starts gave a full model WORSE than its nested r=1 model)
  llf <- function(cv, kc, a) fit_model(L, pr, a, run_mult = cv, kap_capsule = kc, iters = PI, keep_post = FALSE)$ll
  for (a in AS[order(abs(AS - 0.03))]) {                           # start at a = 0.03, then the bracket ends
    sweeps <- if (a == 0.03 || length(AS) == 1) 2 else 1
    for (sw in seq_len(sweeps)) {
      l <- gmap(KC, function(k) llf(cv, k, a)); kc <- KC[which.max(l)]
      for (b in bt) { l <- gmap(CG, function(c) { z <- cv; z[b] <- c; llf(z, kc, a) }); cv[b] <- CG[which.max(l)] }
      lg(m, "a", a, "sweep", sw, ": capsule", kc, " runs", paste(names(cv), cv, sep = "=", collapse = " ")) }
    f <- fit_model(L, pr, a, run_mult = cv, kap_capsule = kc, iters = IT)
    f0 <- fit_model(L, pr, a, run_mult = NULL, kap_capsule = kc, iters = IT, keep_post = FALSE)             # n_pcr only, no run multipliers
    fq <- fit_model(L, pr, a, run_mult = cv, kap_capsule = kc, iters = PI * 2, qfix = .999, keep_post = FALSE)
    fr <- fit_model(L, pr, a, run_mult = cv, kap_capsule = kc, iters = PI * 2, rfix = .999, keep_post = FALSE)
    fits[[paste(m, a)]] <- f
    summ[[paste(m, a)]] <- data.frame(marker = m, a = a, species_observed = n, species_est_total = f$Ntot, kappa_capsule = kc, ll = f$ll, ll_no_run_multipliers = f0$ll,
                                      dll_run_multipliers = f$ll - f0$ll, n_run_params = length(cv), BIC_supports = (f$ll - f0$ll) > 0.5 * log(n) * length(cv))
    runm[[paste(m, a)]] <- data.frame(marker = m, a = a, batch = c(REF_BATCH, names(cv)), run_multiplier = c(1, cv)) %>%
      mutate(platform = ifelse(batch %in% names(PLATFORM), "illumina", "nanopore"))
    levl[[paste(m, a)]] <- data.frame(marker = m, a = a, model = c("no bag level (q=1)", "no filter level (r=1)"), dll_vs_full = c(fq$ll, fr$ll) - f$ll)
    cal[[paste(m, a)]] <- bind_rows(lapply(seq_len(nrow(f$units)), function(i) { u <- f$units[i, ]
      data.frame(marker = m, a = a, site = u$site, visit_id = u$visit_id,
                 observed = n_distinct(f$pres$taxon[f$pres$library_id %in% f$L$library_id[f$L$visit_id == u$visit_id & f$L$site == u$site]]), expected = pred_unit(f, u)) }))
    lg(m, "a", a, "final ll", round(f$ll, 1), " vs n_pcr only", round(f0$ll, 1), " Ntot", round(f$Ntot, 1))
  }
}
# PORTED ROBUSTNESS (2026-10-26): save the fits object FIRST. It is the only
# expensive-to-recreate artifact (the whole EM fit); the four t_*.csv tables are
# cheap summaries built from `fits` in memory. Writing fits.rds before the CSVs
# means a glitch while assembling a CSV cannot throw away hours of fitting -- the
# original order (CSVs then saveRDS) did exactly that once. Numerically identical:
# same object, same tables, just reordered writes.
fp <- paste0(OUT, "fits", SFX, ".rds"); if (file.exists(fp)) stop("exists: ", fp); saveRDS(fits, fp)
wcsv(bind_rows(summ), paste0("t_fit_summary", SFX)); wcsv(bind_rows(runm), paste0("t_run_multipliers", SFX)); wcsv(bind_rows(levl), paste0("t_level_tests", SFX)); wcsv(bind_rows(cal), paste0("t_calibration", SFX))
print(bind_rows(summ)); print(bind_rows(levl)); print(bind_rows(cal) %>% group_by(marker, a, site) %>% summarise(units = n(), obs = mean(observed), exp = mean(expected), .groups = "drop"), n = 50)
