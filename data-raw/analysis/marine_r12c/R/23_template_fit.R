# 23_template_fit.R -- peer (paper chat) 2026-10-11 evening: tube-level template model. 18 found the reactions of one filter strongly dependent
# given the species type (j of n singles U-shaped). Here template amount varies among filters: given the filter carries the species, a class G
# (mean 1) multiplies every reaction's hazard -log(1 - p). Main form (peer, 2026-10-11): G ~ Gamma(alpha, rate alpha), so a filter's miss
# probability over a pooled library is (1 + beta m_eff(k))^(-alpha), beta = -log(1 - p) / alpha, in closed form (model.R tmpl_gamma).
# Alternative: two classes, a share w_low of filters at rho x the hazard of the rest (tmpl_two_class).
# a = 0.03 throughout. Steps per marker:
#   (0) checks: tmpl = NULL reproduces the saved fit; one class G = 1 through the new code gives the same likelihood;
#   (1) profile alpha and the two-class grid with capsule and run multipliers held at the saved fit;
#   (2) COI without the May 2024 run (JVB3506), whose tap-water blanks are as rich as samples (22_single_reaction_detections.R);
#   (3) at the best alpha, re-choose capsule and run multipliers (one sweep on 02's grids), re-profile alpha, final fit (3000 iterations);
#   (4) j-of-n check on the single-reaction runs under the final fit (compare 18).
# Writes fits_template.rds (FITDIR; names as fits.rds) for 03 (env FITS), and t_template_* tables.
source("common.R"); source("model.R")
OUT <- FITDIR; B <- readRDS(paste0(FITDIR, "build.rds")); fits <- readRDS(paste0(FITDIR, "fits.rds")); a <- 0.03; PI <- 300; IT <- 1500; ITF <- 3000
KAP <- c(Inf, 8, 4, 2, 1, 0.5, 0.25, 0.1)          # alpha (gamma shape); Inf = no heterogeneity
TWO <- expand.grid(w_low = c(0.3, 0.5, 0.7), rho = c(0.02, 0.1, 0.3))
KC <- c(.125, .25, .5, 1, 2); CG <- c(.125, .25, .5, 1, 2, 4, 8)
prof <- list(); summ <- list(); runm <- list(); jchk <- list(); out <- list()
for (m in MARKERS) { f0 <- fits[[paste(m, a)]]; L <- B$L[B$L$marker == m, ]; pr <- B$pres[[m]]
  fit <- function(..., L. = L, pr. = pr, cv = f0$run_mult, kc = f0$kap_capsule, iters = IT) fit_model(L., pr., a, run_mult = cv, kap_capsule = kc, iters = iters, ...)
  t0 <- Sys.time(); c0 <- fit(iters = ITF, keep_post = FALSE); lg(m, "check: tmpl NULL ll", round(c0$ll, 4), "saved", round(f0$ll, 4), "(", round(difftime(Sys.time(), t0, units = "secs")), "s)")
  stopifnot(abs(c0$ll - f0$ll) < 1e-6)
  b1 <- fit(keep_post = FALSE)$ll; t0 <- Sys.time(); c1 <- fit(keep_post = FALSE, tmpl = list(G = 1, w = 1))$ll
  lg(m, "check: one class through new code", round(c1 - b1, 8), "(", round(difftime(Sys.time(), t0, units = "secs")), "s)"); stopifnot(abs(c1 - b1) < 1e-6)
  P1 <- bind_rows(lapply(KAP, function(k) { t0 <- Sys.time(); g <- fit(keep_post = FALSE, tmpl = tmpl_gamma(k))
    lg(m, "held: alpha", k, "ll", round(g$ll, 2), "U", round(g$unseen, 1), "(", round(difftime(Sys.time(), t0, units = "secs")), "s)")
    data.frame(marker = m, step = "alpha, multipliers held", data = "all", alpha = k, w_low = NA, rho = NA, ll = g$ll, U = g$unseen) }))
  P2 <- bind_rows(lapply(seq_len(nrow(TWO)), function(i) { g <- fit(keep_post = FALSE, tmpl = tmpl_two_class(TWO$w_low[i], TWO$rho[i]))
    lg(m, "held: two-class", TWO$w_low[i], TWO$rho[i], "ll", round(g$ll, 2))
    data.frame(marker = m, step = "two-class, multipliers held", data = "all", alpha = NA, w_low = TWO$w_low[i], rho = TWO$rho[i], ll = g$ll, U = g$unseen) }))
  P3 <- NULL
  if (m == "COI") { L2 <- L[L$batch != "JVB3506", ]; pr2 <- pr[pr$library_id %in% L2$library_id, ]
    P3 <- bind_rows(lapply(KAP, function(k) { g <- fit(L. = L2, pr. = pr2, keep_post = FALSE, tmpl = tmpl_gamma(k)); lg(m, "without JVB3506: alpha", k, "ll", round(g$ll, 2))
      data.frame(marker = m, step = "alpha, multipliers held", data = "without JVB3506", alpha = k, w_low = NA, rho = NA, ll = g$ll, U = g$unseen) })) }
  kb <- P1$alpha[which.max(P1$ll)]; tb <- tmpl_gamma(kb); cv <- f0$run_mult; kc <- f0$kap_capsule
  if (is.finite(kb)) {                                                                              # (3) re-choose multipliers at kb
    llf <- function(cv, kc) fit(cv = cv, kc = kc, iters = PI, keep_post = FALSE, tmpl = tb)$ll
    l <- sapply(KC, function(k) llf(cv, k)); kc <- KC[which.max(l)]
    for (b in names(cv)) { l <- sapply(CG, function(c) { z <- cv; z[b] <- c; llf(z, kc) }); cv[b] <- CG[which.max(l)] }
    lg(m, "re-chosen at alpha", kb, ": capsule", kc, " runs", paste(names(cv), cv, sep = "=", collapse = " ")) }
  P4 <- bind_rows(lapply(KAP, function(k) { g <- fit(cv = cv, kc = kc, keep_post = FALSE, tmpl = tmpl_gamma(k)); lg(m, "re-chosen: alpha", k, "ll", round(g$ll, 2))
    data.frame(marker = m, step = "alpha, multipliers re-chosen", data = "all", alpha = k, w_low = NA, rho = NA, ll = g$ll, U = g$unseen) }))
  kf <- P4$alpha[which.max(P4$ll)]
  if (is.finite(kf)) { lr <- optimize(function(lk) fit(cv = cv, kc = kc, keep_post = FALSE, tmpl = tmpl_gamma(exp(lk)))$ll, log(kf) + c(-0.7, 0.7), maximum = TRUE, tol = 0.02)
    kf <- exp(lr$maximum); lg(m, "refined alpha", round(kf, 3), "ll", round(lr$objective, 2)) }
  ff <- fit(cv = cv, kc = kc, iters = ITF, tmpl = tmpl_gamma(kf)); out[[paste(m, a)]] <- ff
  lg(m, "final: alpha", kf, "ll", round(ff$ll, 2), "vs saved", round(f0$ll, 2), " U", round(ff$unseen, 1), "vs", round(f0$unseen, 1))
  prof[[m]] <- bind_rows(P1, P2, P3, P4)
  summ[[m]] <- data.frame(marker = m, alpha = kf, ll_template = ff$ll, ll_saved = f0$ll, dll = ff$ll - f0$ll, extra_params = 1,
    p_value_chisq1 = pchisq(2 * (ff$ll - f0$ll), 1, lower.tail = FALSE), best_two_class = paste(P2$w_low[which.max(P2$ll)], P2$rho[which.max(P2$ll)]),
    dll_two_class_held = max(P2$ll) - P1$ll[is.infinite(P1$alpha)], dll_gamma_held = max(P1$ll) - P1$ll[is.infinite(P1$alpha)],
    species_est_total = ff$Ntot, unseen = ff$unseen, unseen_saved = f0$unseen, kappa_capsule = kc)
  runm[[m]] <- data.frame(marker = m, batch = c(REF_BATCH, names(cv)), run_multiplier_template = c(1, cv), run_multiplier_saved = c(1, f0$run_mult[names(cv)]))
  # (4) j of n on the single-reaction runs (one library = one reaction; M = c_b), expected from the final fit's species posteriors
  for (bt in c("JVB3506", "JVB3105")) { Ls <- L %>% filter(batch == bt, k == 1); if (!nrow(Ls)) next
    ce <- ff$cells; TGc <- ff$tmpl; if (!is.null(TGc) && !identical(TGc$type, "gamma")) next; cb <- if (bt %in% names(cv)) cv[[bt]] else 1
    det <- Ls %>% inner_join(pr, by = "library_id", relationship = "many-to-many") %>% count(Barcode, taxon, name = "j") %>% left_join(Ls %>% count(Barcode, name = "n"), by = "Barcode")
    ex <- bind_rows(lapply(seq_len(nrow(det)), function(i) { n <- det$n[i]; w <- ff$post[det$taxon[i], ]
      up <- sort(unique(ce$p)); hp <- -cb * log(1 - up)
      tab <- if (is.null(TGc)) t(sapply(up, function(p) dbinom(1:n, n, 1 - (1 - p)^cb))) else t(sapply(hp, function(h) sapply(1:n, function(x) choose(n, x) * gamma_pattern(rep(h, x), (n - x) * h, TGc$alpha))))
      pj <- tab[match(ce$p, up), , drop = FALSE]; pj <- pj / rowSums(pj); data.frame(n = n, j = 1:n, expected = as.vector(crossprod(w, pj))) }))
    G <- ex %>% group_by(n, j) %>% summarise(expected = sum(expected), .groups = "drop") %>% left_join(det %>% count(n, j, name = "observed"), by = c("n", "j")) %>%
      mutate(observed = coalesce(observed, 0L))
    jchk[[paste(m, bt)]] <- G %>% summarise(marker = m, batch = bt, alpha = kf, pairs = sum(observed), obs_j1 = sum(observed[j == 1]), exp_j1 = sum(expected[j == 1]),
      obs_jn = sum(observed[j == n]), exp_jn = sum(expected[j == n])) } }
fp <- paste0(FITDIR, "fits_template.rds"); if (file.exists(fp)) stop("exists: ", fp); saveRDS(out, fp)
wcsv(bind_rows(prof), "t_template_profile"); wcsv(bind_rows(summ), "t_template_summary"); wcsv(bind_rows(runm), "t_template_run_multipliers"); wcsv(bind_rows(jchk), "t_template_j_check")
print(bind_rows(summ)); print(bind_rows(jchk)); print(bind_rows(prof) %>% group_by(marker, step, data) %>% mutate(dll = ll - max(ll)) %>% as.data.frame())
