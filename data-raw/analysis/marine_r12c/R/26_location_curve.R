# 26_location_curve.R -- peer (paper chat) 2026-10-11: expected target species against the number of coastline locations (1-12) on a 1,000 km
# stretch straddling Pt. Conception, from the fitted distance decay (20_turnover.R) with its step across the point; exact expectations, no simulation.
# Positions (km along the coast, Pt. Conception at 0, north negative) in the spreading order: centre, ends, then fill in by bisection:
#   0, -500, 500, -250, 250, -375, -125, 125, 375, -437.5, -312.5, -187.5. The centre location counts as south of the point.
# Presence model, per species type c (fitted mixture, a = 0.03): presence at a location on a visit has marginal psi = phi_COAST(c). Along the
# coast it is a two-state Markov chain between adjacent locations: P(present at B | present at A) = pi(d_AB) = [pi_inf + (pi0 - pi_inf)
# exp(-d/Lc)] * step^cross (capped at 1), and P(present at B | absent at A) = psi (1 - pi) / (1 - psi), which keeps the marginal at psi.
# Detection given presence: one visit, B_LOC bags, 1 disc filter, k PCRs (12S 4, COI 12: the inner design of the Pt. Conception rule design),
# through the fit's own levels (incl. the template classes when FITS = fits_template.rds). P(never detected) by the forward recursion: exact
# under this chain. Check: with pi = psi (no spatial structure) it reduces to 1 - (1 - psi D)^n.
# Marginal species per dollar of the n-th location (cost: travel_location + B_LOC x (sample + filter + library + k pcr)) is compared with that of
# ANOTHER VISIT: repeating the visit at all n locations, visits independent given the type (as phi is per visit in the fit), so the gain is
# S(n, 2 visits) - S(n, 1 visit) for n x the location cost.
# Decay parameters: the 20_turnover.R fit (coastal distance, Event3) and its species bootstrap (200 draws, same seed) for the range.
Sys.setenv(TURNOVER_DEFS_ONLY = "1"); source("20_turnover.R")
FN <- Sys.getenv("FITS", "fits.rds"); DF <- readRDS(paste0(FITDIR, FN)); SFX <- if (FN == "fits.rds") "" else paste0("_", sub("^fits_(.*)\\.rds$", "\\1", FN))
C <- read_costs(); B_LOC <- 3; KM <- c(`12S` = 4, COI = 12)
POS <- c(0, -500, 500, -250, 250, -375, -125, 125, 375, -437.5, -312.5, -187.5)
curve <- function(th, f, k) { ce <- f$cells; w <- colSums(f$post); psi <- ce$phi_COAST
  D <- 1 - (1 - ce$q * filt_det(ce, a, k, "disc", f$kap_capsule))^B_LOC
  out <- list()
  for (n in seq_along(POS)) { x <- sort(POS[seq_len(n)]); a0 <- 1 - psi; a1 <- psi * (1 - D)
    if (n > 1) for (i in 2:n) { d <- x[i] - x[i - 1]; cr <- as.numeric((x[i - 1] < 0) != (x[i] < 0))
      p11 <- pmin(1, pi_fn(th, d, cr)); p01 <- ifelse(psi < 1, pmin(1, psi * (1 - p11) / pmax(1 - psi, 1e-12)), 0)
      n0 <- a0 * (1 - p01) + a1 * (1 - p11); n1 <- (a0 * p01 + a1 * p11) * (1 - D); a0 <- n0; a1 <- n1 }
    pnd <- a0 + a1; out[[n]] <- data.frame(n_locations = n, S = sum(w * (1 - pnd)), S_two_visits = sum(w * (1 - pnd^2)), S_all = sum(ce$nu * (1 - pnd)),
      S_independent = sum(w * (1 - (1 - psi * D)^n))) }
  bind_rows(out) }
# check: pi = psi reproduces independence (pi_fn returns pi; set pi0 = pi_inf = psi is type-specific, so check with a per-type override)
chk <- local({ f <- DF[[paste("12S", a)]]; ce <- f$cells; psi <- ce$phi_COAST; D <- 1 - (1 - ce$q * filt_det(ce, a, 4, "disc", f$kap_capsule))^B_LOC
  a0 <- 1 - psi; a1 <- psi * (1 - D); for (i in 2:5) { p11 <- psi; p01 <- psi; n0 <- a0 * (1 - p01) + a1 * (1 - p11); n1 <- (a0 * p01 + a1 * p11) * (1 - D); a0 <- n0; a1 <- n1 }
  max(abs((a0 + a1) - (1 - psi * D)^5)) }); lg("independence check, max abs difference:", signif(chk, 3)); stopifnot(chk < 1e-12)
RES <- list(); CROSS <- list(); PAR <- list()
for (m in MARKERS) { f <- DF[[paste(m, a)]]; k <- KM[[m]]; cl <- C[["travel_location"]] + B_LOC * (C[["sample"]] + C[["filter_disc"]] + C[["library"]] + k * C[["pcr"]])
  R <- build_pairs(m, PC_EVENT)$R %>% mutate(taxon = paste(marker, taxon)); full <- fit1(R, "dist_coast")
  set.seed(20261011); sp <- unique(R$taxon); idx <- split(seq_len(nrow(R)), R$taxon)
  BT <- t(sapply(seq_len(NBOOT), function(b) { s <- sample(sp, replace = TRUE); fit1(R[unlist(idx[s], use.names = FALSE), ], "dist_coast", start = full$th)$th }))
  PAR[[m]] <- data.frame(marker = m, draw = c(0, seq_len(NBOOT)), rbind(full$th, BT)) %>% setNames(c("marker", "draw", "logit_pi_inf", "logit_pi0", "log_Lc", "log_step")) %>%
    mutate(Lc_km = exp(log_Lc), step = exp(log_step))
  cv <- function(th, draw) { z <- curve(th, f, k); z$marker <- m; z$draw <- draw; z$cost_per_location <- cl
    z$marginal_per_dollar_location <- c(z$S[1], diff(z$S)) / cl
    z$marginal_per_dollar_another_visit <- (z$S_two_visits - z$S) / (z$n_locations * cl); z }
  Z0 <- cv(full$th, 0); ZB <- bind_rows(lapply(seq_len(NBOOT), function(b) cv(BT[b, ], b)))
  RES[[m]] <- bind_rows(Z0, ZB)
  # the location count at which the NEXT location is worth less per dollar than another visit to the current set
  nstar <- function(z) { i <- which(z$marginal_per_dollar_location[-1] < z$marginal_per_dollar_another_visit[-nrow(z)])[1]; if (is.na(i)) NA else i }
  CROSS[[m]] <- data.frame(marker = m, n_star_fit = nstar(Z0), n_star_boot = sapply(split(ZB, ZB$draw), nstar)) }
RS <- bind_rows(RES); CR <- bind_rows(CROSS)
SUM <- RS %>% group_by(marker, n_locations) %>% summarise(position_added_km = POS[first(n_locations)], S_fit = S[draw == 0],
  S_boot_lo = quantile(S[draw > 0], .025), S_boot_hi = quantile(S[draw > 0], .975), S_if_independent = S_independent[draw == 0],
  S_two_visits_fit = S_two_visits[draw == 0], marginal_location_fit = marginal_per_dollar_location[draw == 0] * 1000,
  marginal_another_visit_fit = marginal_per_dollar_another_visit[draw == 0] * 1000, .groups = "drop")
NS <- CR %>% group_by(marker) %>% summarise(n_star_fit = first(n_star_fit), n_star_boot_lo = quantile(n_star_boot, .025, na.rm = TRUE, type = 1),
  n_star_boot_hi = quantile(n_star_boot, .975, na.rm = TRUE, type = 1), boot_never = sum(is.na(n_star_boot)), .groups = "drop") %>%
  left_join(bind_rows(PAR) %>% group_by(marker) %>% summarise(Lc_fit = Lc_km[draw == 0], Lc_boot_lo = quantile(Lc_km[draw > 0], .025), Lc_boot_hi = quantile(Lc_km[draw > 0], .975),
    step_fit = step[draw == 0], .groups = "drop"), by = "marker")
wcsv(SUM, paste0("t_location_curve", SFX)); wcsv(NS, paste0("t_location_curve_crossing", SFX)); wcsv(bind_rows(PAR), paste0("t_location_curve_decay_draws", SFX))
print(SUM %>% mutate(across(where(is.double), ~ round(.x, 2))) %>% as.data.frame()); print(as.data.frame(NS))
