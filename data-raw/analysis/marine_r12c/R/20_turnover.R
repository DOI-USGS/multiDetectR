# 20_turnover.R -- peer (paper chat) 2026-10-11: turnover between locations against coastal distance, corrected for detection with the
# fitted inner levels (q, r, p; a = 0.03 fit), and whether there is a step across Pt. Conception.
# Location table: every location's habitat class, coordinates, samples (bags), filters and visits.
# Units: the 21 coastline locations (one visit each) and the three Pt. Conception sites on ONE event (default Event3, 2024-05-29, the
# event nearest the summer-2024 expedition; env PC_EVENT; every event is also run as a sensitivity).
# Detection: d[u, s] = P(unit u's actual effort detects species s | s present at u) = sum over types c of post[s, c] * D_u(c), with
#   D_u(c) = 1 - prod_bags (1 - q (1 - prod_filters (1 - r' (1 - (1 - p)^M_filter)))), M = run multiplier x m_eff(k, a) summed over libraries.
# Model (conditional, per marker): for an ordered pair of units (A, B) and a species observed at A (so present there; no false positives),
#   P(observed at B) = d[B, s] * pi(A, B),   pi = [pi_inf + (pi_0 - pi_inf) exp(-dist / Lc)] * exp(beta * cross),  capped at 0.999,
#   pi = P(present at B | present at A); cross = the pair lies on opposite sides of Pt. Conception. Unseen species are not included.
# Intervals: composite-likelihood profile for Lc (too narrow: pairs share species and units), species bootstrap (200), and
# leave-one-location-out. Distances: along the Natural Earth 1:10m mainland coastline (islands: great circle to everything); great circle
# alone as a sensitivity.
source("common.R"); source("model.R")
suppressPackageStartupMessages(library(sf)); sf_use_s2(FALSE)
OUT <- FITDIR; B <- readRDS(paste0(FITDIR, "build.rds")); fits <- readRDS(paste0(FITDIR, "fits.rds")); a <- 0.03; set.seed(20261011)
PC_EVENT <- Sys.getenv("PC_EVENT", "Event3"); NBOOT <- 200
# ---------- location table ----------
S0 <- arc("samples")
LOC <- S0 %>% group_by(location = Site, area = Area) %>%
  summarise(habitat_class = paste(unique(habitat), collapse = "/"), latitude = mean(Latitude), longitude = mean(Longitude),
            visits = n_distinct(visit_id), samples_bags = n_distinct(paste(visit_id, bag)), filters_tubes = n(),
            filter_types = paste(sort(unique(filter_type[nzchar(filter_type)])), collapse = "/"), first_date = min(Date), last_date = max(Date), .groups = "drop")
for (m in MARKERS) LOC[[paste0("libraries_", m)]] <- as.vector(table(factor(B$L$Site[B$L$marker == m], levels = LOC$location)))
# ---------- distances ----------
cl <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf")
bb <- st_as_sfc(st_bbox(c(xmin = -123.5, xmax = -112, ymin = 26, ymax = 38.5), crs = 4326))
x <- st_cast(st_cast(suppressWarnings(st_intersection(st_geometry(cl), bb)), "MULTILINESTRING"), "LINESTRING")
x <- st_transform(x, 3310); main <- x[which.max(as.numeric(st_length(x)))]            # the mainland Pacific coast, 38.5 N to 26 N
V <- st_coordinates(main)[, 1:2]; seglen <- sqrt(rowSums(diff(V)^2)); cum <- c(0, cumsum(seglen))
snap <- function(pt) { A <- V[-nrow(V), ]; D <- diff(V); t <- pmin(1, pmax(0, ((pt[1] - A[, 1]) * D[, 1] + (pt[2] - A[, 2]) * D[, 2]) / pmax(seglen^2, 1e-9)))
  px <- A[, 1] + t * D[, 1]; py <- A[, 2] + t * D[, 2]; dd <- sqrt((pt[1] - px)^2 + (pt[2] - py)^2); i <- which.min(dd)
  c(pos_km = (cum[i] + t[i] * seglen[i]) / 1000, off_km = dd[i] / 1000) }
P3310 <- st_coordinates(st_transform(st_as_sf(LOC, coords = c("longitude", "latitude"), crs = 4326), 3310))
sn <- t(apply(P3310, 1, snap)); LOC$coast_pos_km <- sn[, "pos_km"]; LOC$off_coast_km <- sn[, "off_km"]; LOC$island <- LOC$off_coast_km > 5
ptc <- snap(st_coordinates(st_transform(st_sfc(st_point(c(-120.4713, 34.4486)), crs = 4326), 3310))[1, ])   # Point Conception light
LOC$side <- ifelse(!LOC$island & LOC$coast_pos_km < ptc["pos_km"], "north of Pt. Conception", "south/east (Bight, islands, Baja)")
GC <- geosphere::distm(cbind(LOC$longitude, LOC$latitude), fun = geosphere::distHaversine) / 1000; dimnames(GC) <- list(LOC$location, LOC$location)
CD <- GC; for (i in seq_len(nrow(LOC))) for (j in seq_len(nrow(LOC))) if (!LOC$island[i] && !LOC$island[j]) CD[i, j] <- abs(LOC$coast_pos_km[i] - LOC$coast_pos_km[j])
LOC$coast_pos_km[LOC$island] <- NA
DEFS_ONLY <- nzchar(Sys.getenv("TURNOVER_DEFS_ONLY"))                                   # 21 sources this file for its definitions only
if (!DEFS_ONLY) wcsv(LOC %>% mutate(across(c(latitude, longitude), ~ round(.x, 5)), across(c(coast_pos_km, off_coast_km), ~ round(.x, 1))), "t_turnover_locations")
# ---------- detection per unit ----------
unit_det <- function(f, Lu) { ce <- f$cells; lm <- 0
  for (bg in unique(Lu$bagkey)) { lf <- 0
    for (bc in unique(Lu$Barcode[Lu$bagkey == bg])) { M <- sum(Lu$M[Lu$Barcode == bc]); ft <- Lu$ftype[Lu$Barcode == bc][1]
      lf <- lf + log(1 - odds_mult(ce$r, if (ft == "capsule") f$kap_capsule else 1) * (1 - (1 - ce$p)^M)) }
    lm <- lm + log(1 - ce$q * (1 - exp(lf))) }
  1 - exp(lm) }
build_pairs <- function(m, ev) { f <- fits[[paste(m, a)]]; Lf <- f$L
  Lf$unit <- ifelse(Lf$site == "COAST", Lf$Site, ifelse(grepl(paste0(" ", ev, "$"), Lf$visit_id), Lf$Site, NA)); Lf <- Lf[!is.na(Lf$unit), ]
  units <- sort(unique(Lf$unit)); post <- f$post / rowSums(f$post)
  d <- sapply(units, function(u) as.vector(post %*% unit_det(f, Lf[Lf$unit == u, ]))); rownames(d) <- f$taxa
  y <- sapply(units, function(u) f$taxa %in% f$pres$taxon[f$pres$library_id %in% Lf$library_id[Lf$unit == u]]); rownames(y) <- f$taxa
  G <- expand.grid(A = units, B = units, stringsAsFactors = FALSE) %>% filter(A != B)
  R <- bind_rows(lapply(seq_len(nrow(G)), function(i) { s <- f$taxa[y[, G$A[i]]]; if (!length(s)) return(NULL)
    data.frame(A = G$A[i], B = G$B[i], taxon = s, yB = y[s, G$B[i]], dB = d[s, G$B[i]]) }))
  R$marker <- m; R$dist_coast <- CD[cbind(R$A, R$B)]; R$dist_gc <- GC[cbind(R$A, R$B)]
  sd <- setNames(LOC$side, LOC$location); R$cross <- as.numeric(sd[R$A] != sd[R$B])
  list(R = R, units = units, d = d, y = y) }
# ---------- fit ----------
pi_fn <- function(th, dist, cross) pmin(0.999, (plogis(th[1]) + (plogis(th[2]) - plogis(th[1])) * exp(-dist / exp(th[3]))) * exp(th[4] * cross))
nll <- function(th, R, dcol, step = TRUE) { if (!step) th[4] <- 0; pr <- pmin(1 - 1e-9, pmax(1e-9, R$dB * pi_fn(th, R[[dcol]], R$cross)))
  -sum(R$yB * log(pr) + (1 - R$yB) * log(1 - pr)) }
fit1 <- function(R, dcol, step = TRUE, Lfix = NULL, start = c(-1, 1, log(100), 0)) {
  if (is.null(Lfix)) { o <- optim(start, nll, R = R, dcol = dcol, step = step, method = "L-BFGS-B", lower = c(-8, -8, log(0.5), -3), upper = c(8, 8, log(20000), 3))
    for (s3 in log(c(5, 50, 500, 5000))) { o2 <- optim(replace(start, 3, s3), nll, R = R, dcol = dcol, step = step, method = "L-BFGS-B", lower = c(-8, -8, log(0.5), -3), upper = c(8, 8, log(20000), 3)); if (o2$value < o$value) o <- o2 }
    return(list(th = o$par, nll = o$value)) }
  g <- function(t3) nll(c(t3[1:2], log(Lfix), t3[3]), R, dcol, step)
  o <- optim(start[c(1, 2, 4)], g, method = "L-BFGS-B", lower = c(-8, -8, -3), upper = c(8, 8, 3)); list(th = c(o$par[1:2], log(Lfix), o$par[3]), nll = o$value) }
LGRID <- exp(seq(log(0.5), log(20000), length.out = 40))
analyse <- function(R, dcol, label, boot = TRUE, jack = TRUE) {
  full <- fit1(R, dcol); nostep <- fit1(R, dcol, step = FALSE); flat <- { o <- optim(0, function(t) nll(c(t, t, 0, 0), R, dcol), method = "Brent", lower = -8, upper = 8); list(nll = o$value) }
  prof <- sapply(LGRID, function(Lc) fit1(R, dcol, Lfix = Lc, start = full$th)$nll); pdrop <- prof - full$nll
  inside <- LGRID[pdrop <= 1.92]
  bt <- if (boot) { sp <- unique(R$taxon); idx <- split(seq_len(nrow(R)), R$taxon)
    t(sapply(seq_len(NBOOT), function(b) { s <- sample(sp, replace = TRUE); Rb <- R[unlist(idx[s], use.names = FALSE), ]; fit1(Rb, dcol, start = full$th)$th })) } else NULL
  jk <- if (jack) { u <- unique(c(R$A, R$B)); sapply(u, function(z) exp(fit1(R[R$A != z & R$B != z, ], dcol, start = full$th)$th[3])) } else NULL
  th <- full$th
  list(summary = data.frame(label = label, distance = dcol, rows = nrow(R), species = n_distinct(R$taxon),
    pi_0 = plogis(th[2]), pi_inf = plogis(th[1]), Lc_km = exp(th[3]), step_multiplier = exp(th[4]),
    Lc_profile_lo = min(inside), Lc_profile_hi = max(inside), Lc_profile_open_hi = max(pdrop[LGRID > exp(th[3])]) <= 1.92,
    Lc_profile_open_lo = max(pdrop[LGRID < exp(th[3])], 0) <= 1.92,
    Lc_boot_lo = if (boot) unname(quantile(exp(bt[, 3]), 0.025)) else NA, Lc_boot_hi = if (boot) unname(quantile(exp(bt[, 3]), 0.975)) else NA,
    boot_at_bound = if (boot) mean(bt[, 3] > log(19000) | bt[, 3] < log(0.6)) else NA,
    step_boot_lo = if (boot) unname(quantile(exp(bt[, 4]), 0.025)) else NA, step_boot_hi = if (boot) unname(quantile(exp(bt[, 4]), 0.975)) else NA,
    Lc_jack_min = if (jack) min(jk) else NA, Lc_jack_max = if (jack) max(jk) else NA,
    dll_decay_vs_flat_nostep = flat$nll - nostep$nll, dll_step = nostep$nll - full$nll),
    profile = data.frame(label = label, distance = dcol, Lc_km = LGRID, drop = pdrop)) }
bins <- function(R, dcol, label) R %>% mutate(bin = cut(.data[[dcol]], c(-1, 15, 100, 300, 700, 1500, 1e5), labels = c("<15", "15-100", "100-300", "300-700", "700-1500", ">1500"))) %>%
  group_by(label = label, distance = dcol, bin, cross) %>% summarise(ordered_pairs = n_distinct(paste(A, B)), species_rows = n(), observed_shared = sum(yB),
    expected_if_all_present = sum(dB), corrected_shared_fraction = sum(yB) / sum(dB), raw_shared_fraction = mean(yB), .groups = "drop")
if (!DEFS_ONLY) {
SUM <- list(); PRO <- list(); BIN <- list(); EVS <- list()
evs <- sort(unique(sub("^.* (Event\\d+)$", "\\1", B$L$visit_id[B$L$site %in% c("GP", "CJ", "JB")])))
for (m in c(MARKERS, "12S+COI")) {
  pp <- if (m == "12S+COI") lapply(MARKERS, function(z) build_pairs(z, PC_EVENT)$R) else list(build_pairs(m, PC_EVENT)$R)
  R <- bind_rows(pp) %>% mutate(taxon = paste(marker, taxon)); lg(m, PC_EVENT, ": rows", nrow(R), " units", n_distinct(R$A))
  for (dc in c("dist_coast", "dist_gc")) { z <- analyse(R, dc, m, boot = dc == "dist_coast", jack = dc == "dist_coast"); SUM[[length(SUM) + 1]] <- z$summary; PRO[[length(PRO) + 1]] <- z$profile
    BIN[[length(BIN) + 1]] <- bins(R, dc, m); lg(m, dc, ": Lc", round(z$summary$Lc_km), " step", round(z$summary$step_multiplier, 2)) }
  isl <- LOC$location[LOC$island]; z <- analyse(R[!R$A %in% isl & !R$B %in% isl, ], "dist_coast", paste(m, "mainland only"), boot = FALSE, jack = FALSE)
  SUM[[length(SUM) + 1]] <- z$summary
  if (m != "12S+COI") for (ev in evs) { Re <- build_pairs(m, ev)$R; if (n_distinct(Re$A[Re$A %in% names(PC_SITES)]) < 1) next
    z <- analyse(Re %>% mutate(taxon = paste(marker, taxon)), "dist_coast", paste(m, ev), boot = FALSE, jack = FALSE)
    EVS[[length(EVS) + 1]] <- z$summary %>% mutate(event = ev, pc_units = paste(intersect(unique(Re$A), names(PC_SITES)), collapse = "/")) } }
SUM <- bind_rows(SUM); wcsv(SUM, "t_turnover_summary"); wcsv(bind_rows(PRO), "t_turnover_profile"); wcsv(bind_rows(BIN), "t_turnover_bins"); wcsv(bind_rows(EVS), "t_turnover_by_pc_event")
print(LOC %>% select(location, habitat_class, latitude, longitude, visits, samples_bags, filters_tubes, island, side, coast_pos_km) %>% as.data.frame())
print(SUM %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())
print(bind_rows(EVS) %>% select(label, Lc_km, Lc_profile_lo, Lc_profile_hi, step_multiplier, dll_step) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())
print(bind_rows(BIN) %>% filter(distance == "dist_coast") %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())
}
