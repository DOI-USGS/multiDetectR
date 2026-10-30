# 25_rule_designs.R -- peer (paper chat) 2026-10-11: every reported Pt. Conception design should come from the RULE, with the exhaustive grid
# (03_designs.R) only as a check. Exact expected values throughout (closed forms per species type, summed over the fitted mixture); no simulation.
# The rule, per marker set (markers are a study-wide choice) and budget B (costs.csv, technician hours at $30/h; a = 0.03):
#  1. Inner levels by the inward pass (paperA S1.1-S1.3), integer counts by rounding in place (argmin of the price at each level):
#     PCR: k_m = argmin_k (c_filter + c_library + c_pcr k) / H_f(k), H_f(k) = -log[1 - r'(1 - nd(k))], nd(k) = miss probability of a pooled
#          library of k reactions: (1 + beta m_eff(k, a))^(-alpha) with the template gamma (beta = -log(1 - p) / alpha), (1 - p)^m_eff without it.
#     filter: f = argmin_f (c_sample + f C_filter) / H_bag(f), H_bag(f) = -log[1 - q (1 - (1 - P_filter)^f)].
#     Done for an effective single target: every fitted species type with weight > 1e-3 proposes its own (f, k); the distinct proposals are
#     carried into steps 2-4 and the one giving the most expected species is kept (as in 09_rule_vs_grid.R).
#  2. Branches (sites GP, CJ, JB): bags added one at a time to the site with the largest marginal expected species per dollar, the marginal test of
#     S1.5 (eq. S1.8) with the concave site hazard H_s(b) = -log[1 - phi_s (1 - (1 - P_bag)^b)] of a species findable at several sites. The first bag
#     at a second or third site also pays travel_site. lambda_t (still undetected) comes from the exact expectation.
#  3. Visits: one integer loop over V; each visit gets B / V - travel_location for steps 2 and 4.
#  4. Fill: money left after the last affordable bag goes inward to PCRs, one PCR per filter at a time to the marker with the larger marginal
#     species per dollar (extract cap: 33 PCRs per filter summed over markers).
# v2 (2026-10-11, after the first check against 27): two parts of the procedure the first version left out.
#  (a) Shared filter: in a two-marker design each marker's library is a branch inside a filter that exists anyway, so its PCR count is priced on
#      its own cost: k_m = argmin_k (c_library + c_pcr k) / H_f(k). These proposals are added to the single-marker ones.
#  (b) Budget-aware enumeration (11_enumerate_small_levels.R, trigger 10): when fewer than 10 bags are affordable in total, the PCR level is
#      enumerated (k = 1-20 per marker) instead of taken from the price.
# Check: the grid optimum from 03 on the same fit (env FITS / COSTTAG as for 03). Env FITS (default fits.rds), GRIDTAG (03 output folder tag);
# RULEVER (default "_v2") is appended to the output names.
source("common.R"); source("model.R")
FN <- Sys.getenv("FITS", "fits.rds"); fits <- readRDS(paste0(FITDIR, FN)); C <- read_costs(); a <- 0.03
SFX <- if (FN == "fits.rds") "" else paste0("_", sub("^fits_(.*)\\.rds$", "\\1", FN))
GRIDDIR <- if (nzchar(Sys.getenv("GRIDTAG"))) paste0(FITDIR, "costs_", Sys.getenv("GRIDTAG"), "/") else OUT
BUDGETS <- c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000); SITES3 <- c("GP", "CJ", "JB"); KMAX <- 60; CAP <- 33
cf <- C[["filter_disc"]]; cl <- C[["library"]]; cp <- C[["pcr"]]; cs <- C[["sample"]]
FF <- list(`12S` = fits[[paste("12S", a)]], COI = fits[[paste("COI", a)]]); W <- lapply(FF, function(f) colSums(f$post))
shared <- intersect(FF$`12S`$taxa, FF$COI$taxa)
# ---- 1. inner pass per species type
nd_of <- function(ce, k, tm) if (is.null(tm)) (1 - ce$p)^m_eff(k, a) else (1 - m_eff(k, a) * log(1 - ce$p) / tm$alpha)^(-tm$alpha)
inner_props <- function(m) { f <- FF[[m]]; ce <- f$cells; tm <- f$tmpl; keep <- which(W[[m]] > 1e-3)
  K <- 1:KMAX; Pf <- sapply(K, function(k) odds_mult(ce$r[keep], 1) * (1 - nd_of(ce[keep, ], k, tm)))          # types x k
  kap <- t(t(1 / -log(1 - Pf)) * (cf + cl + cp * K)); kk <- K[apply(kap, 1, which.min)]
  ksh <- K[apply(t(t(1 / -log(1 - Pf)) * (cl + cp * K)), 1, which.min)]                                     # (a) shared filter
  Pk <- Pf[cbind(seq_along(keep), kk)]; Cf <- cf + cl + cp * kk
  ff <- sapply(seq_along(keep), function(i) { fs <- 1:6; Hb <- -log(1 - ce$q[keep[i]] * (1 - (1 - Pk[i])^fs)); fs[which.min((cs + fs * Cf[i]) / Hb)] })
  data.frame(marker = m, type = keep, weight = W[[m]][keep], k = kk, f = ff, k_shared = ksh) }
PROP <- bind_rows(lapply(MARKERS, inner_props))
# ---- exact expected species of a design
site_lnot <- function(m, f, k, b) { ce <- FF[[m]]$cells; pb <- ce$q * (1 - (1 - filt_det(ce, a, k, "disc", FF[[m]]$kap_capsule))^f)
  sapply(SITES3, function(s) log(pmax(1 - ce[[paste0("phi_", s)]] * (1 - (1 - pb)^b[[s]]), 1e-300))) }        # cells x site, one visit
S_design <- function(mk, V, b, f, k) { ms <- strsplit(mk, "\\+")[[1]]; Sm <- 0; Dm <- list()
  for (m in ms) { P <- 1 - exp(V * rowSums(site_lnot(m, f, k[[m]], b))); Sm <- Sm + sum(W[[m]] * P); Dm[[m]] <- (FF[[m]]$post %*% P)[match(shared, FF[[m]]$taxa)] }
  if (length(ms) == 2) Sm <- Sm - sum(Dm[[1]] * Dm[[2]]); Sm }
visit_cost <- function(mk, b, f, k) { ms <- strsplit(mk, "\\+")[[1]]; nb <- sum(b); ns <- sum(b > 0)
  C[["travel_location"]] + C[["travel_site"]] * max(0, ns - 1) + nb * (cs + f * (cf + sum(sapply(ms, function(m) cl + cp * k[[m]])))) }
# ---- 2-4 for one inner design (f, k), all V
rule_for_inner <- function(mk, B, f, k) { ms <- strsplit(mk, "\\+")[[1]]; best <- NULL
  for (V in seq_len(floor(B / (C[["travel_location"]] + cs + f * (cf + length(ms) * (cl + cp)))))) {
    Bv <- B / V; b <- c(GP = 0, CJ = 0, JB = 0); S0 <- 0; order_added <- character(0)
    repeat { cand <- lapply(SITES3, function(s) { b2 <- b; b2[s] <- b2[s] + 1; cst <- visit_cost(mk, b2, f, k); if (cst > Bv) return(NULL)
        list(s = s, b = b2, gain = (S_design(mk, V, b2, f, k) - S0) / (V * (cst - (if (sum(b) > 0) visit_cost(mk, b, f, k) else 0)))) })
      cand <- Filter(Negate(is.null), cand); if (!length(cand)) break
      j <- which.max(sapply(cand, `[[`, "gain")); b <- cand[[j]]$b; S0 <- S_design(mk, V, b, f, k)
      if (!cand[[j]]$s %in% order_added) order_added <- c(order_added, cand[[j]]$s) }
    if (sum(b) == 0) next
    kk <- k; repeat { left <- Bv - visit_cost(mk, b, f, kk); step <- sum(b) * f * cp; if (left < step) break                     # 4. fill
      opts <- Filter(function(m) sum(unlist(kk[ms])) < CAP, ms); if (!length(opts)) break
      g <- sapply(opts, function(m) { k2 <- kk; k2[[m]] <- k2[[m]] + 1; S_design(mk, V, b, f, k2) }); m <- opts[which.max(g)]; kk[[m]] <- kk[[m]] + 1 }
    S <- S_design(mk, V, b, f, kk)
    if (is.null(best) || S > best$S) best <- list(V = V, b = b, f = f, k = kk, S = S, cost = V * visit_cost(mk, b, f, kk), site_order = paste(order_added, collapse = ">")) }
  best }
res <- list()
for (mk in c("12S", "COI", "12S+COI")) { ms <- strsplit(mk, "\\+")[[1]]
  props <- if (length(ms) == 1) PROP %>% filter(marker == ms) %>% distinct(f, k) %>% mutate(k12 = if (ms == "12S") k else NA, kCOI = if (ms == "COI") k else NA) else
    { p1 <- PROP %>% filter(marker == "12S") %>% distinct(f, k12 = k); p2 <- PROP %>% filter(marker == "COI") %>% distinct(f, kCOI = k)
      s1 <- PROP %>% filter(marker == "12S") %>% distinct(k12 = k_shared); s2 <- PROP %>% filter(marker == "COI") %>% distinct(kCOI = k_shared)
      bind_rows(merge(p1 %>% select(k12), p2, by = NULL), merge(p1, p2 %>% select(kCOI), by = NULL),
                merge(merge(s1, s2, by = NULL), data.frame(f = unique(c(p1$f, p2$f))), by = NULL)) %>% distinct(f, k12, kCOI) }
  lg(mk, ": distinct inner proposals", nrow(props))
  for (Bd in BUDGETS) { best <- NULL
    nb_aff <- floor((Bd - C[["travel_location"]]) / (cs + cf + length(ms) * (cl + cp)))                          # (b) bags affordable at 1 PCR
    pr_b <- if (nb_aff < 10) { e <- expand.grid(f = 1:2, k12 = if ("12S" %in% ms) 1:20 else NA, kCOI = if ("COI" %in% ms) 1:20 else NA)
      bind_rows(props, e) %>% distinct(f, k12, kCOI) } else props
    props_all <- props; props <- pr_b
    for (i in seq_len(nrow(props))) { k <- list(`12S` = props$k12[i], COI = props$kCOI[i]); r <- rule_for_inner(mk, Bd, props$f[i], k)
      if (!is.null(r) && (is.null(best) || r$S > best$S)) best <- r }
    props <- props_all
    res[[paste(mk, Bd)]] <- data.frame(markers = mk, budget = Bd, pcr_enumerated = nb_aff < 10, V = best$V, GP = best$b[["GP"]], CJ = best$b[["CJ"]], JB = best$b[["JB"]], f = best$f,
      k12 = if ("12S" %in% ms) best$k$`12S` else NA, kCOI = if ("COI" %in% ms) best$k$COI else NA, S_rule = best$S, cost_rule = best$cost, site_order = best$site_order)
    lg(mk, Bd, ": rule V", best$V, paste(best$b, collapse = "/"), "f", best$f, "k", paste(unlist(best$k[ms]), collapse = "|"), "S", round(best$S, 2)) } }
R <- bind_rows(res)
# ---- check against the grid (03 on the same fit)
P <- readRDS(paste0(GRIDDIR, sprintf("pc_points_a%s_disc.rds", a)))
G <- bind_rows(lapply(seq_len(nrow(R)), function(i) { X <- P[P$markers == R$markers[i] & P$cost <= R$budget[i], ]; b <- X[which.max(X$S), ]
  data.frame(grid_design = sprintf("V%d %d/%d/%d f%d k%s|%s", b$V, b$GP, b$CJ, b$JB, b$f, ifelse(is.na(b$k12), "", b$k12), ifelse(is.na(b$kCOI), "", b$kCOI)),
             S_grid = b$S, cost_grid = b$cost) }))
R <- bind_cols(R, G) %>% mutate(rule_design = sprintf("V%d %d/%d/%d f%d k%s|%s", V, GP, CJ, JB, f, ifelse(is.na(k12), "", k12), ifelse(is.na(kCOI), "", kCOI)),
  same_design = rule_design == grid_design, loss_species = S_grid - S_rule,
  rule_outside_grid = !(V %in% c(1:6, 8, 10, 12, 16)) | GP > 12 | CJ > 6 | JB > 6 | f > 3 |
    (!is.na(k12) & !(k12 %in% c(1, 2, 3, 4, 6, 8, 10, 12, 15, 20, 25, 30))) | (!is.na(kCOI) & !(kCOI %in% c(1, 2, 3, 4, 6, 8, 10, 12, 15, 20, 25, 30))))
RV <- Sys.getenv("RULEVER", "_v2"); wcsv(R, paste0("t_rule_designs", SFX, RV)); wcsv(PROP, paste0("t_rule_inner_proposals", SFX, RV))
print(R %>% select(markers, budget, rule_design, S_rule, grid_design, S_grid, loss_species, rule_outside_grid, site_order) %>% mutate(across(where(is.double), ~ round(.x, 2))) %>% as.data.frame(), right = FALSE)
print(R %>% group_by(markers) %>% summarise(same = sum(same_design), n = n(), max_loss = max(loss_species), mean_loss = mean(loss_species), rule_beats_grid = sum(loss_species < -1e-9)) %>% as.data.frame())
