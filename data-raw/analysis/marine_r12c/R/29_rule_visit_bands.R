# 29_rule_visit_bands.R -- Kevin 2026-10-11 (ZFF report): the 12S+COI rule (25 v2, template fit, a = 0.03, hour costs) at the budgets where it
# chooses V = 1, 2, 3, 4, 6, 12 visits. 25 ran only 10 budgets (it jumps 4 -> 7 visits between $5,000 and $7,500), so this sweeps the budget:
# coarse ($250 steps, $500-$18,000), then bisection to $5 (the cost unit) at every change in the chosen V. For each budget it keeps the rule's best
# design at EVERY V (the rule's choice is the V with the most expected species), so the band of budgets choosing each V is known.
# Reported design per V: the rule design at the TOP of its band (largest budget still choosing V), which is how the existing 1-4 visit rows were
# picked from 25's budget list ($1,500, $3,000, $4,000 -> cost $3,975, $5,000 were the largest of those budgets choosing 1-4 visits).
# The rule procedure (steps 1-4, v2 additions a-b) is copied verbatim from 25_rule_designs.R; only the V loop keeps every V. Exact, no simulation.
# Env FITS (default fits_template.rds). Fork workers: NCORE (default 6); each worker is small (fit 5 MB serialized).
source("common.R"); source("model.R")
library(parallel)
FN <- Sys.getenv("FITS", "fits_template.rds"); fits <- readRDS(paste0(FITDIR, FN)); C <- read_costs(); a <- 0.03
SFX <- if (FN == "fits.rds") "" else paste0("_", sub("^fits_(.*)\\.rds$", "\\1", FN)); NCORE <- as.integer(Sys.getenv("NCORE", "6"))
SITES3 <- c("GP", "CJ", "JB"); KMAX <- 60; CAP <- 33; MK <- "12S+COI"; ms <- c("12S", "COI"); TARGET_V <- c(1, 2, 3, 4, 6, 12)
cf <- C[["filter_disc"]]; cl <- C[["library"]]; cp <- C[["pcr"]]; cs <- C[["sample"]]
FF <- list(`12S` = fits[[paste("12S", a)]], COI = fits[[paste("COI", a)]]); W <- lapply(FF, function(f) colSums(f$post))
shared <- intersect(FF$`12S`$taxa, FF$COI$taxa)
# ---- verbatim from 25: inner pass, exact expectation, visit cost
nd_of <- function(ce, k, tm) if (is.null(tm)) (1 - ce$p)^m_eff(k, a) else (1 - m_eff(k, a) * log(1 - ce$p) / tm$alpha)^(-tm$alpha)
inner_props <- function(m) { f <- FF[[m]]; ce <- f$cells; tm <- f$tmpl; keep <- which(W[[m]] > 1e-3)
  K <- 1:KMAX; Pf <- sapply(K, function(k) odds_mult(ce$r[keep], 1) * (1 - nd_of(ce[keep, ], k, tm)))
  kap <- t(t(1 / -log(1 - Pf)) * (cf + cl + cp * K)); kk <- K[apply(kap, 1, which.min)]
  ksh <- K[apply(t(t(1 / -log(1 - Pf)) * (cl + cp * K)), 1, which.min)]
  Pk <- Pf[cbind(seq_along(keep), kk)]; Cf <- cf + cl + cp * kk
  ff <- sapply(seq_along(keep), function(i) { fs <- 1:6; Hb <- -log(1 - ce$q[keep[i]] * (1 - (1 - Pk[i])^fs)); fs[which.min((cs + fs * Cf[i]) / Hb)] })
  data.frame(marker = m, type = keep, weight = W[[m]][keep], k = kk, f = ff, k_shared = ksh) }
PROP <- bind_rows(lapply(MARKERS, inner_props))
site_lnot <- function(m, f, k, b) { ce <- FF[[m]]$cells; pb <- ce$q * (1 - (1 - filt_det(ce, a, k, "disc", FF[[m]]$kap_capsule))^f)
  sapply(SITES3, function(s) log(pmax(1 - ce[[paste0("phi_", s)]] * (1 - (1 - pb)^b[[s]]), 1e-300))) }
S_design <- function(mk, V, b, f, k) { ms <- strsplit(mk, "\\+")[[1]]; Sm <- 0; Dm <- list()
  for (m in ms) { P <- 1 - exp(V * rowSums(site_lnot(m, f, k[[m]], b))); Sm <- Sm + sum(W[[m]] * P); Dm[[m]] <- (FF[[m]]$post %*% P)[match(shared, FF[[m]]$taxa)] }
  if (length(ms) == 2) Sm <- Sm - sum(Dm[[1]] * Dm[[2]]); Sm }
visit_cost <- function(mk, b, f, k) { ms <- strsplit(mk, "\\+")[[1]]; nb <- sum(b); ns <- sum(b > 0)
  C[["travel_location"]] + C[["travel_site"]] * max(0, ns - 1) + nb * (cs + f * (cf + sum(sapply(ms, function(m) cl + cp * k[[m]])))) }
# ---- 25's rule_for_inner, returning the result for every V instead of only the best V
rule_for_inner_allV <- function(mk, B, f, k) { ms <- strsplit(mk, "\\+")[[1]]; out <- list()
  for (V in seq_len(floor(B / (C[["travel_location"]] + cs + f * (cf + length(ms) * (cl + cp)))))) {
    Bv <- B / V; b <- c(GP = 0, CJ = 0, JB = 0); S0 <- 0; order_added <- character(0)
    repeat { cand <- lapply(SITES3, function(s) { b2 <- b; b2[s] <- b2[s] + 1; cst <- visit_cost(mk, b2, f, k); if (cst > Bv) return(NULL)
        list(s = s, b = b2, gain = (S_design(mk, V, b2, f, k) - S0) / (V * (cst - (if (sum(b) > 0) visit_cost(mk, b, f, k) else 0)))) })
      cand <- Filter(Negate(is.null), cand); if (!length(cand)) break
      j <- which.max(sapply(cand, `[[`, "gain")); b <- cand[[j]]$b; S0 <- S_design(mk, V, b, f, k)
      if (!cand[[j]]$s %in% order_added) order_added <- c(order_added, cand[[j]]$s) }
    if (sum(b) == 0) next
    kk <- k; repeat { left <- Bv - visit_cost(mk, b, f, kk); step <- sum(b) * f * cp; if (left < step) break
      opts <- Filter(function(m) sum(unlist(kk[ms])) < CAP, ms); if (!length(opts)) break
      g <- sapply(opts, function(m) { k2 <- kk; k2[[m]] <- k2[[m]] + 1; S_design(mk, V, b, f, k2) }); m <- opts[which.max(g)]; kk[[m]] <- kk[[m]] + 1 }
    out[[V]] <- data.frame(V = V, GP = b[["GP"]], CJ = b[["CJ"]], JB = b[["JB"]], f = f, k12 = kk$`12S`, kCOI = kk$COI, S = S_design(mk, V, b, f, kk),
      cost = V * visit_cost(mk, b, f, kk), site_order = paste(order_added, collapse = ">")) }
  bind_rows(out) }
p1 <- PROP %>% filter(marker == "12S") %>% distinct(f, k12 = k); p2 <- PROP %>% filter(marker == "COI") %>% distinct(f, kCOI = k)
s1 <- PROP %>% filter(marker == "12S") %>% distinct(k12 = k_shared); s2 <- PROP %>% filter(marker == "COI") %>% distinct(kCOI = k_shared)
PROPS <- bind_rows(merge(p1 %>% select(k12), p2, by = NULL), merge(p1, p2 %>% select(kCOI), by = NULL),
                   merge(merge(s1, s2, by = NULL), data.frame(f = unique(c(p1$f, p2$f))), by = NULL)) %>% distinct(f, k12, kCOI)
lg(MK, ": distinct inner proposals", nrow(PROPS), "(25 v2 had 1286)")
# rule at budget B: best design per V over the proposals; the rule's choice = the V with the largest S (ties: first, as 25's strict >)
rule_allV <- function(B) { nb_aff <- floor((B - C[["travel_location"]]) / (cs + cf + 2 * (cl + cp)))
  pr <- if (nb_aff < 10) bind_rows(PROPS, expand.grid(f = 1:2, k12 = 1:20, kCOI = 1:20)) %>% distinct(f, k12, kCOI) else PROPS
  X <- bind_rows(lapply(seq_len(nrow(pr)), function(i) { r <- rule_for_inner_allV(MK, B, pr$f[i], list(`12S` = pr$k12[i], COI = pr$kCOI[i])); if (nrow(r)) r$prop <- i; r }))
  X <- X %>% group_by(V) %>% slice(which.max(S)) %>% ungroup() %>% mutate(budget = B, pcr_enumerated = nb_aff < 10)
  X$chosen <- seq_len(nrow(X)) == which.max(X$S); X }
if (nzchar(Sys.getenv("TEST_B"))) { print(as.data.frame(rule_allV(as.numeric(Sys.getenv("TEST_B"))))); quit(save = "no") }   # check vs 25 v2
CACHE <- new.env()
run_budgets <- function(Bs) { Bs <- setdiff(Bs, as.numeric(ls(CACHE))); if (!length(Bs)) return(invisible())
  t0 <- Sys.time(); R <- mclapply(Bs, rule_allV, mc.cores = NCORE, mc.preschedule = FALSE)
  bad <- which(!sapply(R, is.data.frame)); if (length(bad)) stop("worker failed at budget ", Bs[bad[1]], ": ", R[[bad[1]]])
  for (i in seq_along(Bs)) assign(as.character(Bs[i]), R[[i]], envir = CACHE)
  lg(length(Bs), "budgets in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min") }
chosenV <- function(B) { x <- get(as.character(B), envir = CACHE); x$V[x$chosen] }
# ---- coarse sweep
COARSE <- seq(500, 18000, by = 250); run_budgets(COARSE)
cv <- sapply(COARSE, chosenV); lg("coarse chosen V:", paste(sprintf("%d:%d", COARSE, cv), collapse = " "))
# ---- bisection to $5 at each change in chosen V between neighbouring coarse budgets
for (i in which(diff(cv) != 0)) { lo <- COARSE[i]; hi <- COARSE[i + 1]
  while (hi - lo > 5) { mids <- unique(round(seq(lo, hi, length.out = NCORE + 2)[-c(1, NCORE + 2)] / 5) * 5); mids <- mids[mids > lo & mids < hi]
    run_budgets(mids); grid <- c(lo, mids, hi); v <- sapply(grid, chosenV); j <- which(v != v[1])[1]; lo <- grid[j - 1]; hi <- grid[j] } }
ALL <- bind_rows(mget(ls(CACHE), envir = CACHE)) %>% arrange(budget, V)
wcsv(ALL, paste0("t_rule_visit_bands_all", SFX))
CH <- ALL %>% filter(chosen) %>% arrange(budget)
BAND <- CH %>% mutate(run = cumsum(c(1, diff(V) != 0))) %>% group_by(run, V) %>% summarise(budget_from = min(budget), budget_to = max(budget), .groups = "drop")
wcsv(BAND, paste0("t_rule_visit_bands", SFX)); print(as.data.frame(BAND))
# reported row per target V: the rule design at the largest searched budget choosing V; if the rule never chooses V, say so
REP <- bind_rows(lapply(TARGET_V, function(v) { x <- CH %>% filter(V == v)
  if (!nrow(x)) return(data.frame(V_target = v, note = "rule never chooses this V in $500-$18,000"))
  b <- x[which.max(x$budget), ]; data.frame(V_target = v, budget = b$budget, band = paste(BAND$budget_from[BAND$V == v], BAND$budget_to[BAND$V == v], sep = "-", collapse = "; "),
    GP = b$GP, CJ = b$CJ, JB = b$JB, filters_per_bag = b$f, pcr_12S = b$k12, pcr_COI = b$kCOI, cost_total = b$cost, cost_per_visit = b$cost / v, S_expected_observed = b$S,
    site_order = b$site_order, note = "") }))
wcsv(REP, paste0("t_rule_visits_fixed_report", SFX)); print(REP %>% mutate(across(where(is.double), ~ round(.x, 2))) %>% as.data.frame(), right = FALSE)
