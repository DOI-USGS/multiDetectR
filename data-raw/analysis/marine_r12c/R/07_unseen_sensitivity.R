# 07_unseen_sensitivity.R -- peer advice 2026-10-10: the observed-species design is the main result; rerun with the UNSEEN species added
# at low / middle / high totals and say whether the design changes.
# The objective is linear in the species weights, so for a total unseen U the score of any design is S_U = S + (U / U_fit) * (S_all - S)
# (S = observed species; S_all - S = the fitted unseen mass nu*P0, with the shape over species types taken from the fit).
# Stated prior on U (the likelihood is flat in U):  low = Chao2 lower bound from site-visit incidence;  middle = the fit's U;  high = 2 x the fit's U.
# Reads the saved design points (03); no refit.
source("common.R")
fits <- readRDS(paste0(FITDIR, "fits.rds")); B <- readRDS(paste0(FITDIR, "build.rds")); a <- as.numeric(Sys.getenv("A", "0.03"))
chao2 <- function(m) { L <- B$L[B$L$marker == m, ]; pr <- B$pres[[m]] %>% left_join(L %>% distinct(library_id, site, visit_id), by = "library_id")
  inc <- pr %>% distinct(taxon, site, visit_id) %>% count(taxon); T <- nrow(distinct(L, site, visit_id)); Q1 <- sum(inc$n == 1); Q2 <- sum(inc$n == 2)
  data.frame(marker = m, units = T, observed = nrow(inc), Q1 = Q1, Q2 = Q2, unseen_chao2 = (T - 1) / T * Q1^2 / (2 * Q2)) }
U <- bind_rows(lapply(MARKERS, chao2)) %>% mutate(unseen_fit = sapply(marker, function(m) fits[[paste(m, a)]]$Ntot - fits[[paste(m, a)]]$Nobs))
sc <- bind_rows(lapply(MARKERS, function(m) { u <- U[U$marker == m, ]
  data.frame(row.names = NULL, marker = m, scenario = c("observed only", "low (Chao2)", "middle (fit)", "high (2 x fit)"), U = c(0, u$unseen_chao2, u$unseen_fit, 2 * u$unseen_fit), U_fit = u$unseen_fit) }))
wcsv(U, sprintf("t_unseen_prior_a%s", a)); print(U)
BUD <- c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000)
P <- readRDS(paste0(OUT, sprintf("pc_points_a%s_disc.rds", a))); CP <- readRDS(paste0(OUT, sprintf("coast_points_a%s_disc.rds", a)))
lam <- function(mk, s) { if (mk == "12S+COI") stop(); x <- sc[sc$marker == mk & sc$scenario == s, ]; x$U / x$U_fit }
out <- list()
for (s in unique(sc$scenario)) for (mk in c("12S", "COI", "12S+COI")) {
  X <- P[P$markers == mk, ]
  # union: split the unseen part by marker using the single-marker increments (shared-name overlap is among OBSERVED species only)
  X$SU <- if (mk == "12S+COI") { s12 <- P[P$markers == "12S", ]; sco <- P[P$markers == "COI", ]
    g12 <- (s12$S_all - s12$S)[match(paste(X$i, X$k12), paste(s12$i, s12$k12))]; gco <- (sco$S_all - sco$S)[match(paste(X$i, X$kCOI), paste(sco$i, sco$kCOI))]
    X$S + lam("12S", s) * g12 + lam("COI", s) * gco } else X$S + lam(mk, s) * (X$S_all - X$S)
  for (Bd in BUD) { Y <- X[X$cost <= Bd, ]; if (!nrow(Y)) next; b <- Y[which.max(Y$SU), ]
    out[[length(out) + 1]] <- data.frame(scope = "Pt. Conception", scenario = s, markers = mk, budget = Bd, b[, c("V", "sites", "GP", "CJ", "JB", "f", "k12", "kCOI", "cost")],
                                         S_observed = b$S, S_scenario = b$SU, N = NA, b_loc = NA, k = NA) } }
for (s in unique(sc$scenario)) for (m in MARKERS) { X <- CP[CP$markers == m, ]; X$SU <- X$S + lam(m, s) * (X$S_all - X$S)
  for (Bd in c(BUD, 15000, 20000, 30000, 40000)) { Y <- X[X$cost <= Bd, ]; if (!nrow(Y)) next; b <- Y[which.max(Y$SU), ]
    out[[length(out) + 1]] <- data.frame(scope = "coastline", scenario = s, markers = m, budget = Bd, b_loc = b$b, N = b$N, f = b$f, k = b$k, cost = b$cost, S_observed = b$S, S_scenario = b$SU) } }
R <- bind_rows(out); R$a <- a; wcsv(R, sprintf("t_unseen_sensitivity_a%s", a))
key <- function(d) with(d, ifelse(scope == "coastline", paste0("N", N, " b", b_loc, " f", f, " k", k), paste0("V", V, " ", sites, " ", GP, "/", CJ, "/", JB, " f", f, " k", ifelse(is.na(k12), "", k12), "|", ifelse(is.na(kCOI), "", kCOI))))
R$design <- key(R); W <- R %>% select(scope, markers, budget, scenario, design) %>% tidyr::pivot_wider(names_from = scenario, values_from = design)
print(as.data.frame(W), right = FALSE)
cat("\nloss in OBSERVED species from designing for the scenario instead (max over budgets):\n")
print(R %>% group_by(scope, markers, budget) %>% mutate(loss = S_observed[scenario == "observed only"] - S_observed) %>% group_by(scope, markers, scenario) %>% summarise(max_loss = round(max(loss), 2), .groups = "drop") %>% as.data.frame())
