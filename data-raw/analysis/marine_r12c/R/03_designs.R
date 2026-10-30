# 03_designs.R -- cost vs expected intertidal species, every design on a grid; Pareto frontiers; best per budget; cheapest per target.
#   Pt. Conception: V visits; bags at GP / CJ / JB (0 = site skipped); f filters per bag; k PCRs per filter per marker (one pooled library).
#   Coastline:      N exchangeable locations (one visit, one site each); b bags; f filters; k PCRs.
# Cost (costs.csv): per visit to a location travel_location + travel_site per extra site + per bag (sample + f x (filter + sum over markers (library + pcr x k))).
# Run for each extract fraction a and filter type. Env: A (one value or all), FTYPE (disc | capsule).
source("common.R"); source("model.R")
fits <- readRDS(paste0(FITDIR, Sys.getenv("FITS", "fits.rds"))); C <- read_costs()   # env FITS: e.g. fits_a0.0833.rds
FTYPE <- Sys.getenv("FTYPE", "disc"); AS <- if (nzchar(Sys.getenv("A"))) as.numeric(Sys.getenv("A")) else A_GRID
cfilter <- C[[paste0("filter_", FTYPE)]]
KG <- c(1, 2, 3, 4, 6, 8, 10, 12, 15, 20, 25, 30); FG <- 1:3; VG <- if (nzchar(Sys.getenv("VG"))) as.numeric(strsplit(Sys.getenv("VG"), ",")[[1]]) else c(1:6, 8, 10, 12, 16)   # env VG: visit grid (was 1,2,3,4,6 before 2026-10-11; widened when hour costs pushed the optimum to the edge)
BUDGETS <- c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000)
MAXB <- c(GP = 12, CJ = 6, JB = 6); MAXCOST <- max(BUDGETS)
frontier <- function(cost, S) { o <- order(cost, -S); keep <- o[S[o] > cummax(c(-Inf, S[o][-length(o)]))]; keep }   # indices on the Pareto front
summ_pc <- list(); summ_coast <- list(); fr_all <- list(); tgt <- list()
for (a in AS) {
  f12 <- fits[[paste("12S", a)]]; fco <- fits[[paste("COI", a)]]; FF <- list(`12S` = f12, COI = fco)
  kc <- setNames(c(f12$kap_capsule, fco$kap_capsule), MARKERS)
  shared <- intersect(f12$taxa, fco$taxa); lg("a", a, FTYPE, ": shared intertidal species names", length(shared))
  # ---- per marker, per (f, k): log P(site-visit not detected) per cell, for b bags at site s  ->  T[[m]][[s]][b+1, cell] for each (f,k)
  fk <- expand.grid(f = FG, k = KG)
  # objective: expected number of the OBSERVED species detected (weights = posterior cell mass of each observed species).
  # The unseen-species mass (nu) is weakly identified (drifts with EM iterations at flat likelihood), so S_all is secondary.
  W <- lapply(FF, function(f) colSums(f$post)); NU <- lapply(FF, function(f) f$cells$nu)
  pre <- lapply(MARKERS, function(m) { ce <- FF[[m]]$cells
    lapply(seq_len(nrow(fk)), function(j) { pb <- ce$q * (1 - (1 - filt_det(ce, a, fk$k[j], FTYPE, kc[[m]]))^fk$f[j])
      c(setNames(lapply(c("GP", "CJ", "JB"), function(s) t(sapply(0:MAXB[[s]], function(b) log(pmax(1 - ce[[paste0("phi_", s)]] * (1 - (1 - pb)^b), 1e-300))))), c("GP", "CJ", "JB")),
        list(pb = pb)) }) }); names(pre) <- MARKERS
  # ---- Pt. Conception grid
  lay <- expand.grid(V = VG, GP = 0:MAXB[["GP"]], CJ = 0:MAXB[["CJ"]], JB = 0:MAXB[["JB"]], f = FG); lay <- lay[lay$GP + lay$CJ + lay$JB > 0, ]
  lay$nsite <- (lay$GP > 0) + (lay$CJ > 0) + (lay$JB > 0); lay$nb <- lay$GP + lay$CJ + lay$JB
  lay$fixed <- with(lay, V * (C[["travel_location"]] + C[["travel_site"]] * (nsite - 1) + nb * (C[["sample"]] + f * cfilter)))
  lay$nlib <- with(lay, V * nb * f)
  lay <- lay[lay$fixed + lay$nlib * (C[["library"]] + C[["pcr"]]) <= MAXCOST, ]
  lg("Pt. Conception layouts:", nrow(lay))
  pts <- list()
  for (i in seq_len(nrow(lay))) { d <- lay[i, ]; jj <- which(fk$f == d$f)
    Sm <- list(); Sa <- list(); Dm <- list()
    for (m in MARKERS) { ce <- FF[[m]]$cells
      P <- sapply(jj, function(j) { T <- pre[[m]][[j]]; 1 - exp(d$V * (T$GP[d$GP + 1, ] + T$CJ[d$CJ + 1, ] + T$JB[d$JB + 1, ])) })   # cells x k
      Sm[[m]] <- colSums(W[[m]] * P); Sa[[m]] <- colSums(NU[[m]] * P); Dm[[m]] <- (FF[[m]]$post %*% P)[match(shared, FF[[m]]$taxa), , drop = FALSE] }
    cl <- C[["library"]] + C[["pcr"]] * KG
    pts[[length(pts) + 1]] <- data.frame(i = i, markers = "12S", k12 = KG, kCOI = NA, S = Sm$`12S`, S_all = Sa$`12S`, cost = d$fixed + d$nlib * cl)
    pts[[length(pts) + 1]] <- data.frame(i = i, markers = "COI", k12 = NA, kCOI = KG, S = Sm$COI, S_all = Sa$COI, cost = d$fixed + d$nlib * cl)
    OV <- crossprod(Dm$`12S`, Dm$COI); S2 <- outer(Sm$`12S`, Sm$COI, "+") - OV; S2a <- outer(Sa$`12S`, Sa$COI, "+") - OV
    pts[[length(pts) + 1]] <- data.frame(i = i, markers = "12S+COI", k12 = rep(KG, length(KG)), kCOI = rep(KG, each = length(KG)), S = as.vector(S2), S_all = as.vector(S2a),
                                         cost = d$fixed + d$nlib * as.vector(outer(cl, cl, "+"))) }
  P <- bind_rows(pts); P <- cbind(P, lay[P$i, c("V", "GP", "CJ", "JB", "f")]); P <- P[P$cost <= MAXCOST, ]; P$a <- a; P$ftype <- FTYPE
  P$sites <- with(P, paste0(ifelse(GP > 0, "GP", ""), ifelse(CJ > 0, "+CJ", ""), ifelse(JB > 0, "+JB", ""))); P$sites <- sub("^\\+", "", P$sites)
  saveRDS(P, paste0(OUT, sprintf("pc_points_a%s_%s.rds", a, FTYPE)))
  for (mk in unique(P$markers)) { X <- P[P$markers == mk, ]
    fr <- X[frontier(X$cost, X$S), ]; fr_all[[paste("PC", a, mk)]] <- fr %>% mutate(scope = "Pt. Conception")
    for (Bd in BUDGETS) { Y <- X[X$cost <= Bd, ]; if (!nrow(Y)) next; b <- Y[which.max(Y$S), ]
      bs <- function(f) { z <- Y$S[f]; if (length(z)) max(z) else NA }
      summ_pc[[paste(a, mk, Bd)]] <- data.frame(scope = "Pt. Conception", a = a, ftype = FTYPE, markers = mk, budget = Bd, b[, c("V", "GP", "CJ", "JB", "f", "k12", "kCOI", "S", "S_all", "cost")],
        S_GP_only = bs(Y$sites == "GP"), S_GP_CJ = bs(Y$sites == "GP+CJ"), S_GP_JB = bs(Y$sites == "GP+JB"), S_all3 = bs(Y$sites == "GP+CJ+JB"),
        S_CJ_only = bs(Y$sites == "CJ"), S_one_visit = bs(Y$V == 1), S_one_filter = bs(Y$f == 1)) }
    for (St in c(10, 20, 30, 40, 50, 75, 100, 125, 150, 175, 200)) { Y <- fr[fr$S >= St, ]; if (!nrow(Y)) next; b <- Y[which.min(Y$cost), ]
      tgt[[paste("PC", a, mk, St)]] <- data.frame(scope = "Pt. Conception", a = a, ftype = FTYPE, markers = mk, target_S = St, b[, c("V", "GP", "CJ", "JB", "f", "k12", "kCOI", "S", "cost")]) } }
  # ---- coastline: N locations x b bags x f filters x k
  cg <- expand.grid(N = 1:40, b = 1:15, f = FG); cg$fixed <- with(cg, N * (C[["travel_location"]] + b * (C[["sample"]] + f * cfilter))); cg$nlib <- with(cg, N * b * f)
  cg <- cg[cg$fixed + cg$nlib * (C[["library"]] + C[["pcr"]]) <= 4 * MAXCOST, ]
  cpts <- list()
  for (m in MARKERS) { ce <- FF[[m]]$cells
    for (j in seq_len(nrow(fk))) { pb <- pre[[m]][[j]]$pb; g <- cg[cg$f == fk$f[j], ]
      lnot <- sapply(1:15, function(b) log(pmax(1 - ce$phi_COAST * (1 - (1 - pb)^b), 1e-300)))      # cells x b
      w <- W[[m]]; S <- sapply(seq_len(nrow(g)), function(r) sum(w * (1 - exp(g$N[r] * lnot[, g$b[r]])))); Sall <- sapply(seq_len(nrow(g)), function(r) sum(ce$nu * (1 - exp(g$N[r] * lnot[, g$b[r]]))))
      cpts[[length(cpts) + 1]] <- data.frame(g[, c("N", "b", "f")], k = fk$k[j], markers = m, S = S, S_all = Sall, cost = g$fixed + g$nlib * (C[["library"]] + C[["pcr"]] * fk$k[j])) } }
  CP <- bind_rows(cpts); CP$a <- a; CP$ftype <- FTYPE; saveRDS(CP, paste0(OUT, sprintf("coast_points_a%s_%s.rds", a, FTYPE)))
  for (m in MARKERS) { X <- CP[CP$markers == m, ]; fr <- X[frontier(X$cost, X$S), ]; fr_all[[paste("COAST", a, m)]] <- fr %>% mutate(scope = "coastline")
    for (Bd in c(BUDGETS, 15000, 20000, 30000, 40000)) { Y <- X[X$cost <= Bd, ]; if (!nrow(Y)) next; b <- Y[which.max(Y$S), ]
      bs <- function(f) { z <- Y$S[f]; if (length(z)) max(z) else NA }
      summ_coast[[paste(a, m, Bd)]] <- data.frame(scope = "coastline", a = a, ftype = FTYPE, markers = m, budget = Bd, b[, c("N", "b", "f", "k", "S", "S_all", "cost")],
        S_N_le_5 = bs(Y$N <= 5), S_N_le_10 = bs(Y$N <= 10), S_N_le_21 = bs(Y$N <= 21), S_b10_f2_k6 = bs(Y$b == 10 & Y$f == 2 & Y$k == 6)) } }
}
TAG <- paste0(FTYPE, if (length(AS) == 1) paste0("_a", AS) else "")
wcsv(bind_rows(summ_pc), paste0("t_pc_best_by_budget_", TAG)); wcsv(bind_rows(tgt), paste0("t_pc_min_cost_by_target_", TAG))
wcsv(bind_rows(summ_coast), paste0("t_coast_best_by_budget_", TAG)); wcsv(bind_rows(fr_all), paste0("t_frontiers_", TAG))
print(bind_rows(summ_pc) %>% filter(markers == "12S+COI"))
