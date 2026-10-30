# 27_rule_check_fine.R -- check of the rule designs (25) against an exhaustive search on EVERY integer design near them. 03's grid skips PCR counts
# (1, 2, 3, 4, 6, 8, 10, 12, 15, ...) and visit counts above 6, so the rule can beat it; this search cannot be beaten by the rule inside its box.
# Box: V 1..max(rule V) + 4; bags GP 0..12, CJ 0..6, JB 0..6; f 1..3; every integer k per marker from 1 to max(26, largest rule k + 6) (v2; v1 used rule k +- 6).
# Exact expectations as in 03 and 25; only the best design per budget is kept. Env FITS as in 25 (reads t_rule_designs<SFX>.csv).
source("common.R"); source("model.R")
FN <- Sys.getenv("FITS", "fits.rds"); fits <- readRDS(paste0(FITDIR, FN)); C <- read_costs(); a <- 0.03
SFX <- if (FN == "fits.rds") "" else paste0("_", sub("^fits_(.*)\\.rds$", "\\1", FN))
RV <- Sys.getenv("RULEVER", "_v2"); RD <- read.csv(paste0(OUT, "t_rule_designs", SFX, RV, ".csv"), stringsAsFactors = FALSE)
FF <- list(`12S` = fits[[paste("12S", a)]], COI = fits[[paste("COI", a)]]); W <- lapply(FF, function(f) colSums(f$post)); shared <- intersect(FF$`12S`$taxa, FF$COI$taxa)
cf <- C[["filter_disc"]]; cl <- C[["library"]]; cp <- C[["pcr"]]; cs <- C[["sample"]]; MAXB <- c(GP = 12, CJ = 6, JB = 6); BUD <- sort(unique(RD$budget))
KR <- lapply(c(`12S` = "k12", COI = "kCOI"), function(col) { z <- RD[[col]][!is.na(RD[[col]])]; 1:max(26, max(z) + 6) })   # v2: every k from 1
VMAX <- max(RD$V) + 4; FG <- 1:3
# per marker, f, k: log P(site-visit not detected) for b bags at each site (cells x (b+1))
T <- lapply(names(FF), function(m) { ce <- FF[[m]]$cells; lapply(FG, function(f) lapply(KR[[m]], function(k) { pb <- ce$q * (1 - (1 - filt_det(ce, a, k, "disc", FF[[m]]$kap_capsule))^f)
  lapply(c(GP = "GP", CJ = "CJ", JB = "JB"), function(s) sapply(0:MAXB[[s]], function(b) log(pmax(1 - ce[[paste0("phi_", s)]] * (1 - (1 - pb)^b), 1e-300)))) })) }); names(T) <- names(FF)
lay <- expand.grid(GP = 0:12, CJ = 0:6, JB = 0:6, f = FG); lay <- lay[lay$GP + lay$CJ + lay$JB > 0, ]; lay$ns <- (lay$GP > 0) + (lay$CJ > 0) + (lay$JB > 0); lay$nb <- lay$GP + lay$CJ + lay$JB
best <- list()
for (mk in c("12S", "COI", "12S+COI")) { ms <- strsplit(mk, "\\+")[[1]]; t0 <- Sys.time()
  bestS <- setNames(rep(-Inf, length(BUD)), BUD); bestD <- vector("list", length(BUD))
  for (i in seq_len(nrow(lay))) { d <- lay[i, ]; fixed1 <- C[["travel_location"]] + C[["travel_site"]] * (d$ns - 1) + d$nb * (cs + d$f * cf)
    # per marker: matrix cells x k of log not-detected per visit
    LN <- lapply(ms, function(m) sapply(seq_along(KR[[m]]), function(j) { X <- T[[m]][[d$f]][[j]]; X$GP[, d$GP + 1] + X$CJ[, d$CJ + 1] + X$JB[, d$JB + 1] })); names(LN) <- ms
    for (V in seq_len(VMAX)) { vc <- function(kc) V * (fixed1 + d$nb * d$f * kc)                                     # kc = sum over markers of (library + pcr k)
      if (vc(length(ms) * (cl + cp)) > max(BUD)) break
      Sm <- list(); Dm <- list()
      for (m in ms) { P <- 1 - exp(V * LN[[m]]); Sm[[m]] <- colSums(W[[m]] * P); if (length(ms) == 2) Dm[[m]] <- (FF[[m]]$post %*% P)[match(shared, FF[[m]]$taxa), , drop = FALSE] }
      if (length(ms) == 1) { S <- Sm[[1]]; cost <- vc(cl + cp * KR[[ms]]); k1 <- KR[[ms]]; k2 <- NA }
      else { S <- as.vector(outer(Sm$`12S`, Sm$COI, "+") - crossprod(Dm$`12S`, Dm$COI)); k1 <- rep(KR$`12S`, length(KR$COI)); k2 <- rep(KR$COI, each = length(KR$`12S`))
        cost <- vc(2 * cl + cp * (k1 + k2)); ok <- k1 + k2 <= 33; S[!ok] <- -Inf }
      for (bi in seq_along(BUD)) { z <- which(cost <= BUD[bi]); if (!length(z)) next; j <- z[which.max(S[z])]
        if (S[j] > bestS[bi]) { bestS[bi] <- S[j]; bestD[[bi]] <- data.frame(V = V, GP = d$GP, CJ = d$CJ, JB = d$JB, f = d$f, k_a = k1[j], k_b = k2[j], cost = cost[j]) } } } }
  lg(mk, "searched in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min")
  for (bi in seq_along(BUD)) { b <- bestD[[bi]]; kk12 <- if ("12S" %in% ms) b$k_a else NA; kkC <- if (mk == "COI") b$k_a else if (mk == "12S+COI") b$k_b else NA
    best[[paste(mk, BUD[bi])]] <- data.frame(markers = mk, budget = BUD[bi], fine_design = sprintf("V%d %d/%d/%d f%d k%s|%s", b$V, b$GP, b$CJ, b$JB, b$f, ifelse(is.na(kk12), "", kk12), ifelse(is.na(kkC), "", kkC)),
      S_fine = bestS[bi], cost_fine = b$cost) } }
X <- RD %>% left_join(bind_rows(best), by = c("markers", "budget")) %>% mutate(loss_vs_fine = S_fine - S_rule, same_as_fine = rule_design == fine_design)
wcsv(X, paste0("t_rule_check_fine", SFX, RV))
print(X %>% select(markers, budget, rule_design, S_rule, fine_design, S_fine, loss_vs_fine, grid_design, S_grid) %>% mutate(across(where(is.double), ~ round(.x, 3))) %>% as.data.frame(), right = FALSE)
print(X %>% group_by(markers) %>% summarise(same = sum(same_as_fine), n = n(), max_loss = max(loss_vs_fine), mean_loss = mean(loss_vs_fine)) %>% as.data.frame())
