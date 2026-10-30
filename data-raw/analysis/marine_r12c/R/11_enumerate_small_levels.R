# 11_enumerate_small_levels.R -- test the paper-helper's proposed integer procedure (2026-10-10):
# "where the relaxed optimum at any level is below about three units, enumerate that level and the outermost count, and apply the rule to the rest".
# Single targets at Government Point (same setup as 09: a -> 0, a = 0.01 fit, levels visit/sample/filter/PCR); compared with exhaustive search.
# Rule levels spend the whole per-unit budget: n_i = floor(Bu / C_i) at the rule's unit cost, and the leftover passes inward.
source("common.R")
fits <- readRDS(paste0(FITDIR, "fits.rds")); C <- read_costs()
cost4 <- unname(c(C[["travel_location"]], C[["sample"]], C[["filter_disc"]] + C[["library"]], C[["pcr"]]))
BUD <- c(750, 1000, 1500, 2000, 3000, 5000, 7500, 10000)
# env GRID=wide: exhaustive search over a wider grid (peer 2026-10-11: widen until no procedure design lies outside it)
# env LEFTOVER=unspent: a level that is not enumerated takes the relaxed count rounded in place (capped by what the budget affords)
#   and leaves the remainder unspent, instead of floor(Bu / C_i) with the leftover passed inward (the default "fill").
WIDE <- Sys.getenv("GRID") == "wide"; UNSPENT <- Sys.getenv("LEFTOVER") == "unspent"
GR <- if (WIDE) list(V = 1:60, b = 1:100, f = 1:12, k = 1:200) else list(V = 1:30, b = 1:15, f = 1:4, k = 1:60)
G <- expand.grid(V = GR$V, b = GR$b, f = GR$f, k = GR$k); G$cost <- with(G, V * (cost4[1] + b * (cost4[2] + f * (cost4[3] + k * cost4[4])))); G <- G[G$cost <= max(BUD), ]
P1 <- function(th, V, b, f, k) { P4 <- 1 - (1 - th[, 4])^k; P3 <- 1 - (1 - th[, 3] * P4)^f; P2 <- 1 - (1 - th[, 2] * P3)^b; 1 - (1 - th[, 1] * P2)^V }
Pch <- function(th, n) { P <- 1; for (i in length(th):1) P <- 1 - (1 - th[i] * P)^n[i]; P }
# relaxed (continuous) pass: u_i per level, returns relaxed counts n_i = u/h and unit costs C_i
relaxed <- function(th, cst) { d <- length(th); h <- -log(1 - th[d]); Cc <- cst[d]; n <- rep(NA, d); Cs <- rep(NA, d); Cs[d] <- Cc
  if (d > 1) for (i in (d - 1):1) { nn <- seq(1, 300, by = 0.05); hi <- -log(1 - th[i] * (1 - exp(-nn * h))); kap <- (cst[i] + nn * Cc) / hi; j <- which.min(kap)
    n[i + 1] <- nn[j]; Cc <- cst[i] + nn[j] * Cc; h <- hi[j]; Cs[i] <- Cc }
  list(n = n, C = Cs) }
# solve levels i..d inside one containing unit with budget Bu; levels in E are enumerated, the rest take floor(Bu / C_i)
BUDGET_AWARE <- nzchar(Sys.getenv("BUDGET_AWARE")); TRIG <- as.numeric(Sys.getenv("TRIGGER", "3"))   # enumerate levels with fewer than TRIG units   # also enumerate a level whose count AT THIS BUDGET is below 3
solve_from <- function(th, cst, Bu, E, i, Cs) { d <- length(th)
  if (BUDGET_AWARE && i <= d && floor(Bu / max(Cs[i], sum(cst[i:d]))) < TRIG) E <- union(E, i)
  if (i > d) return(list(n = integer(0), P = 1))
  cand <- if (i %in% E) seq_len(max(0, min(max(8, 2 * TRIG), floor(Bu / Cs[i] * 3), floor(Bu / cst[i])))) else if (UNSPENT) max(1, min(round(attr(Cs, "n")[i]), floor(Bu / sum(cst[i:d])))) else { floor(Bu / max(Cs[i], sum(cst[i:d]))) }
  best <- list(n = rep(0, d - i + 1), P = 0)
  for (m in cand[cand >= 1]) { inner <- solve_from(th, cst, Bu / m - cst[i], E, i + 1, Cs); if (any(inner$n < 1) && i < d) next
    P <- 1 - (1 - th[i] * inner$P)^m; if (P > best$P) best <- list(n = c(m, inner$n), P = P) }
  best }
proc <- function(th, cst, B) { rl <- relaxed(th, cst); E <- which(rl$n < TRIG); E <- E[!is.na(E)]; attr(rl$C, "n") <- rl$n
  best <- list(n = rep(0, 4), P = 0)
  for (n1 in seq_len(floor(B / cst[1]))) { inner <- solve_from(th, cst, B / n1 - cst[1], E, 2, rl$C); if (any(inner$n < 1)) next
    P <- 1 - (1 - th[1] * inner$P)^n1; if (P > best$P) best <- list(n = c(n1, inner$n), P = P) }
  c(best, list(E = paste(c("visit", "sample", "filter", "pcr")[c(1, E)], collapse = "+"))) }
res <- list()
for (mk in MARKERS) { f <- fits[[paste(mk, 0.01)]]; ce <- f$cells; W <- colSums(f$post); ok <- ce$phi_GP > 0 & W > 1e-3; ce <- ce[ok, ]; W <- W[ok]
  TH <- cbind(ce$phi_GP, ce$q, ce$r, ce$p)
  for (j in seq_len(nrow(TH))) for (B in BUD) { g <- G[G$cost <= B, ]; pg <- P1(TH[rep(j, nrow(g)), , drop = FALSE], g$V, g$b, g$f, g$k); i <- which.max(pg)
    pr <- proc(TH[j, ], cost4, B)
    res[[length(res) + 1]] <- data.frame(marker = mk, type = j, weight = W[j], budget = B, enumerated = pr$E, grid_design = paste(g$V[i], g$b[i], g$f[i], g$k[i], sep = "/"), P_grid = pg[i],
                                         proc_design = paste(pr$n, collapse = "/"), P_proc = pr$P,
                                         proc_outside_grid = pr$n[1] > max(GR$V) || pr$n[2] > max(GR$b) || pr$n[3] > max(GR$f) || pr$n[4] > max(GR$k)) } }
R <- bind_rows(res); wcsv(R, paste0("t_enumerate_small_levels", if (BUDGET_AWARE) "_budget_aware" else "", if (TRIG != 3) paste0("_trigger", TRIG) else "", if (UNSPENT) "_unspent" else "", if (WIDE) "_widegrid" else ""))
print(R %>% group_by(marker) %>% summarise(n = n(), same = mean(grid_design == proc_design), wmean_loss = weighted.mean(P_grid - P_proc, weight), max_loss = max(P_grid - P_proc),
                                           n_gt_.01 = sum(P_grid - P_proc > .01)) %>% as.data.frame())
print(R %>% count(enumerated)); print(R %>% arrange(desc(P_grid - P_proc)) %>% head(6))
