# 09_rule_vs_grid.R -- compare Kevin's analytical rule (paperA methods draft 7 + S1.5, Aug 2026) with exhaustive integer search,
# using the fitted Pt. Conception detection parameters. Government Point alone, one marker, four nested levels:
#   level 1 visit (theta1 = phi_GP, c1 = travel_location), 2 sample (theta2 = q, c2 = sample), 3 filter (theta3 = r, c3 = disc filter + library),
#   4 PCR (theta4 = p, c4 = pcr).  The rule assumes independent replicates, so both methods use a -> 0 (m_eff = k) on the a = 0.01 fit;
#   the finite-extract correction is shown separately (part D).
# A. single target (each fitted species type): rule (inward pass, integer argmin per level = rounding in place, then n1 = floor(B / C1))
#    vs exhaustive integer search for P1 at the same budget.
# B. assemblage (expected number of observed species, = the detection goal summed over targets): exhaustive search vs the rule applied to
#    one representative target (weighted-mean theta; weighted-median species type) and vs the paper's "effective single target".
# C. branches (GP, CJ, JB): S1.5 small-budget criterion (weighted average rate) vs the site order the grid picked (03).
# D. finite extract: the grid optimum under a = 0.01, 0.03, 0.05, 0.0833 vs a -> 0.
source("common.R")
fits <- readRDS(paste0(FITDIR, "fits.rds")); C <- read_costs()
cost4 <- c(visit = C[["travel_location"]], sample = C[["sample"]], filter = C[["filter_disc"]] + C[["library"]], pcr = C[["pcr"]])
BUD <- c(750, 1000, 1500, 2000, 3000, 5000, 7500, 10000)
NG <- list(V = 1:30, b = 1:15, f = 1:4, k = 1:60)
G <- expand.grid(V = NG$V, b = NG$b, f = NG$f, k = NG$k); G$cost <- with(G, V * (cost4[["visit"]] + b * (cost4[["sample"]] + f * (cost4[["filter"]] + k * cost4[["pcr"]]))))
G <- G[G$cost <= max(BUD), ]
# detection of a cascade for a vector of targets th (rows = targets; cols = theta1..4) and design n = (V, b, f, k)
P1 <- function(th, V, b, f, k, a = 0) { m <- if (a > 0) m_eff(k, a) else k
  P4 <- 1 - (1 - th[, 4])^m; P3 <- 1 - (1 - th[, 3] * P4)^f; P2 <- 1 - (1 - th[, 2] * P3)^b; 1 - (1 - th[, 1] * P2)^V }
# the rule: inward pass, integer n at each level minimising kappa_i = (c_i + n C_(i+1)) / h_i  (the paper's rounding in place),
# then n1 = floor(B / C1). Feasibility (paper: needed at tight budgets): if B buys no outermost unit, fix n1 = 1 and solve the chain
# below it with the remaining budget B - c1 (the paper's "n1 fixed by the question" case), recursively.
pass_in <- function(th, cst, nmax = 200) { d <- length(th); h <- -log(1 - th[d]); Cc <- cst[d]; n <- rep(NA, d)
  if (d > 1) for (i in (d - 1):1) { nn <- 1:nmax; hi <- -log(1 - th[i] * (1 - exp(-nn * h))); kap <- (cst[i] + nn * Cc) / hi; j <- which.min(kap)
    n[i + 1] <- nn[j]; Cc <- cst[i] + nn[j] * Cc; h <- hi[j] }
  list(n = n, C1 = Cc, kappa1 = Cc / h) }
rule <- function(th, cst, B) { th <- unname(th); cst <- unname(cst); d <- length(th); x <- pass_in(th, cst); n1 <- floor(B / x$C1)
  if (n1 >= 1 || d == 1) { n <- x$n; n[1] <- n1 } else { sub <- rule(th[-1], cst[-1], B - cst[1]); n <- c(1, unname(sub[seq_len(d - 1)])); if (any(n < 1)) n[] <- 0 }
  n <- unname(n); c(V = n[1], b = n[2], f = n[3], k = n[4], kappa1 = x$kappa1)[seq_len(if (d == 4) 5 else d)] }
# rule + outermost enumeration: for each affordable n1, give each level-1 unit B/n1 - c1 and solve the chain inside it by the rule
# (recursively, with the same feasibility step), keep the n1 with the best P1. One integer loop; everything inside is analytical.
Pchain <- function(th, n) { P <- 1; for (i in length(th):1) P <- 1 - (1 - th[i] * P)^n[i]; P }   # n = counts per containing unit
rule_inner <- function(th, cst, Bu) {                    # best counts inside one unit with budget Bu for its contents
  d <- length(th); x <- pass_in(th, cst); n1 <- floor(Bu / x$C1); if (n1 >= 1 || d == 1) { n <- x$n; n[1] <- max(n1, 0); return(n) }
  if (Bu < cst[1]) return(rep(0, d)); c(1, rule_inner(th[-1], cst[-1], Bu - cst[1])) }
rule_n1 <- function(th, cst, B) { th <- unname(th); cst <- unname(cst); best <- c(V = 0, b = 0, f = 0, k = 0); bp <- -1
  for (n1 in seq_len(floor(B / cst[1]))) { inner <- rule_inner(th[-1], cst[-1], B / n1 - cst[1]); if (any(inner < 1)) next
    p <- Pchain(th, c(n1, inner)); if (p > bp) { bp <- p; best <- c(V = n1, b = inner[1], f = inner[2], k = inner[3]) } }
  best }
fit <- fits[["COI 0.01"]]; res <- list()
for (mk in MARKERS) { fit <- fits[[paste(mk, 0.01)]]; ce <- fit$cells; W <- colSums(fit$post)
  ok <- ce$phi_GP > 0 & W > 1e-3; ce <- ce[ok, ]; W <- W[ok]; TH <- cbind(ce$phi_GP, ce$q, ce$r, ce$p)
  # ---- A: single target
  A <- bind_rows(lapply(seq_len(nrow(TH)), function(j) bind_rows(lapply(BUD, function(B) { g <- G[G$cost <= B, ]
    pg <- P1(TH[rep(j, nrow(g)), , drop = FALSE], g$V, g$b, g$f, g$k); best <- which.max(pg); r <- rule(TH[j, ], cost4, B)
    rn <- rule_n1(TH[j, ], cost4, B)
    pr <- if (r[["V"]] >= 1) P1(TH[j, , drop = FALSE], r[["V"]], r[["b"]], r[["f"]], r[["k"]]) else 0
    data.frame(marker = mk, type = j, weight = W[j], phi = TH[j, 1], q = TH[j, 2], r = TH[j, 3], p = TH[j, 4], budget = B,
               grid_design = paste(g$V[best], g$b[best], g$f[best], g$k[best], sep = "/"), P_grid = pg[best],
               rule_design = paste(r[["V"]], r[["b"]], r[["f"]], r[["k"]], sep = "/"), P_rule = pr, P_continuous = 1 - exp(-B / r[["kappa1"]]),
               rule_n1_design = paste(rn[["V"]], rn[["b"]], rn[["f"]], rn[["k"]], sep = "/"), P_rule_n1 = if (rn[["V"]] >= 1) P1(TH[j, , drop = FALSE], rn[["V"]], rn[["b"]], rn[["f"]], rn[["k"]]) else 0) }))))
  res[[paste("A", mk)]] <- A
  # ---- B: assemblage
  Bt <- bind_rows(lapply(BUD, function(B) { g <- G[G$cost <= B, ]
    S <- sapply(seq_len(nrow(g)), function(i) sum(W * P1(TH, g$V[i], g$b[i], g$f[i], g$k[i]))); best <- which.max(S)
    Sd <- function(r) if (r[["V"]] >= 1) sum(W * P1(TH, r[["V"]], r[["b"]], r[["f"]], r[["k"]])) else 0
    rm <- rule(colSums(TH * W) / sum(W), cost4, B)                                          # weighted-mean theta
    o <- order(-log(1 - TH[, 4]) * TH[, 3] * TH[, 2] * TH[, 1]); cw <- cumsum(W[o]) / sum(W); jm <- o[which(cw >= .5)[1]]
    rd <- rule(TH[jm, ], cost4, B)                                                         # weighted-median species type (by detectability)
    # effective single target: the species type whose own rule design scores best for the assemblage
    eff <- sapply(seq_len(nrow(TH)), function(j) Sd(rule(TH[j, ], cost4, B))); je <- which.max(eff); re <- rule(TH[je, ], cost4, B)
    eff2 <- sapply(seq_len(nrow(TH)), function(j) Sd(rule_n1(TH[j, ], cost4, B))); je2 <- which.max(eff2); re2 <- rule_n1(TH[je2, ], cost4, B)
    rm2 <- rule_n1(colSums(TH * W) / sum(W), cost4, B)
    data.frame(marker = mk, budget = B, grid_design = paste(g$V[best], g$b[best], g$f[best], g$k[best], sep = "/"), S_grid = S[best],
               rule_mean_design = paste(rm[["V"]], rm[["b"]], rm[["f"]], rm[["k"]], sep = "/"), S_rule_mean = Sd(rm),
               rule_median_design = paste(rd[["V"]], rd[["b"]], rd[["f"]], rd[["k"]], sep = "/"), S_rule_median = Sd(rd),
               rule_effective_design = paste(re[["V"]], re[["b"]], re[["f"]], re[["k"]], sep = "/"), S_rule_effective = eff[je],
               effective_theta = paste(round(TH[je, ], 3), collapse = "/"),
               rule_n1_mean_design = paste(rm2[["V"]], rm2[["b"]], rm2[["f"]], rm2[["k"]], sep = "/"), S_rule_n1_mean = Sd(rm2),
               rule_n1_effective_design = paste(re2[["V"]], re2[["b"]], re2[["f"]], re2[["k"]], sep = "/"), S_rule_n1_effective = eff2[je2],
               effective_theta_n1 = paste(round(TH[je2, ], 3), collapse = "/"), S_max_possible = sum(W)) }))
  res[[paste("B", mk)]] <- Bt
  # ---- D: finite extract (grid optimum for the assemblage at each a, using that a's fit)
  res[[paste("D", mk)]] <- bind_rows(lapply(c(0.01, 0.03, 0.05), function(a) { f2 <- fits[[paste(mk, a)]]; c2 <- f2$cells; w2 <- colSums(f2$post)
    k2 <- c2$phi_GP > 0 & w2 > 1e-3; T2 <- cbind(c2$phi_GP, c2$q, c2$r, c2$p)[k2, ]; w2 <- w2[k2]
    bind_rows(lapply(c(1500, 5000), function(B) { g <- G[G$cost <= B, ]
      S0 <- sapply(seq_len(nrow(g)), function(i) sum(w2 * P1(T2, g$V[i], g$b[i], g$f[i], g$k[i], a = 0)))
      Sa <- sapply(seq_len(nrow(g)), function(i) sum(w2 * P1(T2, g$V[i], g$b[i], g$f[i], g$k[i], a = a)))
      i0 <- which.max(S0); ia <- which.max(Sa)
      data.frame(marker = mk, a = a, budget = B, design_independent = paste(g$V[i0], g$b[i0], g$f[i0], g$k[i0], sep = "/"),
                 design_finite_extract = paste(g$V[ia], g$b[ia], g$f[ia], g$k[ia], sep = "/"), S_finite_at_independent_design = Sa[i0], S_finite_optimum = Sa[ia]) })) }))
  lg(mk, "done") }
A <- bind_rows(res[grep("^A", names(res))]); Bt <- bind_rows(res[grep("^B", names(res))]); D <- bind_rows(res[grep("^D", names(res))])
wcsv(A, "t_rule_vs_grid_single_target_v3"); wcsv(Bt, "t_rule_vs_grid_assemblage_v3"); if (!file.exists(paste0(OUT, "t_rule_vs_grid_finite_extract.csv"))) wcsv(D, "t_rule_vs_grid_finite_extract")
cat("A. single target: share of (type, budget) where rule = grid design; P loss (grid - rule)\n")
print(A %>% group_by(marker) %>% summarise(n = n(), same_design = mean(grid_design == rule_design), max_loss = max(P_grid - P_rule), wmean_loss = weighted.mean(P_grid - P_rule, weight),
                                           n_loss_gt_.01 = sum(P_grid - P_rule > .01), same_design_n1 = mean(grid_design == rule_n1_design), max_loss_n1 = max(P_grid - P_rule_n1),
                                           wmean_loss_n1 = weighted.mean(P_grid - P_rule_n1, weight), n_loss_n1_gt_.01 = sum(P_grid - P_rule_n1 > .01)) %>% as.data.frame())
print(A %>% mutate(loss_n1 = P_grid - P_rule_n1) %>% arrange(desc(loss_n1)) %>% select(marker, phi, q, r, p, budget, grid_design, P_grid, rule_n1_design, P_rule_n1) %>% head(8) %>% as.data.frame())
print(A %>% filter(P_grid - P_rule > .01) %>% arrange(desc(P_grid - P_rule)) %>% head(10) %>% as.data.frame())
cat("\nB. assemblage\n"); print(Bt %>% mutate(across(starts_with("S_"), ~ round(.x, 2))) %>% as.data.frame())
cat("\nD. finite extract\n"); print(D %>% mutate(across(starts_with("S_"), ~ round(.x, 2))) %>% as.data.frame())
