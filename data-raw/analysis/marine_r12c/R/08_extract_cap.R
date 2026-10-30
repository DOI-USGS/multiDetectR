# 08_extract_cap.R -- Jonah Ventures (2026-10-10): an archived extract "can typically be tested about 12 times".
# Two readings, each with its own extract fraction a and cap on PCRs drawn from ONE filter's extract (both markers together):
#   test = one marker library of 6 pooled PCRs  ->  a = 1/72 (~0.014; nearest fit a = 0.01); cap 72 PCRs per filter (never binding on the grid)
#   test = one PCR                               ->  a = 1/12 (fit a = 0.0833);              cap 12 PCRs per filter (k12 + kCOI <= 12)
# Best design per budget under each reading, from the saved design points (03). Disc filters.
source("common.R")
RD <- data.frame(reading = c("test = 6-PCR library", "test = one PCR"), a = c(0.01, 0.0833), cap = c(72, 12))
BUD <- c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000); out <- list()
for (j in seq_len(nrow(RD))) { rd <- RD[j, ]
  P <- readRDS(paste0(OUT, sprintf("pc_points_a%s_disc.rds", rd$a))); CP <- readRDS(paste0(OUT, sprintf("coast_points_a%s_disc.rds", rd$a)))
  P$ktot <- rowSums(cbind(P$k12, P$kCOI), na.rm = TRUE)
  for (mk in c("12S", "COI", "12S+COI")) { X <- P[P$markers == mk, ]
    for (Bd in BUD) for (capped in c(FALSE, TRUE)) { Y <- X[X$cost <= Bd & (!capped | X$ktot <= rd$cap), ]; if (!nrow(Y)) next; b <- Y[which.max(Y$S), ]
      out[[length(out) + 1]] <- data.frame(reading = rd$reading, a = rd$a, cap = if (capped) rd$cap else NA, scope = "Pt. Conception", markers = mk, budget = Bd,
        b[, c("V", "sites", "GP", "CJ", "JB", "f", "k12", "kCOI", "cost", "S")]) } }
  for (m in MARKERS) { X <- CP[CP$markers == m, ]
    for (Bd in c(BUD, 20000, 40000)) for (capped in c(FALSE, TRUE)) { Y <- X[X$cost <= Bd & (!capped | X$k <= rd$cap), ]; if (!nrow(Y)) next; b <- Y[which.max(Y$S), ]
      out[[length(out) + 1]] <- data.frame(reading = rd$reading, a = rd$a, cap = if (capped) rd$cap else NA, scope = "coastline", markers = m, budget = Bd,
        N = b$N, b_loc = b$b, f = b$f, k = b$k, cost = b$cost, S = b$S) } } }
R <- bind_rows(out); wcsv(R, "t_extract_cap_designs")
R$design <- with(R, ifelse(scope == "coastline", paste0("N", N, " b", b_loc, " f", f, " k", k),
                           paste0("V", V, " ", sites, " ", GP, "/", CJ, "/", JB, " f", f, " k", ifelse(is.na(k12), "", k12), "|", ifelse(is.na(kCOI), "", kCOI))))
print(R %>% mutate(S = round(S, 1), capped = !is.na(cap)) %>% select(reading, scope, markers, budget, capped, design, S) %>%
        tidyr::pivot_wider(names_from = capped, values_from = c(design, S)) %>% as.data.frame(), right = FALSE)
