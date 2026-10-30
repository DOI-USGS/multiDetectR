# 13_extract_cap33.R -- Jonah Ventures (email relayed by Kevin Lafferty, 2026-10-11): 100 uL of extract, 3 uL per high-sensitivity PCR,
# "at least 33 replicates. No multiplexing. So 1 assay 33 times, or 33 assays 1 time."  ->  a = 0.03 and at most 33 PCRs per filter's
# extract, summed over markers (k12 + kCOI <= 33). Best design per budget with and without the cap, from the saved a = 0.03 design points (03).
source("common.R")
CAP <- 33; a <- 0.03; BUD <- c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000); out <- list()
P <- readRDS(paste0(OUT, sprintf("pc_points_a%s_disc.rds", a))); P$ktot <- rowSums(cbind(P$k12, P$kCOI), na.rm = TRUE)
for (mk in c("12S", "COI", "12S+COI")) { X <- P[P$markers == mk, ]
  for (Bd in BUD) for (capped in c(FALSE, TRUE)) { Y <- X[X$cost <= Bd & (!capped | X$ktot <= CAP), ]; if (!nrow(Y)) next; b <- Y[which.max(Y$S), ]
    out[[length(out) + 1]] <- data.frame(a = a, cap = if (capped) CAP else NA, markers = mk, budget = Bd, b[, c("V", "sites", "GP", "CJ", "JB", "f", "k12", "kCOI", "cost", "S")]) } }
R <- bind_rows(out); wcsv(R, "t_extract_cap33_designs")
R$design <- with(R, paste0("V", V, " ", GP, "/", CJ, "/", JB, " f", f, " k", ifelse(is.na(k12), "", k12), "|", ifelse(is.na(kCOI), "", kCOI)))
print(R %>% mutate(S = round(S, 1), capped = !is.na(cap)) %>% select(markers, budget, capped, design, S) %>%
        tidyr::pivot_wider(names_from = capped, values_from = c(design, S)) %>% as.data.frame(), right = FALSE)
