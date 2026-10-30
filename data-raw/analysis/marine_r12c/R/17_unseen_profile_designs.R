# 17_unseen_profile_designs.R -- designs at the ends of the profile-likelihood interval for the number of unseen species U (15).
# 03 was run on fits_unseen_lo.rds / fits_unseen_hi.rds (env FITS, COSTTAG = hours_2026-10-11_unseen_lo / _hi). For each end, budget and
# marker set: the design that maximises expected species INCLUDING unseen (S_all) at that end, against the main design (best for observed
# species under the main fit), both scored by S_all at that end. Same layouts and costs, so design points match by (i, markers, k12, kCOI).
source("common.R")
a <- 0.03; BUD <- c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000)
key <- function(d) paste(d$i, d$markers, d$k12, d$kCOI); lab <- function(d) sprintf("V%d %d/%d/%d f%d k%s|%s", d$V, d$GP, d$CJ, d$JB, d$f, ifelse(is.na(d$k12), "", d$k12), ifelse(is.na(d$kCOI), "", d$kCOI))
M <- readRDS(paste0(OUT, sprintf("pc_points_a%s_disc.rds", a))); out <- list()
for (end in c("lo", "hi")) { E <- readRDS(sub("/$", sprintf("_unseen_%s/", end), OUT) %>% paste0(sprintf("pc_points_a%s_disc.rds", a)))
  stopifnot(identical(key(E), key(M)), isTRUE(all.equal(E$cost, M$cost)))
  for (mk in c("12S", "COI", "12S+COI")) for (Bd in BUD) { ix <- which(M$markers == mk & M$cost <= Bd); if (!length(ix)) next
    m0 <- ix[which.max(M$S[ix])]; e0 <- ix[which.max(E$S_all[ix])]
    out[[length(out) + 1]] <- data.frame(end = end, markers = mk, budget = Bd, main_design = lab(M[m0, ]), end_design = lab(E[e0, ]), same = m0 == e0,
      S_all_end_design = E$S_all[e0], S_all_main_design = E$S_all[m0], loss_species_using_main = E$S_all[e0] - E$S_all[m0],
      S_observed_main = M$S[m0], S_observed_end_design = M$S[e0]) } }
R <- bind_rows(out); wcsv(R, "t_unseen_profile_designs")
print(R %>% group_by(end, markers) %>% summarise(same = sum(same), budgets = n(), max_loss = round(max(loss_species_using_main), 2), .groups = "drop") %>% as.data.frame())
print(R %>% filter(!same) %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% select(end, markers, budget, main_design, end_design, loss_species_using_main), right = FALSE)
