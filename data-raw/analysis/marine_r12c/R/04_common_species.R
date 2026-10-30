# 04_common_species.R -- (1) candidate "common" intertidal species: detected at the most Pt. Conception site-visits;
# (2) for a chosen set (target_species.csv: marker,taxon; default = top 5 per marker), the Pt. Conception design that maximises
#     the expected number of those species detected, per budget, and each species' detection probability under it.
source("common.R"); source("model.R")
fits <- readRDS(paste0(FITDIR, "fits.rds")); B <- readRDS(paste0(FITDIR, "build.rds")); C <- read_costs()
a <- as.numeric(Sys.getenv("A", "0.03")); FTYPE <- "disc"; cfilter <- C[["filter_disc"]]
cand <- bind_rows(lapply(MARKERS, function(m) { L <- B$L[B$L$marker == m, ]; pr <- B$pres[[m]] %>% left_join(L %>% select(library_id, site, visit_id, location), by = "library_id")
  units_pc <- L %>% filter(site != "COAST") %>% distinct(site, visit_id) %>% nrow(); n_loc <- n_distinct(L$location[L$site == "COAST"])
  pr %>% group_by(taxon) %>% summarise(marker = m, pc_site_visits = n_distinct(paste(site, visit_id)[site != "COAST"]), pc_site_visits_total = units_pc,
                                      GP_visits = n_distinct(visit_id[site == "GP"]), coast_locations = n_distinct(location[site == "COAST"]), coast_locations_total = n_loc, .groups = "drop") })) %>%
  arrange(marker, desc(pc_site_visits), desc(coast_locations)) %>% group_by(marker) %>% mutate(rank = row_number()) %>% ungroup()
wcsv(cand, "t_common_species_candidates"); print(as.data.frame(cand %>% filter(rank <= 12)))
tf <- paste0(ROOT, "target_species.csv")
T <- if (file.exists(tf)) utils::read.csv(tf, stringsAsFactors = FALSE) else cand %>% filter(rank <= 5) %>% select(marker, taxon)
lg("target set:", paste0(T$taxon, " [", T$marker, "]", collapse = "; "))
FF <- list(`12S` = fits[[paste("12S", a)]], COI = fits[[paste("COI", a)]]); KG <- c(1, 2, 3, 4, 6, 8, 10, 12, 15, 20, 25, 30); VG <- c(1, 2, 3, 4, 6)
lay <- expand.grid(V = VG, GP = 0:12, CJ = 0:6, JB = 0:6, f = 1:3); lay <- lay[lay$GP + lay$CJ + lay$JB > 0, ]
lay$fixed <- with(lay, V * (C[["travel_location"]] + C[["travel_site"]] * ((GP > 0) + (CJ > 0) + (JB > 0) - 1) + (GP + CJ + JB) * (C[["sample"]] + f * cfilter)))
lay$nlib <- with(lay, V * (GP + CJ + JB) * f); lay <- lay[lay$fixed + lay$nlib * (C[["library"]] + C[["pcr"]]) <= 10000, ]
mk_used <- unique(T$marker); res <- list()
for (i in seq_len(nrow(lay))) { d <- lay[i, ]; D <- list()
  for (m in mk_used) { f <- FF[[m]]; tx <- T$taxon[T$marker == m]; ix <- match(tx, f$taxa)
    D[[m]] <- sapply(KG, function(k) { P <- det_pc(f$cells, a, d$V, c(GP = d$GP, CJ = d$CJ, JB = d$JB), d$f, k, kap_capsule = f$kap_capsule); as.vector(f$post[ix, , drop = FALSE] %*% P) }) }
  # expected number of target species detected; markers chosen independently (a species is counted under the marker it is listed with)
  kk <- expand.grid(lapply(mk_used, function(m) seq_along(KG))); names(kk) <- mk_used
  E <- rowSums(sapply(mk_used, function(m) colSums(matrix(D[[m]], ncol = length(KG)))[kk[[m]]]))
  cost <- d$fixed + d$nlib * rowSums(sapply(mk_used, function(m) C[["library"]] + C[["pcr"]] * KG[kk[[m]]]))
  res[[i]] <- data.frame(i = i, kk, E = E, cost = cost, check.names = FALSE) }
R <- bind_rows(res); R <- R[R$cost <= 10000, ]
out <- list()
for (Bd in c(500, 750, 1000, 1500, 2000, 3000, 5000, 10000)) { Y <- R[R$cost <= Bd, ]; if (!nrow(Y)) next; b <- Y[which.max(Y$E), ]; d <- lay[b$i, ]
  pr <- unlist(lapply(mk_used, function(m) { f <- FF[[m]]; tx <- T$taxon[T$marker == m]; k <- KG[b[[m]]]
    setNames(as.vector(f$post[match(tx, f$taxa), , drop = FALSE] %*% det_pc(f$cells, a, d$V, c(GP = d$GP, CJ = d$CJ, JB = d$JB), d$f, k, kap_capsule = f$kap_capsule)), tx) }))
  out[[as.character(Bd)]] <- data.frame(budget = Bd, a = a, d[, c("V", "GP", "CJ", "JB", "f")], k = paste(sapply(mk_used, function(m) paste0(m, ":", KG[b[[m]]])), collapse = " "),
                                        cost = b$cost, expected_targets_detected = b$E, n_targets = nrow(T), t(round(pr, 3)), check.names = FALSE) }
wcsv(bind_rows(out), sprintf("t_common_species_designs_a%s", a)); print(bind_rows(out))
