# 16_pooling_penalty.R -- peer (paper chat) 2026-10-11: within-run evidence on pooling vs summing replicates.
# 12S, May 2024 Government Point (JVB3506, block B): per tube, ten separately sequenced PCR replicates (ids 13-22) and two pools (ids 23, 24),
# each pool built from 3 uL aliquots of those same ten PCR products (Schulte et al. 2026), all in one run. Pool membership is inferred
# (archive README: from ESVs each pool shares only with its block); id 12 is absent; the crosswalk awaits confirmation (Nick).
# The pools mix AMPLIFIED product after PCR, so this measures dilution in the library, not pooling of template before PCR.
# Per tube x species (species-rank calls; reads > 0 = detected), with j = singles detecting it out of n present:
#   recovery  = P(pool detects | at least one single does)  (the "sum" is the union of the singles)
#   M_equiv   = how many of its own source singles the pool is worth: the M at which the union of M singles drawn at random from the
#               tube's n (hypergeometric: 1 - choose(n - j, M) / choose(n, M)) detects as often as the pool, over pairs with j >= 1.
#               Same conditioning as the pool rate, so the union of all n gives exactly 1; interpolated between integer M.
# Reads only the archive; writes to the fit folder (cost-free).
source("common.R"); OUT <- FITDIR
L <- arc("libraries") %>% filter(batch == "JVB3506", marker == "12S")
R <- arc("12S_reads") %>% filter(library_id %in% L$library_id, rank == "species", reads > 0) %>% distinct(library_id, Barcode, taxon)
tg <- arc("species_habitat"); tcol <- intersect(c("rocky_intertidal_class", "habitat_class"), names(tg))[1]
targets <- unique(tg$species[tg[[tcol]] %in% HABITAT_KEEP])
sing <- L %>% filter(n_pcr == 1, replicate %in% 13:22); pools <- L %>% filter(n_pcr == 0, replicate %in% 23:24)
tubes <- intersect(unique(sing$Barcode), unique(pools$Barcode))
lg("tubes with block-B singles and pools:", length(tubes), "; singles per tube:", paste(range(table(sing$Barcode[sing$Barcode %in% tubes])), collapse = "-"))
rows <- list()
for (b in tubes) { s <- sing$library_id[sing$Barcode == b]; n <- length(s); tx <- unique(R$taxon[R$Barcode == b])
  for (t in tx) { j <- sum(s %in% R$library_id[R$taxon == t])
    for (p in pools$library_id[pools$Barcode == b]) rows[[length(rows) + 1]] <- data.frame(Barcode = b, taxon = t, n_singles = n, j = j, pool = p,
      pool_detects = p %in% R$library_id[R$taxon == t]) } }
D <- bind_rows(rows); D$target <- D$taxon %in% targets
D <- D %>% left_join(L %>% select(pool = library_id, pool_reads = reads_all_esvs), by = "pool")
depth <- L %>% filter(Barcode %in% tubes) %>% mutate(kind = ifelse(n_pcr == 0, "pool", "single")) %>% group_by(kind) %>%
  summarise(libraries = n(), median_reads = median(reads_all_esvs), .groups = "drop")
meff <- function(d) { d <- d[d$j >= 1, ]; obs <- mean(d$pool_detects); Ms <- 1:max(d$n_singles)                     # a tube with n < M singles uses all n
  g <- sapply(Ms, function(M) { m <- pmin(M, d$n_singles); mean(1 - choose(d$n_singles - d$j, m) / choose(d$n_singles, m)) })
  if (obs <= g[1]) return(1); approx(g, Ms, xout = obs)$y }
summ <- function(d, set) data.frame(set = set, tube_species_pool_rows = nrow(d), union_detections = sum(d$j >= 1),
  pool_detections = sum(d$pool_detects), pool_only = sum(d$pool_detects & d$j == 0),
  recovery = mean(d$pool_detects[d$j >= 1]), mean_single_rate = mean(d$j[d$j >= 1] / d$n_singles[d$j >= 1]),
  M_equiv = meff(d), n_singles_typical = median(d$n_singles))
S <- bind_rows(summ(D, "all species-rank calls"), summ(D[D$target, ], "intertidal targets"))
byj <- D %>% filter(j >= 1) %>% group_by(j) %>% summarise(rows = n(), pool_detect_rate = mean(pool_detects), .groups = "drop")
agree <- D %>% group_by(Barcode, taxon) %>% filter(n() == 2) %>% summarise(both = all(pool_detects), one = sum(pool_detects) == 1, j = first(j), .groups = "drop") %>%
  filter(both | one) %>% summarise(pairs_detected_by_a_pool = n(), detected_by_both_pools = sum(both), detected_by_one = sum(one))
wcsv(D, "t_pooling_penalty_rows"); wcsv(S, "t_pooling_penalty_summary"); wcsv(byj, "t_pooling_penalty_by_j")
print(depth); print(S); print(byj, n = 20); print(agree)
