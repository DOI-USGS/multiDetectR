# 22_single_reaction_detections.R -- peer (paper chat) 2026-10-11 evening: before building a tube-level template model, check the detections
# made by only ONE of a tube's separately sequenced reactions (j = 1). Are they low-read, and do blanks or index-hopping flags explain them?
# Runs with ~10 singles per tube: May 2024 (JVB3506; 12S block B ids 13-22, COI 1-11) and Jan 2024 12S (JVB3105). Target (intertidal)
# species-rank calls, reads > 0 = detected, as in the fit (01_build.R).
# Per tube x species x detecting library: reads of the species (summed over its ESVs), reads as a share of the library, j, hop flags of those
# ESV x library rows (archive hop_flag; blank = not assessed), the ESVs' contaminant_flag, and whether any of the species' ESVs DETECTED in
# that library also had reads in a blank of the same run (field, tap/RO or lab blank; from the workflow's blank-inclusive read tables,
# CalIntertidal_reads_long_hopflagged_<marker>.rds; blanks judged contaminated are included here) or in any blank at all.
# v2 (2026-10-11, after the CalIntertidal workflow session's check): overlap with a blank mostly measures FIELD DNA spilling into the blanks.
# Each shared ESV is weighted by its field:blank read ratio within the run (field reads summed over the run's non-blank libraries, blank reads
# over its blank libraries). A detection is a contamination SUSPECT only if one of its ESVs has a ratio <= 1 (10 as a looser bound).
source("common.R"); OUT <- FITDIR
B <- readRDS(paste0(FITDIR, "build.rds")); ESVD <- Sys.getenv("HOPDIR", paste0(ROOT, "hop_out/"))  # PORTED: was Kevin's Mac path; hop-flagged reads dir (see hop_wrapper_r12c.R)
# r2 (2026-10-13): HOPDIR = the hop-flagged reads rebuilt from the deposit (raw_inputs.zip + code zip); byte-identical to the workflow files
SETS <- list(c("12S", "JVB3506"), c("COI", "JVB3506"), c("12S", "JVB3105")); rows <- list()
for (m in unique(sapply(SETS, `[`, 1))) {
  rd <- arc(sprintf("%s_reads", m)); ev <- arc(sprintf("%s_ESV", m)) %>% select(esv_id, contaminant_flag)
  bl <- readRDS(paste0(ESVD, sprintf("CalIntertidal_reads_long_hopflagged_%s.rds", m))) %>% filter(is_blank | is_lab_blank, n_reads > 0) %>%
    distinct(esv_id = observation_id, batch, blank = Barcode, event_id, blank_reads = n_reads)
  fr <- readRDS(paste0(ESVD, sprintf("CalIntertidal_reads_long_hopflagged_%s.rds", m))) %>% filter(!is_blank, !is_lab_blank, n_reads > 0) %>%
    group_by(esv_id = observation_id, batch) %>% summarise(field_reads = sum(n_reads), .groups = "drop")
  ratio <- bl %>% group_by(esv_id, batch) %>% summarise(blank_reads = sum(blank_reads), .groups = "drop") %>% left_join(fr, by = c("esv_id", "batch")) %>%
    mutate(field_reads = coalesce(field_reads, 0), ratio = field_reads / blank_reads)
  for (st in SETS[sapply(SETS, `[`, 1) == m]) { bt <- st[2]
    L <- B$L %>% filter(marker == m, batch == bt, k == 1)
    r <- rd %>% filter(library_id %in% L$library_id, rank == "species", reads > 0, taxon %in% B$pres[[m]]$taxon) %>% left_join(ev, by = "esv_id")
    lib_reads <- setNames(L$reads_all_esvs, L$library_id)
    D <- r %>% group_by(Barcode, taxon, library_id) %>%
      summarise(reads = sum(reads), n_esvs = n_distinct(esv_id), hop_not_valid = any(nzchar(hop_flag) & hop_flag != "valid"),
                hop_not_assessed = all(!nzchar(hop_flag) | is.na(hop_flag)), contaminant_flagged = any(contaminant_flag != "no_control_evidence"),
                in_blank_same_run = any(esv_id %in% bl$esv_id[bl$batch == bt]), in_any_blank = any(esv_id %in% bl$esv_id),
                min_field_blank_ratio = suppressWarnings(min(ratio$ratio[ratio$batch == bt & ratio$esv_id %in% esv_id])), .groups = "drop") %>%
      group_by(Barcode, taxon) %>% mutate(j = n()) %>% ungroup() %>%
      mutate(n_singles = as.vector(table(L$Barcode)[Barcode]), share = reads / lib_reads[library_id], marker = m, batch = bt,
             j_class = cut(j, c(0, 1, 4, 100), labels = c("j = 1", "j = 2-4", "j >= 5")))
    rows[[paste(m, bt)]] <- D } }
D <- bind_rows(rows)
S <- D %>% group_by(marker, batch, j_class) %>% summarise(detections = n(), tube_species = n_distinct(paste(Barcode, taxon)),
  reads_q10 = quantile(reads, .1), reads_q25 = quantile(reads, .25), reads_median = median(reads), reads_q75 = quantile(reads, .75), reads_q90 = quantile(reads, .9),
  under_10_reads = mean(reads < 10), under_100_reads = mean(reads < 100), share_median = median(share),
  hop_not_valid = sum(hop_not_valid), hop_not_assessed = sum(hop_not_assessed), contaminant_flagged = sum(contaminant_flagged),
  in_blank_same_run = sum(in_blank_same_run), in_any_blank = sum(in_any_blank),
  suspect_ratio_le_1 = sum(min_field_blank_ratio <= 1), suspect_ratio_le_10 = sum(min_field_blank_ratio <= 10), .groups = "drop")
# species behind the j = 1 detections: are they species that elsewhere detect in many reactions?
SP <- D %>% group_by(marker, batch, taxon) %>% summarise(tubes = n_distinct(Barcode), tubes_j1 = n_distinct(Barcode[j == 1]),
  tubes_j_all = n_distinct(Barcode[j == n_singles]), median_reads_j1 = median(reads[j == 1]), median_reads_j_ge2 = median(reads[j >= 2]), .groups = "drop") %>%
  arrange(marker, batch, desc(tubes_j1))
wcsv(D, "t_single_reaction_detections_rows_v2"); wcsv(S, "t_single_reaction_detections_summary_v2"); wcsv(SP, "t_single_reaction_detections_by_species_v2")
print(D %>% filter(min_field_blank_ratio <= 10) %>% group_by(marker, batch, taxon) %>% summarise(detections = n(), min_ratio = min(min_field_blank_ratio), j = paste(sort(unique(j)), collapse = ","), .groups = "drop") %>% as.data.frame())
print(S %>% mutate(across(where(is.double), ~ signif(.x, 3))) %>% as.data.frame())
print(SP %>% group_by(marker, batch) %>% slice_head(n = 8) %>% as.data.frame())
