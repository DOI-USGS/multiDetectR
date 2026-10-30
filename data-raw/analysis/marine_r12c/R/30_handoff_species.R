# 30_handoff_species.R -- paper chat 2026-10-12: inputs for exact single-species and assemblage numbers on the TEMPLATE fit (a = 0.03).
# Writes to the handoff folder (never overwrites; identical content is left alone):
#  1. template_species_types_mixture: the fitted mixture, same columns as species_types_mixture.csv, plus template shape alpha and the run
#     multiplier used for designs (1: designs are in JVB6097 = current-protocol units).        per-batch multipliers: template_run_multipliers.csv already in the handoff (checked, not rewritten).
#  2. template_species_type_weights: posterior weight of each type for each observed species (rows of fit$post; sum to 1 per species x marker).
#  3. detection_histories_Hjenkinsi_Cviolaceus: one row per library used in the fit (both markers, every site), reads of each species (0 = none).
#  4. abalone_limpet_calls (archive calls at any rank in Haliotidae/Haliotis/Lottiidae/Lottia) + abalone_limpet_candidates (the three species as
#     plausible alternatives in ESV_site_calls) + abalone_limpet_COI_reference (TaxaLikely NCBI reference cache for the run).
#  5. effort_to_date_by_site (site x marker) and effort_to_date_by_filter (one row per filter: what an exact P(no detection) needs).
source("common.R"); source("model.R")
HO <- Sys.getenv("R12C_HANDOFF_OUT", paste0(ROOT, "handoff_out/"))  # PORTED: was Kevin's Mac path
hw <- function(x, name) { p <- paste0(HO, name, ".csv"); tf <- tempfile(fileext = ".csv"); utils::write.csv(x, tf, row.names = FALSE, na = "")
  if (file.exists(p)) { if (unname(tools::md5sum(p)) == unname(tools::md5sum(tf))) { lg("unchanged", name); return(invisible()) }; stop("exists and differs, not overwritten: ", p) }
  file.copy(tf, p); lg("wrote", name, nrow(x), "rows") }
a <- 0.03; fits <- readRDS(paste0(FITDIR, "fits_template.rds")); B <- readRDS(paste0(FITDIR, "build.rds")); L <- B$L
REFCACHE <- Sys.getenv("R12C_REFCACHE", paste0(ROOT, "cache_reference/"))  # PORTED: was Kevin's Mac path; skipped when SKIP_REFCACHE=1

# ---- 1. mixture
types <- bind_rows(lapply(MARKERS, function(m) { f <- fits[[paste(m, a)]]; ce <- f$cells
  data.frame(marker = m, a = a, type = seq_len(nrow(ce)), weight_observed = colSums(f$post), weight_total_nu = ce$nu, prob_never_detected = ce$P0,
             phi_GP = ce$phi_GP, phi_CJ = ce$phi_CJ, phi_JB = ce$phi_JB, phi_COAST = ce$phi_COAST, theta_bag = ce$q, theta_filter_disc = ce$r,
             theta_filter_capsule = odds_mult(ce$r, f$kap_capsule), theta_pcr = ce$p, template_alpha = f$tmpl$alpha, run_multiplier_design = 1) }))
hw(types, "template_species_types_mixture")
rmt <- bind_rows(lapply(MARKERS, function(m) { f <- fits[[paste(m, a)]]; b <- sort(unique(L$batch[L$marker == m]))
  data.frame(marker = m, a = a, batch = b, run_multiplier = ifelse(b %in% names(f$run_mult), f$run_mult[b], 1),
             platform = ifelse(b %in% names(PLATFORM), "illumina", "nanopore"), reference_batch = b == REF_BATCH) }))
# PORTED: 23_template_fit.R writes this as t_template_run_multipliers.csv to FITDIR (the TOP of
# out_r12c/, because 23 sets OUT <- FITDIR), not as template_run_multipliers.csv in the handoff
# folder (Kevin's handoff folder was his output dir, so the two coincided for him; in the port
# they differ). Read it from FITDIR; fall back to an HO copy if one is already present.
.trm <- if (file.exists(paste0(HO, "template_run_multipliers.csv"))) paste0(HO, "template_run_multipliers.csv") else paste0(FITDIR, "t_template_run_multipliers.csv")
ex <- utils::read.csv(.trm) %>% left_join(rmt, by = c("marker", "batch"))   # from 23_template_fit.R
stopifnot(nrow(ex) == nrow(rmt), all(abs(ex$run_multiplier_template - ex$run_multiplier) < 1e-12)); lg("template_run_multipliers.csv agrees with the fit")

# ---- 2. per-species posterior weights over types
w <- bind_rows(lapply(MARKERS, function(m) { P <- fits[[paste(m, a)]]$post
  stopifnot(all(abs(rowSums(P) - 1) < 1e-4))   # types with nu < 1e-7 of the total are pruned by fit_model: rows sum to 1 within 1e-5 on r2;
                                                # r12c COI (640 types) has one row at 1 - 1.7e-5, so the tolerance is 1e-4 since 2026-10-25
  data.frame(marker = m, a = a, taxon = rep(rownames(P), times = ncol(P)), type = rep(seq_len(ncol(P)), each = nrow(P)), weight = as.vector(P)) })) %>%
  filter(weight > 0) %>% arrange(marker, taxon, type)
hw(w, "template_species_type_weights")

# ---- 3. detection histories (libraries used in the fit; reads at species rank, summed over ESVs)
SP <- c("Hypsoblennius jenkinsi", "Cebidichthys violaceus")
rd <- bind_rows(lapply(MARKERS, function(m) arc(sprintf("%s_reads", m)) %>% filter(taxon %in% SP, rank == "species", reads > 0, library_id %in% L$library_id[L$marker == m]) %>%
  mutate(marker = m) %>% group_by(marker, library_id, taxon) %>% summarise(reads = sum(reads), n_esvs = n_distinct(esv_id), .groups = "drop")))
lib <- L %>% transmute(library_id, marker, location, site, Site, visit_id, bag = bagkey, filter = Barcode, filter_type = ftype, volume_ml, batch, platform,
                       pcr_count = k, library_reads_all_esvs = reads_all_esvs)
dh <- bind_rows(lapply(SP, function(s) lib %>% mutate(taxon = s) %>% left_join(rd %>% filter(taxon == s) %>% select(-taxon), by = c("marker", "library_id"), relationship = "one-to-one"))) %>%
  mutate(reads = coalesce(reads, 0L), n_esvs = coalesce(n_esvs, 0L), detected = reads > 0,
         in_presence_data = paste(marker, library_id, taxon) %in% c(paste("12S", B$pres$`12S`$library_id, B$pres$`12S`$taxon), paste("COI", B$pres$COI$library_id, B$pres$COI$taxon))) %>%
  arrange(taxon, marker, location, visit_id, bag, filter, library_id)
lg("detections:", paste(capture.output(print(dh %>% group_by(taxon, marker) %>% summarise(libs = n(), det = sum(detected), .groups = "drop"))), collapse = "\n"))
hw(dh, "detection_histories_Hjenkinsi_Cviolaceus")

# ---- 4. abalone and owl limpet
tax <- arc("taxa"); hit <- tax %>% filter(genus %in% c("Haliotis", "Lottia") | family %in% c("Haliotidae", "Lottiidae")) %>% pull(taxon)
calls <- bind_rows(lapply(c("12S", "18S", "COI"), function(m) arc(sprintf("%s_reads", m)) %>% filter(taxon %in% hit, reads > 0) %>% transmute(marker = m, esv_id, library_id, Barcode, Site, taxon, rank, reads))) %>%
  left_join(arc("libraries") %>% distinct(marker, library_id, batch), by = c("marker", "library_id"), relationship = "many-to-one") %>%
  group_by(marker, taxon, rank) %>% summarise(n_libraries = n_distinct(library_id), n_filters = n_distinct(Barcode), n_esvs = n_distinct(esv_id), reads = sum(reads),
    n_libraries_in_design_fit = n_distinct(library_id[paste(marker, library_id) %in% paste(L$marker, L$library_id)]), sites = paste(sort(unique(Site)), collapse = "; "),
    batches = paste(sort(unique(batch)), collapse = "/"), .groups = "drop") %>% arrange(marker, rank, taxon)
hw(calls, "abalone_limpet_calls")
TGT <- c("Haliotis cracherodii", "Haliotis rufescens", "Lottia gigantea")
sc <- arc("ESV_site_calls"); cand <- bind_rows(lapply(TGT, function(s) { x <- sc[grepl(s, sc$plausible_top, fixed = TRUE), ]
  if (!nrow(x)) return(data.frame(species = s, marker = NA, n_esv_sites = 0L))
  pr <- as.numeric(sub(paste0(".*", s, " ([0-9.]+).*"), "\\1", x$plausible_top))
  x %>% mutate(prob = pr) %>% group_by(marker) %>% summarise(species = s, n_esv_sites = n(), n_esvs = n_distinct(esv_id), max_posterior_as_candidate = max(prob),
    calls_made = paste(names(sort(table(paste0(consensus_taxon, " (", consensus_rank, ")")), decreasing = TRUE)), collapse = "; "),
    sites = paste(sort(unique(Site)), collapse = "; "), .groups = "drop") })) %>% select(species, everything())
hw(cand, "abalone_limpet_candidates")
# Section 4 reads the reference cache, not the data deposit. Env SKIP_REFCACHE=1 (r2, 2026-10-13): skip it; the file is carried forward
# unchanged and labelled "sourced: reference cache, not deposit" (it does not depend on the data release).
if (Sys.getenv("SKIP_REFCACHE") != "1") {
ref <- bind_rows(lapply(c("Haliotis", "Lottia"), function(g) { x <- readRDS(list.files(REFCACHE, paste0("^", g, "_COI_Leray_.*_meta\\.rds$"), full.names = TRUE))
  x %>% filter(species %in% TGT) %>% group_by(species) %>% summarise(n_accessions = n(), n_in_barcode_range = sum(in_barcode_range), length_min = min(slen), length_max = max(slen),
    accessions = paste(acc_version, collapse = ";"), .groups = "drop") }))
scr <- readRDS(paste0(REFCACHE, "reference_accession_cache.rds"))
ref <- ref %>% rowwise() %>% mutate(n_screened = sum(scr$accession %in% strsplit(accessions, ";")[[1]] | sub("\\.\\d+$", "", scr$accession) %in% sub("\\.\\d+$", "", strsplit(accessions, ";")[[1]])),
  screen_cols = paste(intersect(c("verdict", "flag", "decision", "status"), names(scr)), collapse = "/")) %>% ungroup() %>%
  mutate(source = "TaxaLikely::fetch_ncbi_reference_sequences per-genus COI (Leray region) cache, CalIntertidal/cache_reference")
hw(ref %>% select(species, n_accessions, n_in_barcode_range, length_min, length_max, n_screened, source, accessions), "abalone_limpet_COI_reference")
}

# ---- 5. effort to date (libraries used in the fit). Effective PCRs = run multiplier x m_eff(k, a), in JVB6097 units, as in the template likelihood
rmv <- setNames(rmt$run_multiplier, paste(rmt$marker, rmt$batch))
stopifnot(!anyDuplicated(paste(L$marker, L$library_id))); fl <- L %>% mutate(M = rmv[paste(marker, batch)] * m_eff(k, a)) %>%
  group_by(marker, location, site, Site, visit_id, bag = bagkey, filter = Barcode, filter_type = ftype, batch, platform) %>%
  summarise(libraries = n(), pcrs = sum(k), pcrs_per_library = paste(k, collapse = "/"), effective_pcrs_reference_units = sum(M), .groups = "drop") %>%
  arrange(location, visit_id, marker, bag, filter)
hw(fl, "effort_to_date_by_filter")
es <- fl %>% group_by(location, site, marker) %>% summarise(visits = n_distinct(visit_id), bags = n_distinct(bag), filters = n(), capsule_filters = sum(filter_type == "capsule"),
  libraries = sum(libraries), pcrs = sum(pcrs), effective_pcrs_reference_units = sum(effective_pcrs_reference_units),
  nanopore_filters = sum(platform == "nanopore"), illumina_filters = sum(platform == "illumina"), .groups = "drop") %>% arrange(site != "COAST", location, marker)
hw(es, "effort_to_date_by_site"); print(as.data.frame(es %>% filter(site != "COAST")))
print(as.data.frame(calls)); print(as.data.frame(cand %>% select(-sites))); if (exists("ref")) print(as.data.frame(ref %>% select(-accessions)))
