# 31_detection_histories_inverts.R -- paper chat 2026-10-12 (revised request): detection histories for four iconic rocky-intertidal COI species,
# same format as detection_histories_Hjenkinsi_Cviolaceus.csv (30_handoff_species.R, section 3). One row per species x library used in the fit
# (both markers, every site); reads at species rank summed over ESVs (0 = not detected). Type weights for these species are already in
# template_species_type_weights.csv (30). Writes to the handoff folder; never overwrites.
source("common.R")
HO <- Sys.getenv("R12C_HANDOFF_OUT", paste0(ROOT, "handoff_out/"))  # PORTED: was Kevin's Mac path
hw <- function(x, name) { p <- paste0(HO, name, ".csv"); tf <- tempfile(fileext = ".csv"); utils::write.csv(x, tf, row.names = FALSE, na = "")
  if (file.exists(p)) { if (unname(tools::md5sum(p)) == unname(tools::md5sum(tf))) { lg("unchanged", name); return(invisible()) }; stop("exists and differs, not overwritten: ", p) }
  file.copy(tf, p); lg("wrote", name, nrow(x), "rows") }
B <- readRDS(paste0(FITDIR, "build.rds")); L <- B$L; stopifnot(!anyDuplicated(paste(L$marker, L$library_id)))
SP <- c("Endocladia muricata", "Pisaster ochraceus", "Pollicipes polymerus", "Strongylocentrotus purpuratus")
rd <- bind_rows(lapply(MARKERS, function(m) arc(sprintf("%s_reads", m)) %>% filter(taxon %in% SP, rank == "species", reads > 0, library_id %in% L$library_id[L$marker == m]) %>%
  mutate(marker = m) %>% group_by(marker, library_id, taxon) %>% summarise(reads = sum(reads), n_esvs = n_distinct(esv_id), .groups = "drop")))
lib <- L %>% transmute(library_id, marker, location, site, Site, visit_id, bag = bagkey, filter = Barcode, filter_type = ftype, volume_ml, batch, platform,
                       pcr_count = k, library_reads_all_esvs = reads_all_esvs)
dh <- bind_rows(lapply(SP, function(s) lib %>% mutate(taxon = s) %>% left_join(rd %>% filter(taxon == s) %>% select(-taxon), by = c("marker", "library_id"), relationship = "one-to-one"))) %>%
  mutate(reads = coalesce(reads, 0L), n_esvs = coalesce(n_esvs, 0L), detected = reads > 0,
         in_presence_data = paste(marker, library_id, taxon) %in% c(paste("12S", B$pres$`12S`$library_id, B$pres$`12S`$taxon), paste("COI", B$pres$COI$library_id, B$pres$COI$taxon))) %>%
  arrange(taxon, marker, location, visit_id, bag, filter, library_id)
stopifnot(all(dh$in_presence_data == dh$detected))                       # all four are intertidal targets: every detection is in the fit's data
w <- utils::read.csv(paste0(HO, "template_species_type_weights.csv")); stopifnot(all(SP %in% w$taxon[w$marker == "COI"]))
hw(dh, "detection_histories_COI_inverts_algae")
print(as.data.frame(dh %>% filter(detected) %>% group_by(taxon, marker, site) %>% summarise(site_visits = n_distinct(visit_id), filters = n_distinct(filter), libraries = n(), reads = sum(reads), .groups = "drop")))
