# 01_build.R -- effort and intertidal-species presence per library, Pt. Conception (3 sites) + the coastline locations.
source("common.R")
S <- arc("samples") %>%
  mutate(site = ifelse(Site %in% names(PC_SITES), unname(PC_SITES[Site]), "COAST"),
         location = ifelse(Site %in% names(PC_SITES), "Pt. Conception", Site),
         bagkey = paste(visit_id, Site, bag),                                  # archive key for a water sample
         ftype = ifelse(filter_type %in% c("Sterivex", "capsule"), "capsule", "disc"))   # R1009f says Sterivex, the r2 deposit says capsule
stopifnot(!anyNA(S$bag), all(nzchar(S$bag)))
L0 <- arc("libraries") %>% filter(marker %in% MARKERS)
L <- L0 %>% filter(n_pcr > 0,                                                  # drop May 2024 12S pools (singles used) and failed libraries whose rerun is kept
                   is.na(extraction) | extraction %in% c("", "A")) %>%         # Jan 2024 triple extractions: keep A (no extraction level in the hierarchy)
  left_join(S %>% select(Barcode, site, location, bagkey, ftype, volume_ml), by = "Barcode") %>%
  mutate(k = n_pcr, platform = ifelse(batch %in% names(PLATFORM), PLATFORM[batch], "nanopore"))
stopifnot(!anyNA(L$bagkey))
H <- arc("species_habitat"); keep_sp <- H$species[H$rocky_intertidal_class %in% HABITAT_KEEP]
pres <- list(); sp_tab <- list()
for (m in MARKERS) {
  r <- arc(sprintf("%s_reads", m)) %>% filter(rank == "species", reads > 0, library_id %in% L$library_id[L$marker == m])
  r$intertidal <- r$taxon %in% keep_sp
  sp_tab[[m]] <- r %>% group_by(taxon) %>% summarise(marker = m, intertidal = first(intertidal),
                   habitat_class = H$rocky_intertidal_class[match(first(taxon), H$species)], n_libraries = n_distinct(library_id),
                   n_locations = n_distinct(L$location[match(library_id, L$library_id)]), .groups = "drop")
  pres[[m]] <- r %>% filter(intertidal) %>% distinct(library_id, taxon)
  lg(m, ": libraries", sum(L$marker == m), " species", n_distinct(r$taxon), " intertidal+partly", n_distinct(pres[[m]]$taxon))
}
eff <- L %>% group_by(marker, location, site, visit_id, batch, platform) %>%
  summarise(bags = n_distinct(bagkey), filters = n_distinct(Barcode), capsule_filters = n_distinct(Barcode[ftype == "capsule"]),
            libraries = n(), k = paste(sort(unique(k)), collapse = "/"), .groups = "drop")
wcsv(eff, "t_effort_by_visit"); wcsv(bind_rows(sp_tab), "t_species_detected")
saveRDS(list(L = L, pres = pres, S = S), paste0(OUT, "build.rds"))
print(as.data.frame(eff %>% group_by(marker, site) %>% summarise(units = n(), bags = sum(bags), filters = sum(filters), libraries = sum(libraries), .groups = "drop")))
print(as.data.frame(eff %>% filter(marker == "COI") %>% count(site, bags_per_unit = bags)))
