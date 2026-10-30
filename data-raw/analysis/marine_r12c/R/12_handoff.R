# 12_handoff.R -- data for Paper A's Pt. Conception worked case (request from the paper-writing chat, 2026-10-10).
# Writes tables to PtConception/Multi_site_multi_marker/handoff/ (never overwrites). Reads only the data deposit (r12c since 2026-10-25; earlier r2, the R1009f archive) and this study's outputs.
source("common.R"); source("model.R")
HO <- Sys.getenv("R12C_HANDOFF_OUT", paste0(ROOT, "handoff_out/")); dir.create(HO, showWarnings = FALSE, recursive = TRUE)  # PORTED: was Kevin's Mac path
# never overwrites: an existing file with identical content is left alone; a different one stops the script
hw <- function(x, name) { p <- paste0(HO, name, ".csv"); tf <- tempfile(fileext = ".csv"); utils::write.csv(x, tf, row.names = FALSE, na = "")
  if (file.exists(p)) { if (unname(tools::md5sum(p)) == unname(tools::md5sum(tf))) { lg("unchanged", name); return(invisible()) }; stop("exists and differs, not overwritten: ", p) }
  file.copy(tf, p); lg("wrote", name) }
fits <- readRDS(paste0(FITDIR, "fits.rds")); B <- readRDS(paste0(FITDIR, "build.rds")); C <- read_costs(); A0 <- 0.03
L <- B$L; S <- B$S; L0 <- arc("libraries") %>% filter(marker %in% MARKERS)

# ---- 1. design as sampled: event x site x marker (all libraries in the archive for the 2 markers, flagged if not used by the fit)
S2 <- S %>% select(Barcode, bagkey, ftype, event, Date)
d1 <- L0 %>% left_join(S2, by = "Barcode") %>% filter(!is.na(bagkey)) %>%
  mutate(used_in_fit = library_id %in% L$library_id, pooled = ifelse(n_pcr == 1, "independent", ifelse(n_pcr == 0, "pool re-sequencing singles (dropped)", "pooled")),
         n_pcr_evidence = sub("^(record|decision|inferred|unresolved).*$", "\\1", n_pcr_status))
design <- d1 %>% group_by(event, visit_id, Site, marker, batch) %>%
  summarise(date_first = min(Date.x), date_last = max(Date.x), bags = n_distinct(bagkey), filters = n_distinct(Barcode),
            filters_per_bag = paste(sort(unique(table(bagkey[!duplicated(Barcode)]))), collapse = "/"),
            filter_types = paste(sort(unique(ftype)), collapse = "/"),
            extracts_per_filter = paste(sort(unique(table(paste(Barcode, ifelse(is.na(extraction) | extraction == "", "A", extraction)))[0] + 0)), collapse = "/"),
            n_extractions = n_distinct(paste(Barcode, ifelse(is.na(extraction) | extraction == "", "A", extraction))),
            libraries = n(), libraries_used_in_fit = sum(used_in_fit), libraries_per_filter = paste(sort(unique(table(Barcode))), collapse = "/"),
            pcrs_per_library = paste(sort(unique(n_pcr)), collapse = "/"), pooling = paste(sort(unique(pooled)), collapse = "; "),
            n_pcr_evidence = paste(sort(unique(n_pcr_evidence)), collapse = "; "), n_pcr_status = paste(sort(unique(n_pcr_status)), collapse = " | "), .groups = "drop") %>%
  mutate(extracts_per_filter = round(n_extractions / filters, 2), platform = ifelse(batch %in% names(PLATFORM), "Illumina", "Nanopore"),
         site_type = ifelse(Site %in% names(PC_SITES), "Pt. Conception site", "coastline location")) %>%
  arrange(date_first, Site, marker)
hw(design, "design_as_sampled")

# ---- 2. which data inform each level (counts of replicated units in the data the fit used)
info <- bind_rows(lapply(MARKERS, function(m) { x <- L[L$marker == m, ]
  data.frame(marker = m,
    level = c("site-visit presence (phi)", "bag (q)", "filter (r)", "PCR (p)", "PCR (p)", "run multipliers"),
    informed_by = c("site-visits (Pt. Conception) and locations (coast)", "site-visits with >= 2 bags", "bags with >= 2 filters",
                    "filters with >= 2 independent single-PCR libraries", "libraries with k > 1 (pools; enter via 1-(1-p)^(c*m_eff(k,a)))", "batches other than JVB6097"),
    n_units = c(n_distinct(paste(x$site, x$visit_id, x$location)),
                sum(table(unique(x[, c("site", "visit_id", "bagkey")])[, c("site", "visit_id")]) >= 2),
                sum(tapply(x$Barcode, x$bagkey, n_distinct) >= 2),
                sum(tapply(x$library_id[x$k == 1], x$Barcode[x$k == 1], length) >= 2, na.rm = TRUE),
                sum(x$k > 1), n_distinct(x$batch) - 1),
    batches = c("all", paste(sort(unique(x$batch[x$bagkey %in% names(which(tapply(x$Barcode, x$bagkey, n_distinct) >= 1))])), collapse = "/"),
                paste(sort(unique(x$batch[x$bagkey %in% names(which(tapply(x$Barcode, x$bagkey, n_distinct) >= 2))])), collapse = "/"),
                paste(sort(unique(x$batch[x$k == 1])), collapse = "/"), paste(sort(unique(x$batch[x$k > 1])), collapse = "/"),
                paste(setdiff(sort(unique(x$batch)), REF_BATCH), collapse = "/"))) }))
info$batches[info$level == "bag (q)"] <- "Government Point visits (7) and coastline locations (21); Cojo and Jalama have 1 bag per visit"
hw(info, "fit_data_by_level")

# ---- 3. per-species parameters (empirical-Bayes posterior over the fitted mixture support), by marker, a = 0.03 (main) and the bracket
wq <- function(v, w, pr) { o <- order(v); v <- v[o]; cw <- cumsum(w[o]) / sum(w); sapply(pr, function(p) v[which(cw >= p - 1e-12)[1]]) }
sp <- bind_rows(lapply(MARKERS, function(m) bind_rows(lapply(A_GRID, function(a) { f <- fits[[paste(m, a)]]; ce <- f$cells; ce$r_capsule <- odds_mult(ce$r, f$kap_capsule)
  pars <- c(phi_GP = "phi_GP", phi_CJ = "phi_CJ", phi_JB = "phi_JB", phi_COAST = "phi_COAST", theta_bag = "q", theta_filter_disc = "r", theta_filter_capsule = "r_capsule", theta_pcr = "p")
  bind_rows(lapply(rownames(f$post), function(t) { w <- f$post[t, ]; row <- data.frame(marker = m, a = a, taxon = t)
    for (nm in names(pars)) { v <- ce[[pars[[nm]]]]; q <- wq(v, w, c(.025, .5, .975))
      row[[paste0(nm, "_mean")]] <- sum(w * v) / sum(w); row[[paste0(nm, "_q025")]] <- q[1]; row[[paste0(nm, "_q50")]] <- q[2]; row[[paste0(nm, "_q975")]] <- q[3] }
    row })) }))))
hw(sp, "species_parameters")
# the fitted mixture itself (support points and weights): the joint distribution across levels, correlations preserved
types <- bind_rows(lapply(MARKERS, function(m) bind_rows(lapply(A_GRID, function(a) { f <- fits[[paste(m, a)]]; ce <- f$cells
  data.frame(marker = m, a = a, type = seq_len(nrow(ce)), weight_observed = colSums(f$post), weight_total_nu = ce$nu, prob_never_detected = ce$P0,
             phi_GP = ce$phi_GP, phi_CJ = ce$phi_CJ, phi_JB = ce$phi_JB, phi_COAST = ce$phi_COAST, theta_bag = ce$q, theta_filter_disc = ce$r,
             theta_filter_capsule = odds_mult(ce$r, f$kap_capsule), theta_pcr = ce$p) }))))
hw(types, "species_types_mixture")
# community-level summary: observed, fitted total, unseen, Chao2 lower bound, log-likelihood
chao2 <- function(m) { x <- L[L$marker == m, ]; pr <- B$pres[[m]] %>% left_join(x %>% distinct(library_id, site, visit_id, location), by = "library_id")
  inc <- pr %>% distinct(taxon, site, visit_id, location) %>% count(taxon); Tn <- nrow(distinct(x, site, visit_id, location)); Q1 <- sum(inc$n == 1); Q2 <- sum(inc$n == 2)
  (Tn - 1) / Tn * Q1^2 / (2 * Q2) }
comm <- bind_rows(lapply(MARKERS, function(m) bind_rows(lapply(A_GRID, function(a) { f <- fits[[paste(m, a)]]
  data.frame(marker = m, a = a, species_observed = f$Nobs, species_total_fitted = f$Ntot, unseen_fitted = f$Ntot - f$Nobs, unseen_chao2_lower_bound = chao2(m),
             loglik = f$ll, kappa_capsule = f$kap_capsule, n_support_types = nrow(f$cells)) }))))
hw(comm, "community_summary")
rm_ <- utils::read.csv(paste0(FITDIR, "t_run_multipliers.csv")); hw(rm_, "run_multipliers")

# ---- 4. species list and target status
tax <- arc("taxa"); H <- arc("species_habitat")
det <- utils::read.csv(if (file.exists(paste0(FITDIR, "t_species_detected.csv"))) paste0(FITDIR, "t_species_detected.csv") else paste0(FITDIR, "superseded_warmstart/t_species_detected.csv"))   # written by 01_build (data only; unaffected by the warm-start fit bug)
spl <- det %>% left_join(tax %>% filter(rank == "species") %>% select(taxon, kingdom, phylum, class, order, family, genus), by = "taxon") %>%
  left_join(H %>% select(species, rocky_intertidal, rocky_intertidal_class, habitat_confidence), by = c("taxon" = "species")) %>%
  mutate(target = intertidal) %>% arrange(marker, desc(target), phylum, taxon)
hw(spl, "species_list_targets")

# ---- 5. costs: the vector used (technician hours, default since 2026-10-11) beside the two earlier vectors
cst <- utils::read.csv(COSTS_FILE, stringsAsFactors = FALSE)
old <- c(travel_location = 400, travel_site = 60, sample = 10, filter_disc = 40, filter_capsule = 55, library = 125, pcr = 5)
chat <- c(travel_location = 135, travel_site = 7.5, sample = 6, filter_disc = 50, filter_capsule = NA, library = 50, pcr = 5)
rec <- data.frame(item = cst$item, unit = cst$unit, hours = cst$hours, usd = cst$cost_usd, source = cst$source, note = cst$note,
  earlier_kevin_2026_10_10_usd = unname(old[cst$item]), earlier_chat_scoping_usd = unname(chat[cst$item]))
rec$earlier_chat_scoping_usd[rec$item == "filter_disc"] <- 50; rec$note[rec$item == "filter_disc"] <- paste(rec$note[rec$item == "filter_disc"], "| the chat's $50 bundles filter + extraction; Kevin's earlier $40 excluded extraction")
hw(rec, "costs_reconciliation")

# ---- 6. comparison tables (copied from this study's outputs) + site complementarity incl. pairs + order of site addition
for (nm in c("t_rule_vs_grid_single_target_v3", "t_rule_vs_grid_assemblage_v3", "t_rule_vs_grid_finite_extract", "t_enumerate_small_levels", "t_enumerate_small_levels_budget_aware",
             "t_branch_test_sites", "t_unseen_sensitivity_a0.03", "t_extract_cap33_designs", "t_pc_best_by_budget_disc", "t_coast_best_by_budget_disc", "t_site_knowledge",
             "t_enumerate_small_levels_budget_aware_trigger6", "t_enumerate_small_levels_budget_aware_trigger10"))
  hw(utils::read.csv(paste0(OUT, nm, ".csv"), check.names = FALSE), sub("^t_", "cmp_", nm))
for (nm in c("t_level_tests", "t_calibration", "t_fit_summary", "t_fit_summary_a0.0833"))   # fit diagnostics: cost-free, from the fit folder
  hw(utils::read.csv(paste0(FITDIR, nm, ".csv"), check.names = FALSE), sub("^t_", "cmp_", nm))
cx <- bind_rows(lapply(MARKERS, function(m) { f <- fits[[paste(m, A0)]]; ce <- f$cells; W <- colSums(f$post); s3 <- c("GP", "CJ", "JB")
  sets <- list("GP", "CJ", "JB", c("GP", "CJ"), c("GP", "JB"), c("CJ", "JB"))
  bind_rows(lapply(sets, function(st) { inS <- Reduce(`|`, lapply(st, function(s) ce[[paste0("phi_", s)]] > 0)); out <- setdiff(s3, st)
    outS <- if (length(out)) Reduce(`|`, lapply(out, function(s) ce[[paste0("phi_", s)]] > 0)) else rep(FALSE, nrow(ce))
    data.frame(marker = m, a = A0, sites = paste(st, collapse = "+"), expected_observed_species_reachable = sum(W[inS]), reachable_only_here = sum(W[inS & !outS])) })) }))
hw(cx, "site_complementarity")
pcb <- utils::read.csv(paste0(OUT, "t_pc_best_by_budget_disc.csv")) %>% filter(a == A0) %>% select(markers, budget, V, GP, CJ, JB, f, k12, kCOI, S, cost)
hw(pcb, "site_order_by_budget")
# where the remaining rule losses sit: the level(s) whose integer count differs between rule and grid
lv <- function(g, r) { g <- as.integer(strsplit(g, "/")[[1]]); r <- as.integer(strsplit(r, "/")[[1]]); if (any(is.na(r)) || length(r) < 4) return("rule infeasible")
  d <- c("visit", "sample", "filter", "pcr")[g != r]; if (!length(d)) "none" else paste(d, collapse = "+") }
for (nm in c("t_rule_vs_grid_single_target_v3", "t_enumerate_small_levels_budget_aware")) { x <- utils::read.csv(paste0(OUT, nm, ".csv"))
  dcol <- if ("rule_n1_design" %in% names(x)) "rule_n1_design" else "proc_design"; pcol <- if ("P_rule_n1" %in% names(x)) "P_rule_n1" else "P_proc"
  x$levels_differing <- mapply(lv, x$grid_design, x[[dcol]]); x$loss <- x$P_grid - x[[pcol]]
  hw(x %>% group_by(marker, levels_differing) %>% summarise(n = n(), weighted_loss_sum = sum(weight * loss), mean_loss = mean(loss), max_loss = max(loss), .groups = "drop") %>%
       arrange(marker, desc(weighted_loss_sum)), paste0("cmp_loss_by_level_", ifelse(dcol == "proc_design", "budget_aware_enumeration", "visit_enumeration"))) }

# ---- 7. as-run design by era (Pt. Conception sites only)
era <- design %>% filter(site_type == "Pt. Conception site") %>%
  mutate(era = case_when(batch == "JVB2844" ~ "1 Oct 2023 Illumina, Sterivex replicate test", batch == "JVB3105" ~ "2 Jan 2024 Illumina, filtration/replication experiment (Schulte)",
                         batch == "JVB3506" ~ "3 May 2024 Illumina, independent PCRs", batch %in% c("JVB3735", "JVB5058", "JVB5059", "JVB5060", "JVB5061") ~ "4 Aug 2024 - mid 2025 Nanopore",
                         TRUE ~ "5 Sep 2025 on Nanopore, new PCR mix")) %>%
  group_by(era, marker, Site) %>% summarise(visits = n_distinct(visit_id), batches = paste(sort(unique(batch)), collapse = "/"), bags_per_visit = paste(sort(unique(bags)), collapse = "/"),
                                            filters_per_bag = paste(sort(unique(filters_per_bag)), collapse = "/"), libraries_per_filter = paste(sort(unique(libraries_per_filter)), collapse = "/"),
                                            pcrs_per_library = paste(sort(unique(pcrs_per_library)), collapse = "/"), pooling = paste(sort(unique(pooling)), collapse = "; "),
                                            filter_types = paste(sort(unique(filter_types)), collapse = "/"), .groups = "drop") %>% arrange(era, Site, marker)
hw(era, "as_run_design_by_era")

# ---- 8. within-visit co-detection among the three sites (rough check of the independence assumption)
cod <- bind_rows(lapply(MARKERS, function(m) { x <- L[L$marker == m & L$site %in% c("GP", "CJ", "JB"), ]
  x$visit_id <- regmatches(x$visit_id, regexpr("Event[0-9]+", x$visit_id))     # visit_id names the site; the event number links the three sites' same-trip samples
  pr <- B$pres[[m]] %>% filter(library_id %in% x$library_id) %>%
    left_join(x %>% distinct(library_id, site, visit_id), by = "library_id") %>% distinct(taxon, site, visit_id)
  sv <- x %>% distinct(site, visit_id)
  bind_rows(lapply(list(c("GP", "CJ"), c("GP", "JB"), c("CJ", "JB")), function(pp) {
    vs <- intersect(sv$visit_id[sv$site == pp[1]], sv$visit_id[sv$site == pp[2]]); if (length(vs) < 3) return(NULL)
    D1 <- table(factor(pr$taxon[pr$site == pp[1] & pr$visit_id %in% vs]), pr$visit_id[pr$site == pp[1] & pr$visit_id %in% vs])
    tx <- sort(unique(pr$taxon[pr$visit_id %in% vs & pr$site %in% pp])); M <- function(s) { z <- matrix(0, length(tx), length(vs), dimnames = list(tx, vs))
      y <- pr[pr$site == s & pr$visit_id %in% vs, ]; z[cbind(match(y$taxon, tx), match(y$visit_id, vs))] <- 1; z }
    X <- M(pp[1]); Y <- M(pp[2]); keep <- rowSums(X) > 0 & rowSums(Y) > 0; X <- X[keep, , drop = FALSE]; Y <- Y[keep, , drop = FALSE]
    same <- sum(X * Y); nv <- length(vs)
    # expected co-detections if the two sites' detections were independent across visits (species-specific rates): sum_t (sum_v X)(sum_v Y)/nv
    expd <- sum(rowSums(X) * rowSums(Y) / nv)
    data.frame(marker = m, sites = paste(pp, collapse = "-"), shared_visits = nv, species_detected_at_both_sometime = nrow(X),
               same_visit_codetections = same, expected_if_independent = expd, ratio_observed_to_expected = same / expd) })) }))
hw(cod, "codetection_within_visit")
writeLines(capture.output(sessionInfo()), paste0(HO, "sessionInfo_12_handoff.txt"))
lg("done")
