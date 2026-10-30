# 28_template_sensitivity.R -- where does the template heterogeneity (23_template_fit.R) come from? alpha profiles (capsule and run multipliers held
# at the saved a = 0.03 fit, 1,500 iterations, as 23's first step) for:
#   12S without the two single-reaction runs (JVB3105, JVB3506): pooled libraries only, one library per filter;
#   COI with a within-run read floor on the May 2024 run (JVB3506, where cross-talk could not be assessed; CalIntertidal workflow session,
#   2026-10-11): a detection there needs >= FLOOR reads of the species in that library (20 and 50; the run's libraries hold ~10,000 reads).
# 23 already gives COI without JVB3506 (no heterogeneity: alpha -> Inf).
source("common.R"); source("model.R")
OUT <- FITDIR; B <- readRDS(paste0(FITDIR, "build.rds")); fits <- readRDS(paste0(FITDIR, "fits.rds")); a <- 0.03; IT <- 1500
KAP <- c(Inf, 8, 4, 2, 1, 0.5, 0.25)
prof <- function(m, L, pr, label) { f0 <- fits[[paste(m, a)]]
  bind_rows(lapply(KAP, function(k) { g <- fit_model(L, pr, a, run_mult = f0$run_mult, kap_capsule = f0$kap_capsule, iters = IT, keep_post = FALSE, tmpl = tmpl_gamma(k))
    lg(m, label, "alpha", k, "ll", round(g$ll, 2)); data.frame(marker = m, data = label, alpha = k, ll = g$ll) })) %>% mutate(dll_vs_no_heterogeneity = ll - ll[is.infinite(alpha)]) }
L <- B$L[B$L$marker == "12S", ]; L <- L[!L$batch %in% c("JVB3105", "JVB3506"), ]
R1 <- prof("12S", L, B$pres$`12S` %>% filter(library_id %in% L$library_id), "12S without JVB3105 and JVB3506 (pooled libraries only)")
L <- B$L[B$L$marker == "COI", ]; rd <- arc("COI_reads") %>% filter(library_id %in% L$library_id[L$batch == "JVB3506"], rank == "species", reads > 0) %>%
  group_by(library_id, taxon) %>% summarise(reads = sum(reads), .groups = "drop")
R2 <- bind_rows(lapply(c(20, 50), function(fl) { drop <- rd %>% filter(reads < fl) %>% select(library_id, taxon)
  pr <- B$pres$COI %>% anti_join(drop, by = c("library_id", "taxon")); lg("COI floor", fl, ": detections dropped", nrow(B$pres$COI) - nrow(pr))
  prof("COI", L, pr, sprintf("COI, JVB3506 detections need >= %d reads", fl)) %>% mutate(detections_dropped = nrow(B$pres$COI) - nrow(pr)) }))
R <- bind_rows(R1, R2); wcsv(R, "t_template_sensitivity"); print(R %>% mutate(across(where(is.double), ~ round(.x, 1))) %>% as.data.frame())
