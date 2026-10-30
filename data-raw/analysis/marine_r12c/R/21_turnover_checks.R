# 21_turnover_checks.R -- follow-up to 20 (2026-10-11): a finer profile for the decay length Lc, and robustness of the Pt. Conception step.
#   profile: Lc on 150 log-spaced values, 100-20,000 km (20's 40-value grid was too coarse to bracket the 1.92 drop).
#   step:    leave one location out; Jalama dropped (the only north site within 100 km of a south one); San Miguel I. counted as north
#            (cold-water island); coastline locations only (no Pt. Conception sites); pairs under 700 km only.
# Coastal distance, Pt. Conception sites on Event3, as in 20. Composite likelihood: the profile interval is too narrow (see 20).
Sys.setenv(TURNOVER_DEFS_ONLY = "1"); source("20_turnover.R")
LF <- exp(seq(log(100), log(20000), length.out = 150)); PR <- list(); ST <- list()
stepfit <- function(R, label) { f <- fit1(R, "dist_coast"); g <- fit1(R, "dist_coast", step = FALSE)
  data.frame(label = label, rows = nrow(R), cross_rows = sum(R$cross), Lc_km = exp(f$th[3]), step_multiplier = exp(f$th[4]), dll_step = g$nll - f$nll) }
for (m in c(MARKERS, "12S+COI")) {
  mk <- if (m == "12S+COI") MARKERS else m
  R <- bind_rows(lapply(mk, function(z) build_pairs(z, PC_EVENT)$R)) %>% mutate(taxon = paste(marker, taxon))
  full <- fit1(R, "dist_coast"); pr <- sapply(LF, function(Lc) fit1(R, "dist_coast", Lfix = Lc, start = full$th)$nll) - full$nll
  ins <- LF[pr <= 1.92]; PR[[m]] <- data.frame(marker = m, Lc_km = LF, drop = pr)
  lg(m, "Lc", round(exp(full$th[3])), "profile", round(min(ins)), "-", round(max(ins)))
  ST[[length(ST) + 1]] <- stepfit(R, "all") %>% mutate(marker = m, Lc_profile_lo = min(ins), Lc_profile_hi = max(ins))
  for (z in sort(unique(R$A))) ST[[length(ST) + 1]] <- stepfit(R[R$A != z & R$B != z, ], paste("without", z)) %>% mutate(marker = m)
  ST[[length(ST) + 1]] <- stepfit(R[R$A != "Jalama Beach" & R$B != "Jalama Beach", ], "Jalama dropped") %>% mutate(marker = m)
  sm <- "Cuyler Harbor, San Miguel I."; side2 <- setNames(LOC$side, LOC$location); side2[sm] <- "north of Pt. Conception"
  R2 <- R; R2$cross <- as.numeric(side2[R2$A] != side2[R2$B]); ST[[length(ST) + 1]] <- stepfit(R2, "San Miguel I. north") %>% mutate(marker = m)
  pcs <- names(PC_SITES); ST[[length(ST) + 1]] <- stepfit(R[!R$A %in% pcs & !R$B %in% pcs, ], "coastline locations only") %>% mutate(marker = m)
  ST[[length(ST) + 1]] <- stepfit(R[R$dist_coast < 700, ], "pairs under 700 km") %>% mutate(marker = m) }
ST <- bind_rows(ST); wcsv(bind_rows(PR), "t_turnover_profile_fine"); wcsv(ST, "t_turnover_step_checks")
print(ST %>% filter(!grepl("^without", label)) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())
print(ST %>% filter(grepl("^without", label)) %>% group_by(marker) %>% summarise(step_min = min(step_multiplier), step_max = max(step_multiplier),
  dll_min = min(dll_step), which_min_dll = label[which.min(dll_step)], Lc_min = min(Lc_km), Lc_max = max(Lc_km)) %>% as.data.frame())
