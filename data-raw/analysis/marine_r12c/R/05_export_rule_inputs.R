# 05_export_rule_inputs.R -- inputs for the analytic allocation rule (MEE draft), per marker and extract fraction a:
#  (1) the FITTED community distribution (incl. unseen species): one row per species type (cell) with its weight nu
#      (expected number of species of that type) and its conditional theta at every level;
#  (2) per observed species: posterior-mean thetas (descriptive only; the rule should use (1), per the peer 2026-10-10);
#  (3) the per-unit marginal costs from costs.csv.
source("common.R"); source("model.R")
fits <- readRDS(paste0(FITDIR, "fits.rds")); C <- read_costs()
lev <- c(theta_site_visit_GP = "phi_GP", theta_site_visit_CJ = "phi_CJ", theta_site_visit_JB = "phi_JB", theta_location_coast = "phi_COAST",
         theta_bag = "q", theta_filter_disc = "r", theta_pcr = "p")
D <- list(); O <- list()
for (k in names(fits)) { f <- fits[[k]]; ce <- f$cells; m <- sub(" .*", "", k)
  x <- data.frame(marker = m, a = f$a, weight_nu = ce$nu, setNames(ce[, lev], names(lev)))
  x$theta_filter_capsule <- odds_mult(ce$r, f$kap_capsule); x$theta_pool6 <- 1 - (1 - ce$p)^m_eff(6, f$a)
  D[[k]] <- x[x$weight_nu > 1e-4, ]
  W <- f$post; O[[k]] <- data.frame(marker = m, a = f$a, taxon = rownames(W), sapply(lev, function(v) as.vector(W %*% ce[[v]]))) }
wcsv(bind_rows(D), "rule_input_species_types"); wcsv(bind_rows(O), "rule_input_observed_species_posterior_means")
wcsv(data.frame(level = names(C), cost_usd = unname(C)), "rule_input_costs")
print(bind_rows(D) %>% group_by(marker, a) %>% summarise(types = n(), species_total = sum(weight_nu), .groups = "drop"))
