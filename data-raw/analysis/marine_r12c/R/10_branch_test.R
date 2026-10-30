# 10_branch_test.R -- paperA S1.5 (branch theorem) applied to the three Pt. Conception sites as branches, against the grid (03).
# A branch's rate against species t = nats per dollar from one sample at that site (1 filter, k PCRs, a -> 0), on one visit:
#   r_ts = -log(1 - phi_ts * q_t * r_t * (1 - (1 - p_t)^k)) / (cost of one sample)
# S1.5(b) small-budget test: the branch with the largest sum_t w_t r_ts takes the first dollars.
# S1.5(c) large-budget test: a branch keeps everything only if it is not beaten on its slowest target; a site with phi = 0 for some
#   species has rate 0 against it, so "places spread". Our objective counts a species ONCE across sites, so a site's distinct value is
#   the species it can reach that the others cannot (or reach far less): reported as the weight of species whose presence there is
#   the highest of the three sites.
source("common.R")
fits <- readRDS(paste0(FITDIR, "fits.rds")); C <- read_costs(); a <- 0.03; k <- 15
cs <- C[["sample"]] + C[["filter_disc"]] + C[["library"]] + k * C[["pcr"]]
out <- bind_rows(lapply(MARKERS, function(m) { f <- fits[[paste(m, a)]]; ce <- f$cells; W <- colSums(f$post)
  inner <- ce$q * ce$r * (1 - (1 - ce$p)^k)
  bind_rows(lapply(c("GP", "CJ", "JB"), function(s) { ph <- ce[[paste0("phi_", s)]]; others <- setdiff(c("GP", "CJ", "JB"), s)
    mo <- pmax(ce[[paste0("phi_", others[1])]], ce[[paste0("phi_", others[2])]])
    data.frame(marker = m, site = s, small_budget_rate_species_per_1000usd = 1000 * sum(W * -log(1 - ph * inner)) / cs,
               species_reachable = sum(W[ph > 0]), species_only_here = sum(W[ph > 0 & mo == 0]), species_best_here = sum(W[ph > mo])) })) }))
wcsv(out, "t_branch_test_sites"); print(out)
