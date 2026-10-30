# 18_a_from_singles.R -- peer (paper chat) 2026-10-11: estimate the extract fraction a from single-reaction data only (one run, ~10
# separately sequenced reactions per tube), with the species-type mixture, and give its profile interval.
# (1) In the model as fitted (02), a enters only through m_eff(k, a) = (1 - (1 - a)^k) / a per LIBRARY; a single has k = 1, so m_eff = 1
#     for every a and singles carry no information on a. Checked by fitting at three values of a.
# (2) Shared extract: the reactions of one tube draw from one extract, so reaction number i (its order among the tube's singles, by
#     replicate id; assumed to be the pipetting order) counts as (1 - a)^(i - 1) reactions. The ten then sum to m_eff(10, a), as a pool
#     of ten would. a is profiled two-sided (a < 0 = later reactions detect MORE); 95% interval = 1.92 drop.
#     "normalised" (the estimate): each tube's weights are rescaled to mean 1, so a changes only the TREND across reactions. Without that
#     ("unnormalised"), a also shrinks every reaction's hazard, which the coarse p grid cannot absorb: the profile then snaps between
#     grid points by several units (first, aborted run 2026-10-11: log_18_aborted_unnormalised.txt). Both are written.
# (3) Independence check: given its type, the j of n singles detecting a species in a tube is Binomial(n, p). Observed vs expected
#     distribution of j (given j >= 1) from the a = 0 fit; excess at j = 1 and j = n would mean dependence between reactions.
# Data: May 2024 (JVB3506; 12S and COI) and, as a replicate, Jan 2024 12S (JVB3105). One run each, so no run multiplier.
source("common.R"); source("model.R")
OUT <- FITDIR; B <- readRDS(paste0(FITDIR, "build.rds")); IT <- 1500
SETS <- list(c("12S", "JVB3506"), c("COI", "JVB3506"), c("12S", "JVB3105"))
AG <- c(-0.06, -0.04, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2)
flat <- list(); prof <- list(); gof <- list(); info <- list()
for (st in SETS) { m <- st[1]; bt <- st[2]
  L <- B$L %>% filter(marker == m, batch == bt, k == 1) %>% group_by(Barcode) %>% mutate(order = rank(replicate, ties.method = "first"), n_singles = n()) %>% ungroup() %>% as.data.frame()
  pr <- B$pres[[m]] %>% filter(library_id %in% L$library_id)
  info[[paste(m, bt)]] <- data.frame(marker = m, batch = bt, tubes = n_distinct(L$Barcode), libraries = nrow(L), species = n_distinct(pr$taxon),
    singles_per_tube = paste(range(L$n_singles), collapse = "-"), visits = paste(unique(L$visit_id), collapse = "; "))
  lg(m, bt, ": tubes", n_distinct(L$Barcode), " libraries", nrow(L), " species", n_distinct(pr$taxon))
  fl <- sapply(c(0.001, 0.03, 0.1), function(a) fit_model(L, pr, a, iters = IT, keep_post = FALSE)$ll)       # (1) no M_lib: a should not matter
  flat[[paste(m, bt)]] <- data.frame(marker = m, batch = bt, a = c(0.001, 0.03, 0.1), ll = fl, ll_minus_first = fl - fl[1])
  P <- bind_rows(lapply(c("normalised", "unnormalised"), function(v) bind_rows(lapply(AG, function(a) {                                # (2)
    w <- (1 - a)^(L$order - 1); if (v == "normalised") w <- w / ave(w, L$Barcode)
    f <- fit_model(L, pr, 0.03, iters = IT, keep_post = FALSE, M_lib = w)
    lg(m, bt, v, "a", a, "ll", round(f$ll, 3)); data.frame(marker = m, batch = bt, version = v, a = a, ll = f$ll, U = f$unseen) }))))
  P <- P %>% group_by(version) %>% mutate(drop = max(ll) - ll) %>% ungroup(); prof[[paste(m, bt)]] <- P
  f0 <- fit_model(L, pr, 0.03, iters = IT, M_lib = rep(1, nrow(L)))                                                # (3) a = 0, with posteriors
  ce <- f0$cells; det <- L %>% inner_join(pr, by = "library_id", relationship = "many-to-many") %>% count(Barcode, taxon, name = "j") %>%
    left_join(L %>% distinct(Barcode, n_singles), by = "Barcode")
  ex <- bind_rows(lapply(seq_len(nrow(det)), function(i) { n <- det$n_singles[i]; w <- f0$post[det$taxon[i], ]
    pj <- sapply(1:n, function(x) dbinom(x, n, ce$p)) / (1 - (1 - ce$p)^n)                                   # cells x j, given j >= 1
    data.frame(n_singles = n, j = 1:n, expected = as.vector(crossprod(w, pj))) }))
  G <- ex %>% group_by(n_singles, j) %>% summarise(expected = sum(expected), .groups = "drop") %>%
    left_join(det %>% count(n_singles, j, name = "observed"), by = c("n_singles", "j")) %>% mutate(observed = coalesce(observed, 0L))
  gof[[paste(m, bt)]] <- G %>% mutate(marker = m, batch = bt) }
I <- bind_rows(info); Fl <- bind_rows(flat); P <- bind_rows(prof); G <- bind_rows(gof)
# profile summary: MLE on the grid (interpolated by a parabola through the best point and its neighbours) and the 1.92 interval
S <- P %>% group_by(marker, batch, version) %>% group_modify(function(d, k) { d <- d[order(d$a), ]; i <- which.max(d$ll)
  ahat <- if (i > 1 && i < nrow(d)) { x <- d$a[(i - 1):(i + 1)]; y <- d$ll[(i - 1):(i + 1)]; cf <- coef(lm(y ~ x + I(x^2))); -cf[2] / (2 * cf[3]) } else d$a[i]
  end <- function(ix) { z <- d[ix, ]; j <- which(z$drop > 1.92)[1]; if (is.na(j)) return(NA_real_); approx(z$drop[(j - 1):j], z$a[(j - 1):j], xout = 1.92)$y }
  data.frame(a_hat = unname(ahat), lower = end(i:1), upper = end(i:nrow(d)), ll_at_0 = d$ll[d$a == 0], ll_max = max(d$ll),
             drop_at_0.03 = d$drop[d$a == 0.03], grid_min = min(d$a), grid_max = max(d$a)) }) %>% ungroup()
G2 <- G %>% group_by(marker, batch) %>% summarise(pairs = sum(observed), obs_j1 = sum(observed[j == 1]), exp_j1 = sum(expected[j == 1]),
  obs_jn = sum(observed[j == n_singles]), exp_jn = sum(expected[j == n_singles]),
  obs_mean_j = sum(observed * j) / sum(observed), exp_mean_j = sum(expected * j) / sum(expected), .groups = "drop")
wcsv(I, "t_a_from_singles_data"); wcsv(Fl, "t_a_from_singles_flat_check"); wcsv(P, "t_a_from_singles_profile"); wcsv(S, "t_a_from_singles_summary")
wcsv(G, "t_a_from_singles_j_distribution"); wcsv(G2, "t_a_from_singles_j_summary")
print(I); print(Fl); print(S %>% mutate(across(where(is.numeric), ~ signif(.x, 4))) %>% as.data.frame()); print(as.data.frame(G2))
print(as.data.frame(G %>% filter(n_singles == max(n_singles)) %>% mutate(expected = round(expected, 1))))
