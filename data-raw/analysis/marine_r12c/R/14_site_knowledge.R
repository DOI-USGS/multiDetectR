# 14_site_knowledge.R -- two versions of choosing sites within a location (Kevin Lafferty, 2026-10-11):
#   known habitat:   the user knows which sites are best and picks them, with any bag split (GP / CJ / JB as fitted branches) = 03's best design.
#   unknown habitat: the user picks only HOW MANY sites (s) and bags per site (b, equal: the sites can't be told apart); which s of the
#                    3 sites they get is random, so expected species = mean of S over the choose(3, s) site subsets (cost is the same for each).
# Markers are a study-wide choice (every filter gets the same markers), as in 03. From the saved a = 0.03 disc design points (03).
source("common.R")
a <- 0.03; BUD <- c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000)
P <- readRDS(paste0(OUT, sprintf("pc_points_a%s_disc.rds", a)))
P$s <- (P$GP > 0) + (P$CJ > 0) + (P$JB > 0); P$bmax <- pmax(P$GP, P$CJ, P$JB)
P$equal <- with(P, (GP == 0 | GP == bmax) & (CJ == 0 | CJ == bmax) & (JB == 0 | JB == bmax))
U <- P[P$equal, ] %>% group_by(markers, V, f, k12, kCOI, s, b = bmax) %>%
  summarise(n_subsets = n(), cost = first(cost), S_min = min(S), S = mean(S), .groups = "drop") %>%
  filter(n_subsets == choose(3, s))                        # every subset of s sites must be on the grid
out <- list()
for (mk in c("12S", "COI", "12S+COI")) for (Bd in BUD) {
  X <- P[P$markers == mk & P$cost <= Bd, ]; Y <- U[U$markers == mk & U$cost <= Bd, ]; if (!nrow(X) || !nrow(Y)) next
  k <- X[which.max(X$S), ]; u <- Y[which.max(Y$S), ]
  out[[length(out) + 1]] <- data.frame(markers = mk, budget = Bd,
    known_design = sprintf("V%d %d/%d/%d f%d k%s|%s", k$V, k$GP, k$CJ, k$JB, k$f, ifelse(is.na(k$k12), "", k$k12), ifelse(is.na(k$kCOI), "", k$kCOI)), known_S = k$S, known_cost = k$cost,
    unknown_design = sprintf("V%d %d site(s) x %d bags f%d k%s|%s", u$V, u$s, u$b, u$f, ifelse(is.na(u$k12), "", u$k12), ifelse(is.na(u$kCOI), "", u$kCOI)),
    unknown_S = u$S, unknown_S_worst_subset = u$S_min, unknown_cost = u$cost, value_of_habitat_knowledge = k$S - u$S) }
R <- bind_rows(out); wcsv(R, "t_site_knowledge"); print(R %>% mutate(across(where(is.numeric), ~ round(.x, 1))), right = FALSE)
