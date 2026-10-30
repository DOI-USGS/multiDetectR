# 32_rule_visits_handoff.R -- the handoff's fixed-visit table and budget bands (12S+COI, template fit, hours), from script 29's outputs.
# Written 2026-10-25 for r12c. With V6RULE=r2 it reproduces the r2 handoff tables (built by hand on 2026-10-13/14) from out_r2 exactly.
# Rows, per V (Kevin 2026-10-11/14):
#   V = 1-4  at the existing budgets $1,500 / $3,000 / $4,000 / $5,000;
#   V = 6    at the lowest searched budget where the rule chooses 6 visits (Kevin 2026-10-25: "yes move V6"). Until then (r2, and the
#            first r12c table) it was $8,250, the R1009f budget, while that still chose 6: env V6RULE=r2 gives that version;
#   V = 12   at the lowest searched budget where the rule chooses 12 visits ("Move the 12 visit to be consistent with calculations");
#   then, for each V, the lowest and highest searched budget choosing it.
# rule_chooses_this_V = the rule's best V at that budget is this V; other_V_within_0.5_species = other V within 0.5 species there.
# Env CHECK_DIR: write there instead of OUT (used to check the reproduction of r2 without writing into out_r2).
source("common.R")
OUT <- paste0(FITDIR, "costs_hours_2026-10-11_template/")
if (nzchar(Sys.getenv("CHECK_DIR"))) OUT <- Sys.getenv("CHECK_DIR")
IN  <- paste0(FITDIR, "costs_hours_2026-10-11_template/")
al  <- utils::read.csv(paste0(IN, "t_rule_visit_bands_all_template.csv"), stringsAsFactors = FALSE)
bd  <- utils::read.csv(paste0(IN, "t_rule_visit_bands_template.csv"), stringsAsFactors = FALSE)
row_at <- function(v, b, def) { x <- al[al$budget == b, ]; r <- x[x$V == v, ]; stopifnot(nrow(r) == 1)
  near <- x$V[x$V != v & x$S >= max(x$S) - 0.5]
  data.frame(visits = v, definition = def, budget = b, GP = r$GP, CJ = r$CJ, JB = r$JB, filters_per_bag = r$f, pcr_12S = r$k12, pcr_COI = r$kCOI,
             cost_total = r$cost, cost_per_visit = r$cost / v, S_expected_observed = round(r$S, 2), rule_chooses_this_V = r$chosen,
             other_V_within_0.5_species = paste(sort(near), collapse = ","), stringsAsFactors = FALSE) }
ch <- al[al$chosen, ]
step250 <- function(v) sort(ch$budget[ch$V == v & ch$budget %% 250 == 0])
s6 <- step250(6)
if (Sys.getenv("V6RULE") == "r2") { b6 <- if (8250 %in% s6) 8250 else min(s6)
  def6 <- sprintf("%s ($250-step budgets choosing V on %s: %s)", if (b6 == 8250) "R1009f budget" else "lowest $250-step budget choosing V", RUN_LABEL, paste(s6, collapse = ", "))
} else { b6 <- min(ch$budget[ch$V == 6]); def6 <- "lowest budget choosing V (Kevin 2026-10-25)" }
fx <- rbind(row_at(1, 1500, "existing budget"), row_at(2, 3000, "existing budget"), row_at(3, 4000, "existing budget"), row_at(4, 5000, "existing budget"),
            row_at(6, b6, def6), row_at(12, min(ch$budget[ch$V == 12]), "lowest budget choosing V (Kevin 2026-10-14)"))
for (v in c(1, 2, 3, 4, 6, 12)) { b <- ch$budget[ch$V == v]
  fx <- rbind(fx, row_at(v, min(b), "lowest searched budget choosing V"), row_at(v, max(b), "highest searched budget choosing V")) }
tag <- paste0("12S_COI_template_", format(Sys.Date()), "_", RUN_LABEL)
w <- function(x, name) { p <- paste0(OUT, name, ".csv"); if (file.exists(p)) stop("exists, not overwritten: ", p); utils::write.csv(x, p, row.names = FALSE, na = ""); lg("wrote", p) }
if (Sys.getenv("V6RULE") != "r2") tag <- paste0(tag, "_V6low")   # the bands do not depend on V6RULE; written once, identical content left alone
w(fx, paste0("rule_visits_fixed_", tag))
bn <- paste0(OUT, "rule_visit_bands_12S_COI_template_", format(Sys.Date()), "_", RUN_LABEL, ".csv"); bx <- bd[, c("V", "budget_from", "budget_to")]
if (file.exists(bn)) stopifnot(identical(utils::read.csv(bn), bx)) else w(bx, sub("\\.csv$", "", basename(bn)))
print(fx)
