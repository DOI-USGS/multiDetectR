# Multiscale occupancy mixture (frame 2). Per species ("cell"):
#   phi_GP, phi_CJ, phi_JB : presence at that Pt. Conception site on a visit (branches; visit turnover included)
#   phi_COAST              : presence at a random coastline location on its (single) visit
#   q : a water sample (bag) carries it, given the site-visit does
#   r : a filter carries it, given its bag does (capsule: odds multiplier)
#   p : one PCR of the REFERENCE run detects it, given the filter carries it (a -> 0 limit; see m_eff in common.R)
# A library of k pooled PCRs from run b detects with 1 - (1-p)^(c_b * m_eff(k, a)).
# Estimator as in ../model.R (zero-truncated Poisson-mixture NPMLE by EM on a grid; unseen species estimated).
PHI <- c(0, .1, .3, .6, 1); QG <- c(.05, .15, .35, .65, .95); RG <- c(.05, .15, .35, .65, .95); PG <- c(.005, .02, .05, .1, .2, .4, .7)
odds_mult <- function(r, k) r * k / (1 - r + r * k)
# L: libraries of one marker (Barcode, site, visit_id, bagkey, ftype, library_id, batch, k); pres: library_id, taxon
fit_model <- function(L, pres, a, run_mult = NULL, kap_capsule = 1, iters = 1500, nu0 = NULL, qfix = NULL, rfix = NULL, keep_post = TRUE,
                      lambda = 0, add_unseen = 0, M_lib = NULL, tmpl = NULL, trace_at = NULL) {
  # lambda > 0: profile likelihood for the unseen total U = sum(nu * P0): maximise ll - lambda * U (fixed point nu = nu * cp / (1 - (1 - lambda) * P0)).
  # add_unseen: after EM, add that many unseen species to the cell with the largest P0 (the cheapest place for extra U; upper end of the profile).
  # M_lib: per-library effective PCR count replacing m_eff(L$k, a) (run multiplier still applied); used to let separate reactions share
  #   one extract (18_a_from_singles.R).
  # tmpl: template amount varying among filters (23_template_fit.R): list(G = class multipliers of the PCR hazard, w = class weights, mean 1),
  #   e.g. tmpl_discrete_gamma(kappa). Given the filter carries the species, its class multiplies every reaction's hazard -log(1-p) by G, so the
  #   reactions of one filter are dependent. NULL (default) = one class with G = 1, exactly the earlier model.
  #   tmpl = tmpl_gamma(alpha): G continuous, Gamma(alpha, rate alpha). P(no reaction detects | carries) = (1 + beta M)^(-alpha), beta = -log(1-p)/alpha,
  #   exact; a filter with several libraries (separately sequenced reactions) needs E[prod_det (1 - e^(-h G)) e^(-H0 G)], done by numerical integration
  #   over log G (cached by pattern).
  # trace_at: EM iterations at which to record ll and U (out$trace). Defaults leave the fit exactly as before.
  G <- expand.grid(q = if (is.null(qfix)) QG else qfix, r = if (is.null(rfix)) RG else rfix, p = PG); nq <- nrow(G)
  cb <- if (is.null(run_mult)) rep(1, nrow(L)) else ifelse(L$batch %in% names(run_mult), run_mult[L$batch], 1)
  L$M <- cb * (if (is.null(M_lib)) m_eff(L$k, a) else M_lib)
  tubes <- L %>% distinct(Barcode, site, visit_id, bagkey, ftype) %>% arrange(visit_id, bagkey, Barcode); tubes$ti <- seq_len(nrow(tubes))
  stopifnot(!anyDuplicated(tubes$Barcode))
  bags <- tubes %>% distinct(bagkey, visit_id, site); bags$bi <- seq_len(nrow(bags)); tubes$bi <- bags$bi[match(tubes$bagkey, bags$bagkey)]
  units <- bags %>% distinct(visit_id, site); units$ui <- seq_len(nrow(units)); bags$ui <- units$ui[match(paste(bags$visit_id, bags$site), paste(units$visit_id, units$site))]
  L$ti <- tubes$ti[match(L$Barcode, tubes$Barcode)]; Msum <- tapply(L$M, L$ti, sum)[as.character(tubes$ti)]
  lq1p <- log(1 - G$p); kt <- ifelse(tubes$ftype == "capsule", kap_capsule, 1)
  TG <- if (is.null(tmpl)) list(G = 1, w = 1) else tmpl; GAM <- identical(tmpl$type, "gamma")
  pnd <- if (GAM) function(M) (1 - M * lq1p / tmpl$alpha)^(-tmpl$alpha) else function(M) { s <- 0; for (g in seq_along(TG$G)) s <- s + TG$w[g] * exp(M * TG$G[g] * lq1p); s }   # P(no reaction detects | carries)
  up <- sort(unique(G$p)); ip <- match(G$p, up); gcache <- new.env()
  LF0 <- t(sapply(seq_len(nrow(tubes)), function(j) { rj <- odds_mult(G$r, kt[j]); log((1 - rj) + rj * pnd(Msum[j])) }))   # tubes x nq
  SB0 <- rowsum(LF0, tubes$bi, reorder = TRUE)
  LB0 <- log(matrix(1 - rep(G$q, each = nrow(SB0)), nrow(SB0)) + matrix(rep(G$q, each = nrow(SB0)), nrow(SB0)) * exp(SB0))
  LU0 <- rowsum(LB0, bags$ui, reorder = TRUE)
  LZ <- lapply(PHI, function(ph) log((1 - ph) + ph * exp(LU0)))
  Zs <- lapply(SITES, function(s) sapply(seq_along(PHI), function(k) colSums(LZ[[k]][units$ui[units$site == s], , drop = FALSE])))
  pr <- pres %>% filter(library_id %in% L$library_id); txs <- sort(unique(pr$taxon)); A <- vector("list", length(txs))
  for (i in seq_along(txs)) {
    det <- L$library_id %in% pr$library_id[pr$taxon == txs[i]]
    Ai <- lapply(seq_along(SITES), function(s) t(Zs[[s]])); dt <- unique(L$ti[det]); lt <- list()
    for (t in dt) { v <- log(odds_mult(G$r, kt[t]))
      if (is.null(tmpl)) { for (j in which(L$ti == t)) v <- v + if (det[j]) log(pmax(1 - exp(L$M[j] * lq1p), 1e-300)) else L$M[j] * lq1p }
      else if (GAM) { li <- which(L$ti == t)
        if (length(li) == 1) v <- v + log(pmax(1 - pnd(L$M[li]), 1e-300)) else {
          Md <- sort(L$M[li][det[li]]); H0 <- sum(L$M[li][!det[li]]); key <- paste(signif(H0, 12), paste(signif(Md, 12), collapse = ","))
          if (is.null(gcache[[key]])) assign(key, sapply(-log(1 - up), function(hp) gamma_pattern(hp * Md, hp * H0, tmpl$alpha)), envir = gcache)
          v <- v + log(pmax(gcache[[key]][ip], 1e-300)) } }
      else { vg <- matrix(log(TG$w), length(TG$G), nq)                                                     # classes x nq, log-sum-exp over classes
        for (j in which(L$ti == t)) { h <- outer(L$M[j] * TG$G, lq1p); vg <- vg + if (det[j]) log(pmax(1 - exp(h), 1e-300)) else h }
        mxg <- apply(vg, 2, max); v <- v + mxg + log(colSums(exp(vg - rep(mxg, each = nrow(vg))))) }
      lt[[as.character(t)]] <- v }
    db <- unique(tubes$bi[dt]); lb <- list()
    for (b in db) { v <- log(G$q) + SB0[b, ]; for (t in dt[tubes$bi[dt] == b]) v <- v - LF0[t, ] + lt[[as.character(t)]]; lb[[as.character(b)]] <- v }
    for (u in unique(bags$ui[db])) { s <- match(units$site[u], SITES); lu1 <- LU0[u, ]
      for (b in db[bags$ui[db] == u]) lu1 <- lu1 - LB0[b, ] + lb[[as.character(b)]]
      for (k in seq_along(PHI)) Ai[[s]][k, ] <- Ai[[s]][k, ] - LZ[[k]][u, ] + (if (PHI[k] == 0) -Inf else log(PHI[k]) + lu1) }
    A[[i]] <- Ai }
  ph <- do.call(expand.grid, rep(list(seq_along(PHI)), length(SITES))); ph <- ph[rowSums(ph > 1) > 0, , drop = FALSE]; names(ph) <- paste0("i", seq_along(SITES))
  cells <- cbind(ph[rep(seq_len(nrow(ph)), each = nq), , drop = FALSE], g = rep(seq_len(nq), times = nrow(ph)))
  LL <- matrix(0, length(txs), nrow(cells)); lP0 <- 0
  for (s in seq_along(SITES)) { ix <- cbind(cells[[paste0("i", s)]], cells$g)
    for (i in seq_along(txs)) LL[i, ] <- LL[i, ] + A[[i]][[s]][ix]
    lP0 <- lP0 + Zs[[s]][cbind(cells$g, cells[[paste0("i", s)]])] }
  P0 <- exp(lP0); mx <- apply(LL, 1, max); LL <- exp(LL - mx)                          # LL now holds the scaled likelihoods
  nu <- if (is.null(nu0)) rep(length(txs) / ncol(LL), ncol(LL)) else nu0
  trace <- NULL
  for (it in seq_len(iters)) { den <- pmax(as.vector(LL %*% nu), 1e-300); cp <- as.vector(crossprod(LL, 1 / den))
    if (it %in% trace_at) trace <- rbind(trace, data.frame(iter = it - 1, ll = sum(log(den) + mx) - sum(nu * (1 - P0)), U = sum(nu * P0), Ntot = sum(nu)))
    nu <- if (lambda == 0) nu * cp + nu * P0 else nu * cp / (1 - (1 - lambda) * P0) }
  if (add_unseen > 0) { j <- which.max(P0); nu[j] <- nu[j] + add_unseen / P0[j] }
  den <- pmax(as.vector(LL %*% nu), 1e-300); ll <- sum(log(den) + mx) - sum(nu * (1 - P0))
  for (s in seq_along(SITES)) cells[[paste0("phi_", SITES[s])]] <- PHI[cells[[paste0("i", s)]]]
  cells <- cbind(cells, G[cells$g, ]); cells$nu <- nu; cells$P0 <- P0
  out <- list(a = a, run_mult = run_mult, kap_capsule = kap_capsule, tmpl = tmpl, trace = trace, ll = ll, nu_full = nu, Nobs = length(txs), Ntot = sum(nu), unseen = sum(nu * P0), P0max = max(P0),
              tubes = tubes, bags = bags, units = units, taxa = txs, L = L, pres = pr)
  if (keep_post) { keep <- which(nu > 1e-7 * sum(nu)); post <- (LL[, keep, drop = FALSE] * rep(nu[keep], each = nrow(LL))) / den; rownames(post) <- txs
    out$cells <- cells[keep, ]; attr(out$cells, "tmpl") <- tmpl; out$post <- post }
  out }

# ---------- design prediction (Nanopore, reference-run units) ----------
# P(one filter detects | its bag carries); k PCRs pooled from one extract, fraction a each
# tmpl: the fit's template classes (attached to its cells by fit_model); NULL = one class, the earlier model
filt_det <- function(ce, a, k, ftype = "disc", kap_capsule = 1, tmpl = attr(ce, "tmpl")) {
  nd <- if (is.null(tmpl)) (1 - ce$p)^m_eff(k, a) else if (identical(tmpl$type, "gamma")) (1 - m_eff(k, a) * log(1 - ce$p) / tmpl$alpha)^(-tmpl$alpha) else {
    s <- 0; for (g in seq_along(tmpl$G)) s <- s + tmpl$w[g] * (1 - ce$p)^(m_eff(k, a) * tmpl$G[g]); s }
  odds_mult(ce$r, if (ftype == "capsule") kap_capsule else 1) * (1 - nd) }
# discrete gamma (Yang 1994, mean of each of K equal-probability classes) for template amount: shape kappa, mean 1; kappa = Inf -> NULL
tmpl_discrete_gamma <- function(kappa, K = 8) { if (!is.finite(kappa)) return(NULL); cut <- qgamma((0:K) / K, kappa, kappa)
  G <- K * diff(pgamma(cut * kappa, kappa + 1)); list(G = G / sum(G / K), w = rep(1 / K, K), kappa = kappa) }
tmpl_gamma <- function(alpha) if (!is.finite(alpha)) NULL else list(type = "gamma", alpha = alpha)
# E over G ~ Gamma(alpha, rate alpha) of prod_j (1 - exp(-hd_j G)) * exp(-h0 G): one filter's libraries, those in hd detecting; integrated over t = log G
gamma_pattern <- function(hd, h0, alpha) {
  f <- function(t) { g <- exp(t); lv <- dgamma(g, alpha, alpha, log = TRUE) + t - h0 * g; for (h in hd) lv <- lv + log(-expm1(-h * g)); exp(lv) }
  lo <- log(qgamma(1e-14, alpha, alpha)) - 5; hi <- log(qgamma(1 - 1e-14, alpha, alpha)) + 5
  z <- tryCatch(integrate(f, lo, hi, rel.tol = 1e-10, subdivisions = 2000L)$value, error = function(e) NA_real_)
  if (is.na(z)) { tt <- seq(lo, hi, length.out = 20001); z <- sum(f(tt)) * (tt[2] - tt[1]) }
  z }
tmpl_two_class <- function(w_low, rho) { hi <- 1 / (w_low * rho + 1 - w_low); list(G = c(rho * hi, hi), w = c(w_low, 1 - w_low), w_low = w_low, rho = rho) }
# Pt. Conception: V visits; b = named bags per site (GP, CJ, JB); f filters per bag; k PCRs per filter
det_pc <- function(ce, a, V, b, f, k, ftype = "disc", kap_capsule = 1) {
  pb <- ce$q * (1 - (1 - filt_det(ce, a, k, ftype, kap_capsule))^f); lm <- 0
  for (s in names(b)) if (b[[s]] > 0) lm <- lm + V * log(pmax(1 - ce[[paste0("phi_", s)]] * (1 - (1 - pb)^b[[s]]), 1e-300))
  1 - exp(lm) }
# Coastline: N exchangeable locations, one visit and one site each; b bags; f filters; k PCRs
det_coast <- function(ce, a, N, b, f, k, ftype = "disc", kap_capsule = 1) {
  pb <- ce$q * (1 - (1 - filt_det(ce, a, k, ftype, kap_capsule))^f)
  1 - (1 - ce$phi_COAST * (1 - (1 - pb)^b))^N }
S_of <- function(.fit, P) sum(.fit$cells$nu * P)                         # expected species (incl. unseen) for a cell-detection vector
taxa_of <- function(.fit, P) setNames(as.vector(.fit$post %*% P), rownames(.fit$post))   # per observed species
