#!/bin/bash
# run_r12c.sh -- rebuild frame 2 on release r12c (Dryad deposit dryad_2026-10-24; Kevin 2026-10-25 "yes, do it"), screened calls
# (common.R SCREEN=1 default). Same methods and run sequence as run_r2.sh, with the unseen-fit merge step that run_r2.sh lacked.
# 01_build.R was run first by hand (COSTTAG=). Script 22 reads the hop-flagged reads rebuilt from the deposit (HOPDIR).
# MEMORY GATE (coordinator 2026-10-25): every job reserves its expected peak (r2 peak x 1.35, for the larger COI set). A job starts
# only when (a) my reservations + it <= CAP_GB (3), (b) all R on the machine (others' live RSS + my reservations + it) <= ALLR_GB (6),
# (c) available memory >= AVAIL_GB (4). Otherwise it waits, and logs "WAIT <reason>" once per wait. Two lanes share the gate.
cd "/Users/lafferty/My Drive/Rscripts/eDNA/CaliforniaIntertidal/studies_2026-10-10/sampling_opt/frame2" || exit 1
unset COSTS COSTTAG FITS GRIDTAG RULEVER A FTYPE VG BUDGET_AWARE TRIGGER GRID LEFTOVER MARKER NCORE TEST_B SCREEN CI_DEPOSIT CI_LABEL
O=out_r12c; D=$O/costs_hours_2026-10-11; T=$O/costs_hours_2026-10-11_template; mkdir -p $D $T
HOPDIR=/private/tmp/claude-501/-Users-lafferty/3861b05b-4942-4975-afc5-5207b58e09cb/scratchpad/hop_out_r12c/
RES=/private/tmp/claude-501/-Users-lafferty/3861b05b-4942-4975-afc5-5207b58e09cb/scratchpad/r12c_res; mkdir -p $RES
CAP_GB=3; ALLR_GB=6; AVAIL_GB=4
st() { echo "$(date +%H:%M:%S) $*" >> $O/run_r12c_stages.log; }
avail_gb() { vm_stat | awk '/Pages free/ {f=$3} /Pages inactive/ {i=$3} /Pages speculative/ {s=$3} END {printf "%.2f", (f+i+s)*16384/1073741824}'; }
mine_res() { cat $RES/*.res 2>/dev/null | awk '{s+=$1} END {printf "%.2f", s+0}'; }
# live RSS of R processes that are NOT descendants of this driver (forked workers of my jobs count as mine)
DRIVER=$$
others_gb() { ps -ax -o pid=,ppid=,rss=,comm= | awk -v d=$DRIVER '{ pp[$1]=$2; rs[$1]=$3; cm[$1]=$4 }
  END { for (p in pp) if (cm[p] ~ /\/R$|^R$/) { q=p; mine=0; for (i=0; i<50 && q>1; i++) { if (q==d) { mine=1; break }; q=pp[q] }
        if (!mine) s+=rs[p] }; printf "%.2f", s/1048576 }'; }
gate() { local est=$1 name=$2 why="" last=""
  while true; do
    if mkdir $RES/.lock 2>/dev/null; then
      local m o a; m=$(mine_res); o=$(others_gb); a=$(avail_gb)
      why=$(awk -v m=$m -v o=$o -v a=$a -v e=$est -v c=$CAP_GB -v r=$ALLR_GB -v v=$AVAIL_GB 'BEGIN {
        if (m > 0 && m + e > c) print "my tree " m "+" e " GB > " c; else if (o + m + e > r) print "all R " o "+" m "+" e " GB > " r;
        else if (a < v) print "available " a " GB < " v; else print "" }')
      if [ -z "$why" ]; then echo "$est" > $RES/$name.res; rmdir $RES/.lock; return 0; fi
      rmdir $RES/.lock
      [ "$why" != "$last" ] && case "$why" in available*|all*) st "WAIT $name: $why (PAUSE: tell the coordinator)";; *) st "WAIT $name: $why";; esac
      last=$why
    fi
    sleep 30
  done }
# usage: J est_gb name log env... Rscript x.R
J() { local est=$1 name=$2 log=$3; shift 3; gate $est $name
  /usr/bin/time -l "$@" > "$log" 2>&1 || st "FAILED $log"; mv $RES/$name.res $RES/$name.res.done; }
st start
# ---- stage 1: fits
( J 2.10 fit02 $O/log_02_fit.txt env COSTTAG= Rscript 02_fit.R ) &
( J 0.30 s16 $O/log_16.txt Rscript 16_pooling_penalty.R; J 0.80 s18 $O/log_18.txt Rscript 18_a_from_singles.R
  J 1.95 fit02a $O/log_02_fit_a0.0833.txt env COSTTAG= A=0.0833 Rscript 02_fit.R ) &
wait; st "stage 1 done"
# ---- stage 2: two lanes through the gate
( J 2.10 s23 $O/log_23.txt Rscript 23_template_fit.R; st "template fit done"
  J 1.40 t03 $T/log_03.txt env FITS=fits_template.rds A=0.03 COSTTAG=hours_2026-10-11_template Rscript 03_designs.R
  J 0.60 t25 $T/log_25_v2.txt env FITS=fits_template.rds COSTTAG=hours_2026-10-11_template Rscript 25_rule_designs.R
  J 0.35 t27 $T/log_27_v2.txt env FITS=fits_template.rds COSTTAG=hours_2026-10-11_template Rscript 27_rule_check_fine.R
  J 1.30 t24 $O/log_24_template.txt env FITS=fits_template.rds Rscript 24_em_convergence.R
  J 0.45 t26 $O/log_26_template.txt env FITS=fits_template.rds Rscript 26_location_curve.R
  st "template designs done; 29 starting"
  J 2.50 t29 $T/log_29.txt env COSTTAG=hours_2026-10-11_template NCORE=7 Rscript 29_rule_visit_bands.R; st "lane A done" ) &
( J 1.55 d03 $D/log_03_disc.txt env FTYPE=disc Rscript 03_designs.R; J 0.70 d13 $D/log_13.txt Rscript 13_extract_cap33.R
  J 0.65 d14 $D/log_14.txt Rscript 14_site_knowledge.R; J 1.10 d07 $D/log_07.txt env A=0.03 Rscript 07_unseen_sensitivity.R
  J 0.50 d04 $D/log_04.txt Rscript 04_common_species.R; J 0.60 d25 $D/log_25_v2.txt Rscript 25_rule_designs.R
  J 0.35 d27 $D/log_27_v2.txt Rscript 27_rule_check_fine.R
  J 0.90 d15a $D/log_15_12S.txt env MARKER=12S Rscript 15_unseen_profile.R; J 1.85 d15b $D/log_15_COI.txt env MARKER=COI Rscript 15_unseen_profile.R
  J 0.50 merge $O/log_unseen_merge.txt Rscript -e 'F="out_r12c/"; for (e in c("lo","hi")) { x=c(readRDS(paste0(F,"fits_unseen_",e,"_12S.rds")), readRDS(paste0(F,"fits_unseen_",e,"_COI.rds"))); p=paste0(F,"fits_unseen_",e,".rds"); stopifnot(!file.exists(p)); saveRDS(x,p) }'
  for e in lo hi; do mkdir -p $O/costs_hours_2026-10-11_unseen_$e
    J 1.35 u03$e $O/costs_hours_2026-10-11_unseen_$e/log_03.txt env FITS=fits_unseen_$e.rds COSTTAG=hours_2026-10-11_unseen_$e A=0.03 FTYPE=disc Rscript 03_designs.R; done
  J 0.75 d17 $D/log_17.txt Rscript 17_unseen_profile_designs.R
  J 1.45 c03 $D/log_03_capsule.txt env A=0.03 FTYPE=capsule Rscript 03_designs.R; J 1.35 c24 $O/log_24.txt Rscript 24_em_convergence.R
  J 0.45 c26 $O/log_26.txt Rscript 26_location_curve.R
  for t in 3 6 10; do J 0.70 w$t $D/log_11_wide_fill_t$t.txt env BUDGET_AWARE=1 GRID=wide TRIGGER=$t Rscript 11_enumerate_small_levels.R
    J 0.70 wu$t $D/log_11_wide_unspent_t$t.txt env BUDGET_AWARE=1 GRID=wide LEFTOVER=unspent TRIGGER=$t Rscript 11_enumerate_small_levels.R; done
  J 0.50 d05 $D/log_05.txt Rscript 05_export_rule_inputs.R; J 0.50 d10 $D/log_10.txt Rscript 10_branch_test.R; J 0.50 d09 $D/log_09.txt Rscript 09_rule_vs_grid.R
  J 0.30 d11 $D/log_11.txt Rscript 11_enumerate_small_levels.R; J 0.30 d11b $D/log_11_ba.txt env BUDGET_AWARE=1 Rscript 11_enumerate_small_levels.R
  J 0.30 d11c $D/log_11_t6.txt env BUDGET_AWARE=1 TRIGGER=6 Rscript 11_enumerate_small_levels.R; J 0.30 d11d $D/log_11_t10.txt env BUDGET_AWARE=1 TRIGGER=10 Rscript 11_enumerate_small_levels.R
  J 1.75 s19 $O/log_19.txt Rscript 19_a_within_run_profile.R; J 0.50 s20 $O/log_20.txt Rscript 20_turnover.R; J 0.40 s21 $O/log_21.txt Rscript 21_turnover_checks.R
  J 1.70 s28 $O/log_28.txt Rscript 28_template_sensitivity.R
  J 0.65 s22 $O/log_22_v2.txt env HOPDIR=$HOPDIR Rscript 22_single_reaction_detections.R; st "lane B done" ) &
wait; st "all done"
