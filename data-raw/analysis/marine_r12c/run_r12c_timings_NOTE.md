# r12c frame-2 run: timings and the memory-gate queue (2026-10-25)

`run_r12c_timings.csv` has one row per job of `../run_r12c.sh`:
- lane;
- reservation (GB);
- start and end (the log file's birth and modification times);
- queued_min: time since the previous job in the lane ended, i.e. time spent waiting at the gate;
- wall time and peak RSS (from `/usr/bin/time -l`);
- the same job's r2 wall time and peak (`../out_r2/`), and the r12c/r2 wall ratio.

Peak RSS is that of the job's main process. For script 29 it excludes the 7 forked workers (see below).

The run: 12:22:02 to 18:33:26 PDT (6 h 11 min). No job failed. No PAUSE waits: the all-R and available-memory gates never bound. All R on the machine peaked at 3.13 GB (15:55). Available memory was never below 13.9 GB, and swap was unchanged.

## Script 29 (visit bands, template, NCORE=7)
- Distinct 12S+COI inner proposals: 1,649 on r12c vs 1,275 on r2 (×1.29). The larger screened COI set gives more distinct proposals.
- Coarse pass (71 budgets): 26.6 min vs 18.1 min (×1.47). The refine chunks reached about 4.2 min each, against about 3.4 min at the end on r2.
- Total wall: 149.6 min vs 123.3 min on r2 (×1.21).
- Memory of its tree (parent plus workers, from the memwatch log, less about 1.1 GB of other sessions' R): about 1.0–1.3 GB in steady state, peaking near 1.9 GB.
- It was reserved at 2.5 GB (r2 tree peak 1.81 GB × 1.35).

## The queue
Under the 3 GB cap on my tree, two jobs queued for long:
- **Script 29 waited 30 min** (15:04:56 → 15:34:59) for lane B's c24 (1.35 GB) to finish.
- **Lane B then waited 145 min** (w3, 15:39:27 → 18:04:41). Every later lane B job needed ≥ 0.70 GB with 29's 2.5 GB held. Once 29 finished, all of lane B's remaining jobs (about 29 min of work) ran in 29 min.
- Other queues: fit02a 53 min (behind fit02), d03 42 min (behind s23), d15b 7.5 min.

So about 2.4 h of the 6.2 h run was lane B idling behind script 29's reservation.

## For the next larger run
- **Reserve 29 at its measured tree peak** (about 2.0 GB), not 2.5. Lane B's 0.30–0.70 GB jobs (scripts 11, 05, 10, 09, 20, 21, 22) then fit beside it under 3 GB.
- **Or order lane B so its small jobs run while 29 runs,** with 19 and 28 (≥ 1.7 GB) after 29.
- **Script 29's cost grows faster than the proposal count** (×1.29 proposals → ×1.47 coarse pass). Scale its time estimate by roughly (proposals ratio)^1.5.
- **One reservation was exceeded:** s28 peaked at 1.83 GB on a 1.70 reservation (harmless: it ran alone). Raise s28 to about 2.0.
- Scripts 19, 28 and 15_COI peaked about 30% above r2; the ×1.35 margin covered them.
