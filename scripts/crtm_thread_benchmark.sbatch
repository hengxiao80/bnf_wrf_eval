#!/bin/bash
#SBATCH --job-name=crtm_threads
#SBATCH --account=cli120
#SBATCH --partition=batch_short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/wolf2/arm/cli120/scratch/hengxiao80/crtm_bench_%j.out
#
# Full-resolution CRTM thread-scaling benchmark on a compute node.
#
# Why a batch job rather than the login node: the login node is shared, and
# the earlier login-node benchmarks were contended by other work (including
# our own concurrent jobs), which both skews timings and is poor etiquette.
# A compute node gives dedicated cores plus far more memory headroom.
#
# `--cpus-per-task` is what SLURM actually reserves; the script below sweeps
# CRTM's own `n_threads` up to that number so the two stay consistent.
#
# MEASURED RESULTS (2026-08-29, job 1023444, cirrus1, full resolution,
# 213,561 profiles, 20250917 rund1). Every row verified bit-identical to the
# 1-thread reference -- CRTM's forward model is column-independent, so
# threading cannot change the answer:
#
#     1 thread    243.6 s   (reference)
#     8 threads    73.4 s   3.3x
#    32 threads    59.0 s   4.1x
#    64 threads    56.6 s   4.3x
#
# Recommendation: `--cpus-per-task=8` with `n_threads=8` for real work --
# it captures 3.3x of the 4.3x ceiling using 1/8th the cores. The sweep here
# keeps 64 only so it can measure the higher thread counts.
#
# The final "NOCHUNK" case (max_profiles_per_chunk huge, i.e. one CRTM call
# for the whole domain) SEGFAULTS -- exit 139, peaking at ~49 GB of the
# 250 GB requested, so NOT an out-of-memory failure. Chunking is therefore
# load-bearing for stability, not merely memory insurance, and since it
# costs only ~5% of runtime it should stay on. The case is kept in the sweep
# only to document that; expect the job to end in FAILED because of it.

set -u

REPO=/gpfs/wolf2/arm/cli120/proj-shared/hengxiao80/bnf_wrf_eval
# Use the MAIN checkout's .venv: pyCRTM is installed only there (it can only
# be built against the local CRTM install, so it is not a declared project
# dependency and `uv run` in a worktree would create an empty venv instead).
PY=$REPO/.venv/bin/python3
# Source tree to import from -- override with SRC=... to test a worktree.
SRC=${SRC:-$REPO/src}

module load intel/20.0.4 openmpi/4.1.4 netcdf-c/4.8.1 netcdf-fortran/4.5.3

# Keep the OpenMP runtime's own view consistent with what we reserved, and
# bind threads so they don't migrate between cores mid-run.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "host=$(hostname) cpus-per-task=${SLURM_CPUS_PER_TASK} src=${SRC}"

$PY - <<'PYEOF'
import os, sys, time
sys.path.insert(0, os.environ["SRC"])
import numpy as np
from bnf_wrf_eval import crtm

WRF = ("/gpfs/wolf2/arm/cli120/proj-shared/hengxiao80/bnf_wrf_eval/satoshi_testruns/"
       "20250917lassobnfwrfhrrr3/rund1/wrfout_d01_2025-09-17_20_00_00")

max_threads = int(os.environ["SLURM_CPUS_PER_TASK"])
# Sweep thread counts, and separately test whether skipping chunking
# entirely (possible with a compute node's memory) helps on top of threads.
ref = None
for n_threads, chunk, label in [
    (1,              4000,  "chunked  1thr"),
    (8,              4000,  "chunked  8thr"),
    (32,             4000,  "chunked 32thr"),
    (max_threads,    4000,  f"chunked {max_threads}thr"),
    (max_threads, 10**9,    f"NOCHUNK {max_threads}thr"),
]:
    if n_threads > max_threads:
        continue
    t0 = time.time()
    lon, lat, tb, vt = crtm.read_simulated_brightness_temperature(
        WRF, stride=1, cache_dir=None, n_threads=n_threads,
        max_profiles_per_chunk=chunk,
    )
    elapsed = time.time() - t0
    if ref is None:
        ref, note = tb, "(reference)"
    else:
        # Threading/chunking must never change the answer -- verify, don't assume.
        note = f"identical={np.array_equal(tb, ref)}"
    print(f"RESULT {label:20s} {elapsed:8.1f}s  profiles={tb.size}  {note}", flush=True)
PYEOF
