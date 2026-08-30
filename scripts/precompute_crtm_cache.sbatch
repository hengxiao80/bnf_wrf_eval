#!/bin/bash
#SBATCH --job-name=crtm_precalc
#SBATCH --account=cli120
#SBATCH --partition=batch_all
#SBATCH --array=0-127%32
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/wolf2/arm/cli120/scratch/hengxiao80/crtm_precalc_%A_%a.out
#
# Pre-compute CRTM-simulated ABI ch.13 brightness temperature for every
# wrfout time of the d1 HRRR3 runs, filling crtm_cache/ so the notebooks
# and plotting helpers just reload .npz files.
#
# WHY THESE NUMBERS (from scripts/crtm_calibration.sbatch, job 1023455 on a
# 256 GB node, 2026-08-30):
#   * 8 threads/task: ~2.9x wall speedup over 1 thread (crtm-only ~3.1x);
#     returns diminish past 4. chunk size barely matters for speed, so
#     8000 is chosen for a good memory/one-CRTM_Init balance (~3.5 GB peak).
#   * 16 tasks/node (128 cores / 8): --cpus-per-task=8 + --mem=14G lets
#     SLURM pack 16 per node (16 x 14 = 224 GB < 256). Peak RSS measured
#     ~3.5 GB, so 14 GB is generous headroom.
#   * %32 throttle = 32 tasks at once = 2 nodes. Under packing the real
#     limiter is GPFS contention on the 3.4 GB wrfout reads (per-snapshot
#     wall rose from ~55 s solo to ~157 s at 16-wide), not CPU or CRTM.
#     Measured throughput ~366 snapshots/hr/node, so 2601 snapshots take
#     roughly:  %32 (2 nodes) ~3.5 h  |  %64 (4 nodes) ~1.8 h  |
#     %128 (8 nodes) ~55 min. Pick the throttle by queue etiquette/urgency;
#     the batch is restartable so a low throttle costs only wall time.
#   * --array=0-127: 128 array tasks; task k processes worklist[k::128]
#     (strided, so each gets an even mix of case days and forecast hours).
#     MUST stay 0-based. Bump the upper bound only together with the
#     stride math below.
#   * Idempotent: each task skips wrfout times already in crtm_cache/, so a
#     re-run (after a timeout, a node failure, or the segfault guard below)
#     just fills the gaps. Safe to resubmit the whole array.
#   * Segfault guard: a rare tbbmalloc SIGSEGV was seen in calibration, but
#     only on very small CRTM calls (<~2000 profiles); every real
#     full-resolution run (chunks >= 2000) was stable across ~55 processes.
#     A hard segfault kills the process, so --retry inside the module can't
#     catch it -- instead this script re-runs the task up to 3 times, and
#     the module's skip-if-cached behaviour means each re-run resumes.
#
# Submit (from the main checkout):
#   sbatch scripts/precompute_crtm_cache.sbatch
# Dry run first (no CRTM, just prints the work list for one task):
#   SRC=$PWD/.claude/worktrees/crtm-precalc/src \
#     .venv/bin/python3 -m bnf_wrf_eval.crtm_precompute \
#     --runs scripts/precompute_crtm_runs.txt --list --task-index 0 --task-count 128
# Calibrate/run a worktree's code instead of the installed src:
#   SRC=/path/to/worktree/src RUNS=/path/to/worktree/scripts/precompute_crtm_runs.txt \
#     sbatch /path/to/worktree/scripts/precompute_crtm_cache.sbatch

set -u

REPO=/gpfs/wolf2/arm/cli120/proj-shared/hengxiao80/bnf_wrf_eval
PY=$REPO/.venv/bin/python3
# Source tree to import bnf_wrf_eval from. Default is the editable install's
# own path (the main checkout's src); override SRC to test a worktree.
SRC=${SRC:-$REPO/src}
RUNS=${RUNS:-$REPO/scripts/precompute_crtm_runs.txt}
CACHE_DIR=${CACHE_DIR:-$REPO/crtm_cache}

THREADS=${THREADS:-8}
CHUNK=${CHUNK:-8000}
# Optional extra filtering, passed straight through if set:
#   EVERY_MINUTES=60  -> hourly only
EXTRA=${EXTRA:-}
[ -n "${EVERY_MINUTES:-}" ] && EXTRA="$EXTRA --every-minutes ${EVERY_MINUTES}"

module load intel/20.0.4 openmpi/4.1.4 netcdf-c/4.8.1 netcdf-fortran/4.5.3

export PYTHONPATH="$SRC${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-$THREADS}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

cd "$REPO"

ARRAY_COUNT=$(( ${SLURM_ARRAY_TASK_MAX:-0} - ${SLURM_ARRAY_TASK_MIN:-0} + 1 ))
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

echo "host=$(hostname) date=$(date -Is) task=${TASK_ID}/${ARRAY_COUNT} src=${SRC}"
echo "threads=${THREADS} chunk=${CHUNK} cache_dir=${CACHE_DIR} extra='${EXTRA}'"

rc=1
for attempt in 1 2 3; do
    "$PY" -m bnf_wrf_eval.crtm_precompute \
        --runs "$RUNS" \
        --cache-dir "$CACHE_DIR" \
        --threads "$THREADS" \
        --chunk "$CHUNK" \
        --retry 1 \
        --task-index "$TASK_ID" \
        --task-count "$ARRAY_COUNT" \
        $EXTRA
    rc=$?
    [ "$rc" -eq 0 ] && break
    echo "task ${TASK_ID} attempt ${attempt} exited ${rc}; re-running (cached times will be skipped)"
    sleep 5
done

echo "task ${TASK_ID} final rc=${rc} ($(date -Is))"
exit "$rc"
