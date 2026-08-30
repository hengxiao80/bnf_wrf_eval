#!/bin/bash
#SBATCH --job-name=crtm_calib
#SBATCH --account=cli120
#SBATCH --partition=batch_short
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/wolf2/arm/cli120/scratch/hengxiao80/crtm_calib_%j.out
#
# CRTM pre-calculation calibration, for a STANDARD 256 GB compute node.
# ---------------------------------------------------------------------
# Drives scripts/crtm_calibration_probe.py to measure what the pre-calc
# plan (prompts/precalc_crtm.md) needs before committing to a job-array
# layout:
#
#   A. Single-thread cost breakdown: per-chunk CRTM_Init+Destroy overhead
#      vs. per-profile forward-model cost (linear fit of chunk time vs.
#      profile count, done offline from the "RESULT scale" lines).
#   A2. Repeated-cycle stability: does the login-node tbbmalloc segfault
#      (CRTM's OpenMP region, ODPS_Predictor_Define.f90) reproduce on a
#      compute node when many CRTM_Init/Destroy cycles run in one process?
#   B. n_threads x max_profiles_per_chunk sweep at full resolution: wall
#      time, WRF-read vs. CRTM split, peak RSS, and a Tb sha256 so every
#      setting is checked bit-identical to the threads=1/chunk=4000
#      reference (threading and chunking must not change the answer).
#   C. NOCHUNK + allocator env vars: can a single whole-domain CRTM call
#      be made to run (it segfaults by default -- see crtm.py), e.g. with
#      KMP_DISABLE_FAST_MEMORY_ALLOCATOR / MALLOC_ARENA_MAX? If so, that
#      is one CRTM_Init per snapshot instead of ~56.
#   D. Per-microphysics-scheme cost: Thompson vs. P3 mp52 (2 ice cats,
#      heaviest) vs. P3 mp53, since profile build and cloud handling
#      differ.
#   E. Concurrency / packing: run N independent full-resolution snapshots
#      at once with threads set so N*threads = 128, for
#      (N,threads) = (16,8) (32,4) (64,2) (128,1). Reports aggregate
#      throughput (snapshots/hr/node) and peak node memory -- the number
#      that decides how many array tasks to pack per node. The 128-wide
#      wave is expected to hit the 256 GB limit; that is a result, not a
#      failure.
#
# Sections run cheap-and-informative first (A, B) so a time-out or node
# problem still leaves the core numbers. Expect the job to end FAILED if
# section C segfaults or section E OOMs -- read the RESULT/CASE lines, not
# the job state.
#
# Submit (from the main checkout):
#   sbatch scripts/crtm_calibration.sbatch
# To calibrate a worktree's crtm.py instead:
#   SRC=/path/to/worktree/src PROBE=/path/to/worktree/scripts/crtm_calibration_probe.py \
#       sbatch /path/to/worktree/scripts/crtm_calibration.sbatch
#
# A companion sbatch for the 512 GB (batch_high_memory) and 1 TB (gpu)
# nodes is intentionally not written yet -- calibrate the standard node
# first, then decide if the others are worth it.

set -u

REPO=/gpfs/wolf2/arm/cli120/proj-shared/hengxiao80/bnf_wrf_eval
# Use the MAIN checkout's .venv: pyCRTM is installed only there (it can only
# be built against the local CRTM install, so it is not a declared project
# dependency and `uv run` in a worktree would create an empty venv instead).
PY=$REPO/.venv/bin/python3
SRC=${SRC:-$REPO/src}
PROBE=${PROBE:-$REPO/scripts/crtm_calibration_probe.py}

# Representative wrfout files, one per microphysics scheme this project's
# d1 runs use (2025-09-17 case, 20:00 UTC -- an active convective time).
TR=$REPO/satoshi_testruns/20250917lassobnfwrfhrrr3
WRF_THOMPSON=$TR/rund1/wrfout_d01_2025-09-17_20_00_00
WRF_P3_MP52=$TR/rund1-mp52/wrfout_d01_2025-09-17_20_00_00
WRF_P3_MP53=$TR/rund1-mp53/wrfout_d01_2025-09-17_20_00_00

module load intel/20.0.4 openmpi/4.1.4 netcdf-c/4.8.1 netcdf-fortran/4.5.3

export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "host=$(hostname) date=$(date -Is) src=$SRC probe=$PROBE"
echo "node memory:"; free -g

# run_probe <label> <n_threads_for_omp_env> <extra_env> -- <probe.py args...>
# <extra_env> is a space-separated "VAR=val VAR=val" string (use "" for
# none), passed to the child via env(1) so it works regardless of export
# state. Runs one probe under /usr/bin/time -v, echoes a CASE line with
# return code, wall seconds and max RSS (KB) from time -v, and never
# aborts the job on a probe crash.
run_probe() {
    local label=$1 omp=$2 extra_env=$3; shift 3
    [ "$1" = "--" ] && shift
    local tf; tf=$(mktemp)
    local t0 t1 rc
    t0=$(date +%s.%N)
    env OMP_NUM_THREADS="$omp" $extra_env \
        /usr/bin/time -v -o "$tf" "$PY" "$PROBE" "$@" --src "$SRC" --label "$label"
    rc=$?
    t1=$(date +%s.%N)
    local maxrss; maxrss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tf")
    echo "CASE label=$label rc=$rc extra_env='${extra_env}'" \
         "wall_s=$(awk "BEGIN{printf \"%.1f\", $t1-$t0}") time_v_maxrss_kb=${maxrss:-NA}"
    rm -f "$tf"
    return 0
}

echo; echo "======== A. single-thread cost breakdown (Thompson) ========"
run_probe A_scale 1 "" -- scale --wrf "$WRF_THOMPSON" \
    --rows-list 1,2,3,4,6,8,12,16,24,32,48

echo; echo "======== A2. repeated-cycle stability (does it segfault?) ========"
run_probe A2_repeat5 1 "" -- scale --wrf "$WRF_THOMPSON" \
    --rows-list 4,8,16 --repeat 5

echo; echo "======== B. n_threads x chunk sweep (Thompson, full resolution) ========"
for chunk in 2000 4000 8000 16000; do
    for thr in 1 2 4 8; do
        run_probe "B_t${thr}_c${chunk}" "$thr" "" -- fullres --wrf "$WRF_THOMPSON" \
            --threads "$thr" --chunk "$chunk"
    done
done

echo; echo "======== D. per-scheme cost (threads=8, chunk=8000) ========"
run_probe D_thompson 8 "" -- fullres --wrf "$WRF_THOMPSON" --threads 8 --chunk 8000
run_probe D_p3_mp52  8 "" -- fullres --wrf "$WRF_P3_MP52"  --threads 8 --chunk 8000
run_probe D_p3_mp53  8 "" -- fullres --wrf "$WRF_P3_MP53"  --threads 8 --chunk 8000

echo; echo "======== C. NOCHUNK + allocator variants (threads=8, Thompson) ========"
run_probe C1_c64000          8 "" -- fullres --wrf "$WRF_THOMPSON" --threads 8 --chunk 64000
run_probe C2_nochunk_default 8 "" -- fullres --wrf "$WRF_THOMPSON" --threads 8 --chunk 100000000
run_probe C3_nochunk_kmpnofast 8 "KMP_DISABLE_FAST_MEMORY_ALLOCATOR=1" \
    -- fullres --wrf "$WRF_THOMPSON" --threads 8 --chunk 100000000
run_probe C4_nochunk_arena1  8 "MALLOC_ARENA_MAX=1" \
    -- fullres --wrf "$WRF_THOMPSON" --threads 8 --chunk 100000000
run_probe C5_nochunk_both    8 "KMP_DISABLE_FAST_MEMORY_ALLOCATOR=1 MALLOC_ARENA_MAX=2" \
    -- fullres --wrf "$WRF_THOMPSON" --threads 8 --chunk 100000000

echo; echo "======== E. concurrency / packing (Thompson, chunk=8000) ========"
for pair in "16 8" "32 4" "64 2" "128 1"; do
    set -- $pair
    N=$1; THR=$2
    echo "-- E wave: N=$N x threads=$THR  (N*threads=$((N*THR)))"
    export OMP_NUM_THREADS=$THR
    # background node-memory sampler for this wave
    ( for i in $(seq 1 120); do echo "MEM $(date +%s) $(awk '/MemAvailable/{print $2}' /proc/meminfo)"; sleep 5; done ) &
    SAMPLER=$!
    declare -a PIDS=()
    wave0=$(date +%s.%N)
    for k in $(seq 0 $((N-1))); do
        c0=$((k*THR)); c1=$((c0+THR-1))
        taskset -c "${c0}-${c1}" "$PY" "$PROBE" fullres --wrf "$WRF_THOMPSON" \
            --threads "$THR" --chunk 8000 --src "$SRC" --label "E_N${N}_t${THR}_k${k}" &
        PIDS+=($!)
    done
    fail=0
    for p in "${PIDS[@]}"; do wait "$p" || fail=$((fail+1)); done
    wave1=$(date +%s.%N)
    kill "$SAMPLER" 2>/dev/null
    echo "CASE label=E_N${N}_t${THR} wave_wall_s=$(awk "BEGIN{printf \"%.1f\", $wave1-$wave0}")" \
         "n_procs=$N n_failed=$fail throughput_snap_per_hr=$(awk "BEGIN{printf \"%.1f\", $N*3600.0/($wave1-$wave0)}")"
    free -g
done

echo; echo "calibration done: $(date -Is)"
