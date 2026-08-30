#!/usr/bin/env python3
"""One CRTM calibration measurement, run as its own process.

The pre-calculation planning (see prompts/precalc_crtm.md) needs numbers this
project does not have yet: how much of a single-threaded CRTM run is fixed
`CRTM_Init`/`CRTM_Destroy` overhead vs. per-profile forward-model cost, how
wall time and peak RSS respond to `n_threads` and `max_profiles_per_chunk`,
how many full-resolution snapshots can run concurrently on one 256 GB node
before memory runs out, and whether the tbbmalloc segfault seen on the login
node (in CRTM's OpenMP region, `ODPS_Predictor_Define.f90`) reproduces on a
compute node or is defeatable with an allocator env var.

Each invocation does exactly ONE scenario and prints one `RESULT ...` line
(plus, for `scale`, one `RESULT scale ...` line per slice). Running one
scenario per process is deliberate: a segfault in a heavy case then only
loses that case, and `scripts/crtm_calibration.sbatch` (which drives this)
can keep going. It also lets the sbatch launch many of these at once for the
concurrency/packing measurement.

Not part of the importable package; a standalone measurement tool like
`scripts/crtm_thread_benchmark.sbatch`. Imports `bnf_wrf_eval.crtm` from
whatever `sys.path`/`--src` points at, so it can measure a worktree's
version of `crtm.py` without installing it.

Usage
-----
    crtm_calibration_probe.py MODE --wrf FILE [options]

MODE:
  scale     Single-thread cost breakdown. Reads the WRF state once, then
            runs `crtm._run_crtm_chunk` on contiguous row-slices of the
            sizes given by --rows-list (default "1,2,4,8,16,32"), all in
            this one process, timing each. A linear fit of run time vs.
            profile count (done by the sbatch, offline) gives intercept =
            per-chunk CRTM_Init+Destroy cost, slope = per-profile forward
            cost. --repeat N re-runs the whole list N times to probe the
            "3rd CRTM cycle segfaults" instability.
  fullres   One full `crtm.read_simulated_brightness_temperature` call
            (cache disabled) at --threads / --chunk, on --wrf. Prints wall
            time, the WRF-read vs. CRTM split, chunk count, peak RSS, and a
            sha256 + summary stats of the Tb field for bit-identity checks
            across thread/chunk settings.

Options:
  --wrf FILE            wrfout file (required)
  --threads N           n_threads for fullres (default 1)
  --chunk N             max_profiles_per_chunk for fullres (default 4000);
                        a huge value (e.g. 100000000) forces a single
                        whole-domain CRTM call ("NOCHUNK")
  --rows-list a,b,c     slice sizes (in grid rows) for scale mode
  --repeat N            scale mode: repeat the whole rows-list N times
  --row-start N         scale mode: first row of the slice band (default 180)
  --coeff-path DIR      override CRTM coefficientPath (e.g. a /dev/shm copy)
  --label STR           free-form tag echoed back in the RESULT line
  --src DIR             prepend to sys.path before importing bnf_wrf_eval
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time


def _peak_rss_gb() -> float:
    """Process high-water-mark RSS (VmHWM) in GiB, from /proc/self/status."""
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) / (1024.0 * 1024.0)
    except OSError:
        pass
    return float("nan")


def _tb_fingerprint(tb) -> str:
    import numpy as np

    finite = np.isfinite(tb)
    digest = hashlib.sha256(np.ascontiguousarray(tb)).hexdigest()[:16]
    return (
        f"sha={digest} nan={int((~finite).sum())} "
        f"min={np.nanmin(tb):.3f} max={np.nanmax(tb):.3f} mean={np.nanmean(tb):.3f}"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["scale", "fullres"])
    ap.add_argument("--wrf", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=4000)
    ap.add_argument("--rows-list", default="1,2,4,8,16,32")
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--row-start", type=int, default=180)
    ap.add_argument("--coeff-path", default=None)
    ap.add_argument("--label", default="")
    ap.add_argument("--src", default=None)
    args = ap.parse_args()

    if args.src:
        sys.path.insert(0, args.src)

    import numpy as np  # noqa: F401  (used via _tb_fingerprint / _peak_rss_gb import path)
    import bnf_wrf_eval.crtm as C

    tag = args.label or args.mode

    if args.mode == "scale":
        t0 = time.time()
        state = C._read_wrf_state(args.wrf)
        t_read = time.time() - t0
        ny, nx = state["lon"].shape
        nz = state["p"].shape[0]
        pyCRTM, profilesCreate = C._get_pycrtm()  # noqa: N806
        print(
            f"RESULT scale_meta label={tag} read_state_s={t_read:.2f} "
            f"ny={ny} nx={nx} nz={nz} peak_rss_gb={_peak_rss_gb():.2f}",
            flush=True,
        )
        rows_list = [int(x) for x in args.rows_list.split(",") if x]
        for rep in range(args.repeat):
            for rows in rows_list:
                end = min(args.row_start + rows, ny)
                sl = C._slice_state(state, slice(args.row_start, end), slice(None))
                nprof = sl["lon"].shape[0] * sl["lon"].shape[1]
                t0 = time.time()
                tb = C._run_crtm_chunk(
                    pyCRTM, profilesCreate, sl, C.SENSOR_ID, C.CHANNEL,
                    C.DEFAULT_SAT_LON, C.DEFAULT_SAT_HEIGHT, args.coeff_path, 1, None,
                )
                dt_ = time.time() - t0
                print(
                    f"RESULT scale label={tag} rep={rep} rows={rows} profiles={nprof} "
                    f"chunk_s={dt_:.3f} ms_per_prof={1000.0 * dt_ / nprof:.3f} "
                    f"peak_rss_gb={_peak_rss_gb():.2f} "
                    f"tb_mean={float(np.nanmean(tb)):.3f}",
                    flush=True,
                )
        return 0

    # fullres
    orig_read = C._read_wrf_state
    orig_chunk = C._run_crtm_chunk
    stats = {"read_s": 0.0, "n_chunks": 0, "crtm_s": 0.0}

    def timed_read(path):
        s = time.time()
        r = orig_read(path)
        stats["read_s"] += time.time() - s
        return r

    def timed_chunk(*a, **k):
        s = time.time()
        r = orig_chunk(*a, **k)
        stats["crtm_s"] += time.time() - s
        stats["n_chunks"] += 1
        return r

    C._read_wrf_state = timed_read
    C._run_crtm_chunk = timed_chunk

    wall0 = time.time()
    lon, lat, tb, vt = C.read_simulated_brightness_temperature(
        args.wrf,
        n_threads=args.threads,
        max_profiles_per_chunk=args.chunk,
        coefficient_path=args.coeff_path,
        cache_dir=None,
    )
    wall = time.time() - wall0

    print(
        f"RESULT fullres label={tag} threads={args.threads} chunk={args.chunk} "
        f"wall_s={wall:.1f} read_state_s={stats['read_s']:.1f} "
        f"crtm_s={stats['crtm_s']:.1f} n_chunks={stats['n_chunks']} "
        f"profiles={tb.size} peak_rss_gb={_peak_rss_gb():.2f} "
        f"snap_per_hr={3600.0 / wall:.1f} {_tb_fingerprint(tb)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
