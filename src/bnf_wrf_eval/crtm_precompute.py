"""Batch pre-computation of CRTM-simulated brightness temperature.

Separates the expensive CRTM forward-model step (~1 min per wrfout time at
8 threads -- see `crtm.read_simulated_brightness_temperature`) from
plotting. Run this once to populate `crtm_cache/` for every wrfout time of
a set of runs; afterwards the notebooks and `plotting.plot_*_tb_*` helpers
just reload the cached `.npz` (~0.02 s) instead of invoking CRTM inline.
Pair with `require_cache=True` on the plotting side (a
`crtm.read_simulated_brightness_temperature` argument) to turn a missing
entry into a loud error rather than a silent minutes-long recompute.

Enumerates each run directory's own `wrfout_<domain>_*` files, parsing the
valid time straight from the filename, so it does not go through the
plotting layer's per-time file lookups. The flat (run, time) work list is
deterministically ordered, so `--task-index`/`--task-count` (or the
`SLURM_ARRAY_TASK_ID`/`SLURM_ARRAY_TASK_COUNT` a job array sets) can slice
it across array tasks -- strided, so each task gets an even mix of case
days and forecast hours regardless of how many tasks there are. See
`scripts/precompute_crtm_cache.sbatch`.

Runs as a module (deliberately no console-script entry, so a bare
`uv sync` can't drop it and no reinstall is needed):

    python -m bnf_wrf_eval.crtm_precompute --runs scripts/precompute_crtm_runs.txt \\
        --cache-dir crtm_cache --threads 8 --chunk 8000

Key options (see `_build_parser` for the full list):
  --runs FILE        newline-separated run directories ('#' comments ok)
  --run DIR          one run directory (repeatable; combine with/instead of --runs)
  --domain d01       wrfout domain tag
  --every-minutes N  keep only times whose minute-of-day is a multiple of N
                     (e.g. 60 for hourly); default: every wrfout file
  --start / --end    ISO datetimes bounding which wrfout times to include
  --threads N        CRTM OpenMP threads (default 8 -- the calibration sweet spot)
  --chunk N          max_profiles_per_chunk (default 8000)
  --recompute        recompute and overwrite existing cache entries
  --retry N          re-attempt a time that raised an exception, up to N extra
                     times (default 1). Does NOT cover a hard segfault -- that
                     kills the process; the sbatch re-runs the whole task for
                     that, relying on this module skipping already-cached times.
  --task-index K --task-count M   process only worklist[K::M]
  --list             print the resolved work list (index, run, time) and exit

STATUS (2026-08-30): new. The per-snapshot compute path it drives
(`crtm.read_simulated_brightness_temperature`) is the same one the
notebooks already use and is verified; this module only adds enumeration,
work-list slicing, skip/retry bookkeeping, and logging around it.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import sys
import time
from pathlib import Path

DEFAULT_DOMAIN = "d01"
DEFAULT_THREADS = 8
DEFAULT_CHUNK = 8000
DEFAULT_CACHE_DIR = "crtm_cache"


def _parse_wrfout_time(path: Path, domain: str) -> dt.datetime | None:
    """Valid time encoded in a `wrfout_<domain>_%Y-%m-%d_%H_%M_%S` filename,
    or None if `path`'s name does not match that pattern (so unrelated
    files in the directory are skipped rather than erroring)."""
    prefix = f"wrfout_{domain}_"
    name = path.name
    if not name.startswith(prefix):
        return None
    stamp = name[len(prefix):]
    try:
        return dt.datetime.strptime(stamp, "%Y-%m-%d_%H_%M_%S")
    except ValueError:
        return None


def list_wrfout_times(
    run_dir: str | Path,
    domain: str = DEFAULT_DOMAIN,
    start: dt.datetime | None = None,
    end: dt.datetime | None = None,
    every_minutes: int | None = None,
) -> list[Path]:
    """Sorted-by-time list of `wrfout_<domain>_*` files in `run_dir`,
    filtered to [`start`, `end`] (inclusive, either bound optional) and, if
    `every_minutes` is given, to times whose minute-of-day is a multiple of
    it (e.g. 60 -> on the hour, 15 -> every quarter hour)."""
    run_dir = Path(run_dir)
    dated: list[tuple[dt.datetime, Path]] = []
    for path in run_dir.iterdir():
        when = _parse_wrfout_time(path, domain)
        if when is None:
            continue
        if start is not None and when < start:
            continue
        if end is not None and when > end:
            continue
        if every_minutes is not None:
            minute_of_day = when.hour * 60 + when.minute
            if minute_of_day % every_minutes != 0:
                continue
        dated.append((when, path))
    dated.sort()
    return [path for _, path in dated]


def build_worklist(
    run_dirs: list[str | Path],
    domain: str = DEFAULT_DOMAIN,
    start: dt.datetime | None = None,
    end: dt.datetime | None = None,
    every_minutes: int | None = None,
) -> list[Path]:
    """Flat, deterministically ordered list of every wrfout file to process,
    across all `run_dirs` (ordered by run directory string, then by time
    within each). Stable ordering is what makes `--task-index`/`--task-count`
    slicing reproducible across array tasks."""
    worklist: list[Path] = []
    for run_dir in sorted(str(d) for d in run_dirs):
        files = list_wrfout_times(run_dir, domain, start, end, every_minutes)
        if not files:
            print(f"WARNING: no matching wrfout files in {run_dir}", file=sys.stderr)
        worklist.extend(files)
    return worklist


def precompute(
    worklist: list[Path],
    *,
    cache_dir: str | Path = DEFAULT_CACHE_DIR,
    n_threads: int = DEFAULT_THREADS,
    max_profiles_per_chunk: int = DEFAULT_CHUNK,
    recompute: bool = False,
    retry: int = 1,
    log=print,
) -> dict:
    """Run `crtm.read_simulated_brightness_temperature` for each file in
    `worklist`, writing/refreshing its `crtm_cache/` entry.

    Idempotent: a file whose cache entry already exists is skipped unless
    `recompute=True`, so re-running after a partial failure (including a
    segfault that killed a previous process) just fills the gaps. Per-file
    exceptions are caught and retried up to `retry` extra times, then
    logged and skipped so one bad file does not abort the batch. Returns a
    summary dict (`ok`, `skipped`, `failed`, `failures`, `seconds`)."""
    from . import crtm  # lazy: pulls in pyCRTM

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    summary = {"ok": 0, "skipped": 0, "failed": 0, "failures": [], "seconds": 0.0}
    t_start = time.time()
    total = len(worklist)

    for i, wrf_file in enumerate(worklist, 1):
        wrf_file = Path(wrf_file)
        tag = f"{wrf_file.parent.parent.name}/{wrf_file.parent.name}/{wrf_file.name}"
        cache_path = crtm._cache_path(
            wrf_file, cache_dir, crtm.SENSOR_ID, crtm.CHANNEL, 1
        )
        if cache_path.exists() and not recompute:
            summary["skipped"] += 1
            log(f"[{i}/{total}] SKIP  {tag}  (cached)")
            continue

        for attempt in range(retry + 1):
            t0 = time.time()
            try:
                _, _, values, valid_time = crtm.read_simulated_brightness_temperature(
                    wrf_file,
                    max_profiles_per_chunk=max_profiles_per_chunk,
                    cache_dir=cache_dir,
                    recompute=recompute,
                    n_threads=n_threads,
                )
                dt_s = time.time() - t0
                summary["ok"] += 1
                log(
                    f"[{i}/{total}] OK    {tag}  {valid_time:%Y-%m-%d_%H:%M}  "
                    f"{dt_s:.1f}s  mean={float(values.mean()):.2f}K"
                )
                break
            except Exception as exc:  # noqa: BLE001 -- one bad file must not abort the batch
                dt_s = time.time() - t0
                if attempt < retry:
                    log(
                        f"[{i}/{total}] RETRY {tag}  attempt {attempt + 1} failed "
                        f"after {dt_s:.1f}s: {exc!r}"
                    )
                    continue
                summary["failed"] += 1
                summary["failures"].append(str(wrf_file))
                log(
                    f"[{i}/{total}] FAIL  {tag}  after {attempt + 1} attempt(s): {exc!r}"
                )

    summary["seconds"] = time.time() - t_start
    return summary


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m bnf_wrf_eval.crtm_precompute",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--runs", type=Path, default=None,
                   help="file with one run directory per line ('#' comments allowed)")
    p.add_argument("--run", action="append", default=[], metavar="DIR",
                   help="a run directory (repeatable); combined with --runs")
    p.add_argument("--domain", default=DEFAULT_DOMAIN)
    p.add_argument("--cache-dir", type=Path, default=Path(DEFAULT_CACHE_DIR))
    p.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    p.add_argument("--chunk", type=int, default=DEFAULT_CHUNK,
                   help="max_profiles_per_chunk")
    p.add_argument("--every-minutes", type=int, default=None,
                   help="keep only times whose minute-of-day is a multiple of this")
    p.add_argument("--start", type=dt.datetime.fromisoformat, default=None)
    p.add_argument("--end", type=dt.datetime.fromisoformat, default=None)
    p.add_argument("--recompute", action="store_true")
    p.add_argument("--retry", type=int, default=1,
                   help="extra attempts per file on a Python-level exception")
    p.add_argument("--task-index", type=int, default=None,
                   help="process only worklist[task_index::task_count]")
    p.add_argument("--task-count", type=int, default=None)
    p.add_argument("--limit", type=int, default=None,
                   help="process at most this many files (after task slicing) -- for smoke tests")
    p.add_argument("--list", action="store_true",
                   help="print the resolved work list and exit")
    return p


def _resolve_run_dirs(args: argparse.Namespace) -> list[str]:
    run_dirs: list[str] = list(args.run)
    if args.runs is not None:
        for line in args.runs.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                run_dirs.append(line)
    if not run_dirs:
        raise SystemExit("no run directories given (use --runs and/or --run)")
    return run_dirs


def _resolve_task_slice(args: argparse.Namespace) -> tuple[int | None, int | None]:
    """--task-index/--task-count, falling back to the SLURM job-array
    environment so the sbatch does not have to pass them explicitly."""
    idx = args.task_index
    cnt = args.task_count
    if idx is None and "SLURM_ARRAY_TASK_ID" in os.environ:
        idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
    if cnt is None and "SLURM_ARRAY_TASK_COUNT" in os.environ:
        cnt = int(os.environ["SLURM_ARRAY_TASK_COUNT"])
    if (idx is None) != (cnt is None):
        raise SystemExit("task-index and task-count must be given together")
    return idx, cnt


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    run_dirs = _resolve_run_dirs(args)

    worklist = build_worklist(
        run_dirs, args.domain, args.start, args.end, args.every_minutes
    )
    idx, cnt = _resolve_task_slice(args)
    if idx is not None:
        worklist = worklist[idx::cnt]
        scope = f" (task {idx}/{cnt})"
    else:
        scope = ""
    if args.limit is not None:
        worklist = worklist[: args.limit]
        scope += f" (limit {args.limit})"

    if args.list:
        for j, path in enumerate(worklist):
            when = _parse_wrfout_time(Path(path), args.domain)
            print(f"{j:5d}  {when:%Y-%m-%d_%H:%M}  {path}")
        print(f"# {len(worklist)} file(s){scope}")
        return 0

    print(
        f"crtm_precompute: {len(worklist)} file(s){scope}  threads={args.threads}  "
        f"chunk={args.chunk}  cache_dir={args.cache_dir}  recompute={args.recompute}",
        flush=True,
    )
    summary = precompute(
        worklist,
        cache_dir=args.cache_dir,
        n_threads=args.threads,
        max_profiles_per_chunk=args.chunk,
        recompute=args.recompute,
        retry=args.retry,
    )
    print(
        f"crtm_precompute done{scope}: ok={summary['ok']} skipped={summary['skipped']} "
        f"failed={summary['failed']}  {summary['seconds'] / 60:.1f} min",
        flush=True,
    )
    for failed in summary["failures"]:
        print(f"  FAILED: {failed}", flush=True)
    return 1 if summary["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
