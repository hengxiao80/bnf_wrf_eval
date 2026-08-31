#!/usr/bin/env python3
"""Bulk-download the GOES ABI ch.13 and MRMS composite-reflectivity files
the plots need, for the three d1 HRRR3 case windows, into ``goes_data/``
and ``mrms_data/``.

`bnf_wrf_eval.plotting`'s wrappers can each fetch one observation file on
demand (``auto_download_*``), but a job array of hundreds of plot tasks all
hitting S3 -- and racing on the same missing path -- is worth avoiding.
Run this once first; then run ``batch_plot`` without ``--auto-download``.

    # hourly only -- enough for the 4-panel Tb / 3-panel refl comparisons
    python scripts/download_obs.py --hourly

    # full cadence -- GOES native ~5 min, MRMS every 15 min -- for the
    # single-panel movies
    python scripts/download_obs.py --which goes,mrms --cadence-mrms 15

Anonymous S3 (NOAA Open Data), no credentials. Idempotent: files already
present are skipped.
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3
from botocore import UNSIGNED
from botocore.config import Config

REPO = Path(__file__).resolve().parent.parent
RUNS_ROOT = REPO / "satoshi_testruns"
DOMAIN = "d01"
CASES = [
    "20250502lassobnfwrfhrrr3",
    "20250520lassobnfwrfhrrr3",
    "20250917lassobnfwrfhrrr3",
]
REF_SCHEME = "rund1"

GOES_BUCKET = "noaa-goes19"
GOES_PRODUCT = "ABI-L2-CMIPC"
GOES_CHANNEL = 13
GOES_SAT = "G19"

MRMS_BUCKET = "noaa-mrms-pds"
MRMS_PRODUCT = "MergedReflectivityQCComposite_00.50"


def _s3():
    return boto3.client("s3", config=Config(signature_version=UNSIGNED))


def _case_span(case: str) -> tuple[dt.datetime, dt.datetime]:
    ref_dir = RUNS_ROOT / case / REF_SCHEME
    prefix = f"wrfout_{DOMAIN}_"
    stamps = []
    for p in ref_dir.iterdir():
        if not p.name.startswith(prefix):
            continue
        try:
            stamps.append(dt.datetime.strptime(p.name[len(prefix):], "%Y-%m-%d_%H_%M_%S"))
        except ValueError:
            pass
    if not stamps:
        raise SystemExit(f"no wrfout_{DOMAIN}_* files under {ref_dir}")
    return min(stamps), max(stamps)


def _goes_scan_time(key: str) -> dt.datetime:
    s_field = key.split("_s")[1].split("_")[0]  # YYYYJJJHHMMSSs
    return dt.datetime.strptime(s_field[:-1], "%Y%j%H%M%S")


def _mrms_scan_time(key: str) -> dt.datetime:
    stem = Path(key).name.removesuffix(".grib2.gz")
    return dt.datetime.strptime(stem.rsplit("_", 1)[1], "%Y%m%d-%H%M%S")


def list_goes_keys(s3, start: dt.datetime, end: dt.datetime) -> list[str]:
    """Every ABI ch.13 CMIP object whose scan start is in [start, end]
    (one S3 list per hour; the product's native cadence is ~5 min)."""
    # List one hour past `end` too, and accept scans up to 10 min past it, so
    # the closing on-the-hour frame (whose nearest scan is ~HH:01) isn't lost.
    keys: list[str] = []
    cutoff = end + dt.timedelta(minutes=10)
    hour = start.replace(minute=0, second=0, microsecond=0)
    last = end.replace(minute=0, second=0, microsecond=0) + dt.timedelta(hours=1)
    while hour <= last:
        prefix = (
            f"{GOES_PRODUCT}/{hour:%Y}/{hour:%j}/{hour:%H}/"
            f"OR_{GOES_PRODUCT}-M6C{GOES_CHANNEL:02d}_{GOES_SAT}_"
        )
        for page in s3.get_paginator("list_objects_v2").paginate(Bucket=GOES_BUCKET, Prefix=prefix):
            for obj in page.get("Contents", []):
                if start <= _goes_scan_time(obj["Key"]) <= cutoff:
                    keys.append(obj["Key"])
        hour += dt.timedelta(hours=1)
    return sorted(set(keys))


def list_mrms_keys(
    s3, start: dt.datetime, end: dt.datetime, cadence_minutes: int
) -> list[str]:
    """One MRMS object per `cadence_minutes` mark in [start, end] -- the
    scan closest to each mark (native cadence is ~2 min)."""
    by_time: dict[dt.datetime, str] = {}
    day = start.date()
    while day <= end.date():
        prefix = f"CONUS/{MRMS_PRODUCT}/{day:%Y%m%d}/MRMS_{MRMS_PRODUCT}_{day:%Y%m%d}-"
        for page in s3.get_paginator("list_objects_v2").paginate(Bucket=MRMS_BUCKET, Prefix=prefix):
            for obj in page.get("Contents", []):
                by_time[_mrms_scan_time(obj["Key"])] = obj["Key"]
        day += dt.timedelta(days=1)
    if not by_time:
        return []

    scans = sorted(by_time)
    chosen: list[str] = []
    mark = start
    step = dt.timedelta(minutes=cadence_minutes)
    while mark <= end:
        nearest = min(scans, key=lambda t: abs((t - mark).total_seconds()))
        if abs((nearest - mark).total_seconds()) <= 300:
            chosen.append(by_time[nearest])
        mark += step
    return sorted(set(chosen))


def download_all(bucket: str, keys: list[str], dest_dir: Path, workers: int) -> tuple[int, int]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    todo = [k for k in keys if not (dest_dir / Path(k).name).exists()]
    print(f"  {bucket}: {len(keys)} keys, {len(keys) - len(todo)} present, {len(todo)} to fetch")

    def _get(key: str) -> None:
        _s3().download_file(bucket, key, str(dest_dir / Path(key).name))

    n_ok = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_get, k): k for k in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            try:
                fut.result()
                n_ok += 1
            except Exception as exc:  # noqa: BLE001
                print(f"  FAIL {futs[fut]}: {exc!r}")
            if i % 100 == 0:
                print(f"  ... {i}/{len(todo)}")
    return n_ok, len(todo) - n_ok


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--which", default="goes,mrms", help="comma list: goes, mrms")
    p.add_argument("--hourly", action="store_true",
                   help="restrict both products to one file per hour")
    p.add_argument("--cadence-mrms", type=int, default=15,
                   help="minutes between MRMS files (ignored with --hourly)")
    p.add_argument("--goes-dir", default=str(REPO / "goes_data"))
    p.add_argument("--mrms-dir", default=str(REPO / "mrms_data"))
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--dry-run", action="store_true", help="list keys, download nothing")
    args = p.parse_args(argv)

    which = {w.strip() for w in args.which.split(",") if w.strip()}
    mrms_cadence = 60 if args.hourly else args.cadence_mrms
    s3 = _s3()

    goes_keys: list[str] = []
    mrms_keys: list[str] = []
    for case in CASES:
        start, end = _case_span(case)
        print(f"{case}: {start:%Y-%m-%d %H:%M} .. {end:%Y-%m-%d %H:%M} UTC")
        if "goes" in which:
            gk = list_goes_keys(s3, start, end)
            if args.hourly:
                gk = _thin_hourly(gk, _goes_scan_time)
            print(f"  GOES ch.{GOES_CHANNEL}: {len(gk)} files")
            goes_keys += gk
        if "mrms" in which:
            mk = list_mrms_keys(s3, start, end, mrms_cadence)
            print(f"  MRMS (every {mrms_cadence} min): {len(mk)} files")
            mrms_keys += mk

    goes_keys = sorted(set(goes_keys))
    mrms_keys = sorted(set(mrms_keys))
    print(f"\ntotal: {len(goes_keys)} GOES, {len(mrms_keys)} MRMS")
    if args.dry_run:
        return 0

    if goes_keys:
        ok, fail = download_all(GOES_BUCKET, goes_keys, Path(args.goes_dir), args.workers)
        print(f"GOES: {ok} downloaded, {fail} failed")
    if mrms_keys:
        ok, fail = download_all(MRMS_BUCKET, mrms_keys, Path(args.mrms_dir), args.workers)
        print(f"MRMS: {ok} downloaded, {fail} failed")
    return 0


def _thin_hourly(keys: list[str], scan_time) -> list[str]:
    """Keep the key closest to each whole hour."""
    by_hour: dict[dt.datetime, tuple[float, str]] = {}
    for k in keys:
        t = scan_time(k)
        hour = t.replace(minute=0, second=0, microsecond=0)
        d = abs((t - hour).total_seconds())
        if hour not in by_hour or d < by_hour[hour][0]:
            by_hour[hour] = (d, k)
    return sorted(k for _, k in by_hour.values())


if __name__ == "__main__":
    sys.exit(main())
