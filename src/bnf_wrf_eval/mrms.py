"""Finding, downloading, and reading MRMS (Multi-Radar/Multi-Sensor)
composite reflectivity, for comparison against WRF/HRRR radar reflectivity.

Unlike the GOES-vs-OLR situation, MRMS composite reflectivity is a direct,
apples-to-apples match for WRF's simulated reflectivity: same physical
quantity (dBZ), same basis as HRRR's `refc` -- see `hrrr.py`. NOAA's own
verification practice (METplus, RRFS-vs-HRRR literature) typically scores
this via threshold/skill metrics rather than raw grid differencing, because
of storm-displacement double-penalty effects, but a direct side-by-side dBZ
panel (what `plotting.plot_refl_comparison` does) is still a reasonable
first look.

Also reads ground-based cloud-to-ground lightning flash density from the
National Lightning Detection Network (NLDN), gridded by MRMS
(`read_nldn_cg_density`, `NLDN_CG_PRODUCT`) -- for comparing against WRF's
dynamic-lightning-scheme CG flash counts (`wrf.read_dyn_lightning_flash_counts`)
and/or `lightning.py`'s grid-binning helpers. Unlike composite reflectivity,
this is a genuinely sparse, mostly-zero, near-instantaneous-count product
rather than a smooth continuous field -- see `read_nldn_cg_density`'s
docstring for its value convention and open caveat.
"""

from __future__ import annotations

import datetime as dt
import gzip
import tempfile
from pathlib import Path

import numpy as np

# Public, unauthenticated AWS Open Data bucket (NOAA's NODD program) --
# same access pattern as the GOES buckets in goes.py, no AWS credentials
# needed. Archived here since October 2020, so all of this project's case
# days are covered; NCEP's own real-time mirror only keeps ~1-2 days.
DEFAULT_BUCKET = "noaa-mrms-pds"

# The composite (column-max), quality-controlled reflectivity mosaic --
# the direct analog of WRF's column-max REFL_10CM and HRRR's `refc`, not
# the 3D MergedReflectivityQC product (33 vertical levels). The "_00.50"
# suffix is MRMS's standard 2D-product naming convention, not a physical
# height.
PRODUCT = "MergedReflectivityQCComposite_00.50"

# Ground-based cloud-to-ground lightning flash density, from Vaisala's
# National Lightning Detection Network (NLDN), gridded onto MRMS's native
# ~0.01-deg (~1 km) CONUS grid. Confirmed present on the public bucket at
# four averaging windows (001/005/015/030 min); 1-min is used here (rather
# than a coarser window) specifically to minimize double-counting when
# summing several samples across a comparison window -- each 1-min product
# is a real time-slice, not an overlapping moving average of a longer one.
# See `read_nldn_cg_density` for the value convention.
NLDN_CG_PRODUCT = "NLDN_CG_001min_AvgDensity_00.00"


def find_mrms_file(
    time: dt.datetime,
    bucket: str = DEFAULT_BUCKET,
    tolerance_minutes: float = 2.0,
    product: str = PRODUCT,
) -> tuple[str, str, dt.datetime]:
    """Find the MRMS composite reflectivity object whose scan time is
    closest to `time`, within `tolerance_minutes` (MRMS's ~2-minute
    cadence, so keep this tight).

    Lists the relevant S3 prefix(es) (the bucket is public; no AWS
    credentials needed) and picks the closest-time match. Doesn't download
    anything.

    Returns (bucket, key, scan_time) for the matching object.

    Raises FileNotFoundError if nothing is found within tolerance.
    """
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config

    s3 = boto3.client(
        "s3",
        config=Config(
            signature_version=UNSIGNED,
            connect_timeout=10,
            read_timeout=30,
            retries={"max_attempts": 2},
        ),
    )

    # Objects are keyed by <day>/<file>, so normally only the target day
    # needs listing; only also check the neighboring day if `time` is
    # close enough to midnight that the closest scan could plausibly fall
    # on the other side of it.
    tolerance = dt.timedelta(minutes=tolerance_minutes)
    days_to_check = {time.date()}
    if time - dt.datetime.combine(time.date(), dt.time.min) < tolerance:
        days_to_check.add((time - dt.timedelta(days=1)).date())
    if dt.datetime.combine(time.date(), dt.time.max) - time < tolerance:
        days_to_check.add((time + dt.timedelta(days=1)).date())

    candidates = []
    for day in days_to_check:
        prefix = f"CONUS/{product}/{day:%Y%m%d}/MRMS_{product}_{day:%Y%m%d}-"
        candidates.extend(_list_mrms_objects(s3, bucket, prefix))

    if not candidates:
        raise FileNotFoundError(f"No {product} files found near {time} in s3://{bucket}")

    def _scan_time(obj) -> dt.datetime:
        # Key looks like .../MRMS_..._YYYYMMDD-HHMMSS.grib2.gz
        stem = Path(obj["Key"]).name.removesuffix(".grib2.gz")
        return dt.datetime.strptime(stem.rsplit("_", 1)[1], "%Y%m%d-%H%M%S")

    best = min(candidates, key=lambda obj: abs((_scan_time(obj) - time).total_seconds()))
    best_time = _scan_time(best)
    if abs((best_time - time).total_seconds()) > tolerance.total_seconds():
        raise FileNotFoundError(
            f"Closest {product} scan to {time} is {best_time} "
            f"({abs((best_time - time).total_seconds()) / 60:.1f} min away, "
            f"outside the {tolerance_minutes}-min tolerance)"
        )

    return bucket, best["Key"], best_time


def _list_mrms_objects(s3, bucket: str, prefix: str, start_after: str | None = None) -> list[dict]:
    """`list_objects_v2` under `prefix`, paginated -- a single call
    silently truncates at 1000 keys, which the ~2-min-cadence composite
    reflectivity product (~720/day) never hits but the ~1-min NLDN product
    (~1440/day) does, so an unpaginated listing can silently miss
    everything past whatever time ~1000 objects reaches (confirmed for
    real: a day's NLDN listing truncated at 16:52 UTC). `start_after`, if
    given, is passed through to skip listing keys lexicographically at or
    before it (object keys sort by their `HHMMSS` timestamp suffix, so this
    is also a real optimization when only a narrow time window is needed).
    """
    objects = []
    kwargs = {"Bucket": bucket, "Prefix": prefix}
    if start_after is not None:
        kwargs["StartAfter"] = start_after
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(**kwargs):
        objects.extend(page.get("Contents", []))
    return objects


def find_mrms_files_in_range(
    time_start: dt.datetime,
    time_end: dt.datetime,
    bucket: str = DEFAULT_BUCKET,
    product: str = NLDN_CG_PRODUCT,
) -> list[tuple[str, str, dt.datetime]]:
    """List every MRMS object of `product` scanned in `(time_start,
    time_end]`, sorted by scan time -- for accumulating a near-continuous
    product (like `NLDN_CG_001min_AvgDensity`) across a comparison window,
    as opposed to `find_mrms_file`'s single closest-match lookup.

    Returns a list of (bucket, key, scan_time) tuples; empty if nothing
    falls in the window (not an error, since a quiet period may genuinely
    have no lightning to report -- though the product itself still writes
    files every ~1 min regardless of activity, so an empty result more
    likely means the window fell entirely outside the archive's coverage).
    """
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config

    s3 = boto3.client(
        "s3",
        config=Config(
            signature_version=UNSIGNED,
            connect_timeout=10,
            read_timeout=30,
            retries={"max_attempts": 2},
        ),
    )

    days_to_check = {time_start.date(), time_end.date()}
    candidates = []
    for day in sorted(days_to_check):
        prefix = f"CONUS/{product}/{day:%Y%m%d}/MRMS_{product}_{day:%Y%m%d}-"
        # Skip listing keys before the window even starts on this day --
        # both a real speedup and (with `_list_mrms_objects`'s pagination)
        # what actually keeps a late-day window from being missed.
        window_start_on_day = max(time_start, dt.datetime.combine(day, dt.time.min))
        start_after = f"{prefix}{window_start_on_day:%H%M%S}"
        candidates.extend(_list_mrms_objects(s3, bucket, prefix, start_after=start_after))

    def _scan_time(obj) -> dt.datetime:
        stem = Path(obj["Key"]).name.removesuffix(".grib2.gz")
        return dt.datetime.strptime(stem.rsplit("_", 1)[1], "%Y%m%d-%H%M%S")

    matches = []
    for obj in candidates:
        scan_time = _scan_time(obj)
        if time_start < scan_time <= time_end:
            matches.append((bucket, obj["Key"], scan_time))
    matches.sort(key=lambda item: item[2])
    return matches


def find_local_mrms_files_in_range(
    time_start: dt.datetime,
    time_end: dt.datetime,
    mrms_dir: str | Path,
    product: str = NLDN_CG_PRODUCT,
) -> list[tuple[str, str, dt.datetime]]:
    """Local counterpart of `find_mrms_files_in_range`: every already-
    downloaded `product` file in `mrms_dir` scanned in `(time_start,
    time_end]`, resolved purely from filenames -- no S3 call.

    Use this on hosts without outbound internet (compute nodes): unlike
    `find_local_mrms_file`'s single closest-match lookup,
    `find_mrms_files_in_range` (needed to accumulate a whole window) always
    lists S3 to discover what exists, even when every file is already
    local -- exactly the historical bug this project has hit before with
    other observation sources. Pre-populate `mrms_dir` with
    `scripts/download_obs.py` from a host that has internet.

    Returns a list of `("", path_as_str, scan_time)` tuples (empty bucket,
    since there's nothing to download) sorted by scan time, in the same
    shape `find_mrms_files_in_range` returns -- so callers can treat the
    two interchangeably.
    """
    mrms_dir = Path(mrms_dir)
    matches = []
    for path in mrms_dir.glob(f"MRMS_{product}_*.grib2.gz"):
        try:
            scan_time = _mrms_scan_time_from_name(path.name)
        except (IndexError, ValueError):
            continue
        if time_start < scan_time <= time_end:
            matches.append(("", str(path), scan_time))
    matches.sort(key=lambda item: item[2])
    return matches


def _mrms_scan_time_from_name(name: str) -> dt.datetime:
    """Scan time parsed from a `MRMS_..._YYYYMMDD-HHMMSS.grib2.gz` filename."""
    stem = name.removesuffix(".grib2.gz")
    return dt.datetime.strptime(stem.rsplit("_", 1)[1], "%Y%m%d-%H%M%S")


def find_local_mrms_file(
    time: dt.datetime,
    mrms_dir: str | Path,
    tolerance_minutes: float = 10.0,
    product: str = PRODUCT,
) -> Path:
    """The already-downloaded MRMS file in `mrms_dir` whose scan time is
    closest to `time`, within `tolerance_minutes` -- resolved purely from
    local filenames, no S3 call.

    Use this on hosts without outbound internet (compute nodes): the
    S3-backed `find_mrms_file` blocks forever on `list_objects_v2` there.
    Pre-populate `mrms_dir` with `scripts/download_obs.py` from a host that
    has internet.

    Raises FileNotFoundError if nothing matches within tolerance.
    """
    mrms_dir = Path(mrms_dir)
    best: Path | None = None
    best_delta = None
    for path in mrms_dir.glob(f"MRMS_{product}_*.grib2.gz"):
        try:
            scan = _mrms_scan_time_from_name(path.name)
        except (IndexError, ValueError):
            continue
        delta = abs((scan - time).total_seconds())
        if best_delta is None or delta < best_delta:
            best, best_delta = path, delta
    if best is None or best_delta > tolerance_minutes * 60:
        raise FileNotFoundError(
            f"No local {product} file in {mrms_dir} within {tolerance_minutes} "
            f"min of {time}"
            + (f" (closest is {best_delta / 60:.1f} min away)" if best is not None else "")
        )
    return best


def download_mrms_file(bucket: str, key: str, dest_dir: str | Path) -> Path:
    """Download one S3 object (as found by `find_mrms_file`) into
    `dest_dir`, keeping just its filename (still gzip-compressed -- see
    `read_composite_reflectivity`). Skips the download if a file of that
    name is already there."""
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / Path(key).name
    if not dest.exists():
        s3 = boto3.client(
        "s3",
        config=Config(
            signature_version=UNSIGNED,
            connect_timeout=10,
            read_timeout=30,
            retries={"max_attempts": 2},
        ),
    )
        s3.download_file(bucket, key, str(dest))
    return dest


def _get_eccodes():
    """Import eccodes lazily -- see hrrr._get_eccodes for why."""
    from ._eccodes_setup import ensure_eccodes_loadable

    ensure_eccodes_loadable()
    import eccodes

    return eccodes


def read_composite_reflectivity(
    mrms_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Read composite reflectivity (dBZ) and lat/lon from a downloaded MRMS
    GRIB2 file.

    Unlike HRRR's GRIB2 files (many fields per file, selected by GRIB2
    metadata) and unlike GOES's ABI fixed grid, an MRMS composite
    reflectivity file holds exactly one field on a plain regular lat/lon
    grid -- no field selection or curvilinear/geostationary reprojection
    needed, just the same eccodes `latitudes`/`longitudes`/`values` array
    read `hrrr.read_field` already does.

    The downloaded file is gzip-compressed (`.grib2.gz`); this decompresses
    it to a temporary file first since eccodes needs a real file handle.

    MRMS uses a large negative sentinel (-999, confirmed against a real
    file) for "no radar coverage" at a pixel, distinct from -99, its
    separate placeholder for "a radar covers this pixel and sees
    essentially no echo" (most of any given CONUS grid, at most times).
    Only the former is masked to NaN here as genuinely missing data;
    -99 is left as-is; since it's well below `plotting.DEFAULT_REFL_LEVELS`
    (which starts at 5 dBZ, with `extend="max"` only), it renders as blank
    -- the normal radar-display convention -- the same way WRF's own
    -35 dBZ no-echo floor already does, without needing special-casing.

    Returns (lon, lat, values, valid_time); lon is wrapped to [-180, 180].
    """
    lon, lat, values, valid_time = _read_mrms_gz_grib2(mrms_file)
    values = np.where(values < -900, np.nan, values)
    return lon, lat, values, valid_time


def _read_mrms_gz_grib2(
    mrms_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Shared low-level reader behind `read_composite_reflectivity` and
    `read_nldn_cg_density`: decompresses a gzip'd single-message GRIB2 file
    and returns its raw lat/lon/values/valid_time, with no sentinel
    handling (each caller applies its own product-specific convention).
    """
    eccodes = _get_eccodes()

    with gzip.open(mrms_file, "rb") as gz:
        data = gz.read()

    with tempfile.NamedTemporaryFile(suffix=".grib2") as tmp:
        tmp.write(data)
        tmp.flush()
        with open(tmp.name, "rb") as fh:
            gid = eccodes.codes_grib_new_from_file(fh)
            if gid is None:
                raise ValueError(f"No GRIB message found in {mrms_file}")
            try:
                ni = eccodes.codes_get(gid, "Ni")
                nj = eccodes.codes_get(gid, "Nj")
                lat = np.asarray(eccodes.codes_get_array(gid, "latitudes")).reshape(nj, ni)
                lon = np.asarray(eccodes.codes_get_array(gid, "longitudes")).reshape(nj, ni)
                lon = np.where(lon > 180, lon - 360, lon)
                values = np.asarray(eccodes.codes_get_array(gid, "values")).reshape(nj, ni)
                valid_time = dt.datetime.strptime(
                    f"{eccodes.codes_get(gid, 'validityDate')}"
                    f"{eccodes.codes_get(gid, 'validityTime'):04d}",
                    "%Y%m%d%H%M",
                )
            finally:
                eccodes.codes_release(gid)

    return lon, lat, values, valid_time


def read_nldn_cg_density(
    mrms_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Read NLDN cloud-to-ground lightning density (see `NLDN_CG_PRODUCT`)
    and lat/lon from a downloaded MRMS GRIB2 file.

    Per NOAA/NWS Warning Decision Training Division documentation, each
    value is a per-~1km-grid-cell average flash count over the product's
    averaging window (1 min here) -- i.e. already normalized to "flashes in
    this cell", not a true areal density needing a separate area
    multiplication. Summing several consecutive 1-min files' values at a
    cell therefore gives an approximate total flash count at that cell over
    the summed period; this hasn't yet been cross-validated against a
    known-count case, so treat it as a working assumption pending a sanity
    check (see `prompts/lightning_evaluation.md`).

    Confirmed against a real downloaded file: **-1 is this product's
    "zero flashes this window" sentinel** (99.997% of the CONUS grid at any
    given time), not a missing-radar-coverage flag like composite
    reflectivity's -999 -- NLDN's ground-sensor coverage is effectively
    complete over CONUS, so there's no equivalent "no coverage" case to
    encode. Mapped to 0 here, not NaN.

    Returns (lon, lat, values, valid_time); lon is wrapped to [-180, 180].
    """
    lon, lat, values, valid_time = _read_mrms_gz_grib2(mrms_file)
    values = np.where(values <= -0.5, 0.0, values)
    return lon, lat, values, valid_time
