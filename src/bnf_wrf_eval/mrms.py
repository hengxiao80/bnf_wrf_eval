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

NOTE: nothing in this module has been run yet -- see the plan this was
written from. The S3 bucket/key layout was confirmed by actually listing
the bucket; the missing-value sentinel handling in
`read_composite_reflectivity` follows the common MRMS convention but hasn't
been checked against a real file here.
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


def find_mrms_file(
    time: dt.datetime,
    bucket: str = DEFAULT_BUCKET,
    tolerance_minutes: float = 2.0,
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

    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))

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
        prefix = f"CONUS/{PRODUCT}/{day:%Y%m%d}/MRMS_{PRODUCT}_{day:%Y%m%d}-"
        resp = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
        candidates.extend(resp.get("Contents", []))

    if not candidates:
        raise FileNotFoundError(f"No {PRODUCT} files found near {time} in s3://{bucket}")

    def _scan_time(obj) -> dt.datetime:
        # Key looks like .../MRMS_..._YYYYMMDD-HHMMSS.grib2.gz
        stem = Path(obj["Key"]).name.removesuffix(".grib2.gz")
        return dt.datetime.strptime(stem.rsplit("_", 1)[1], "%Y%m%d-%H%M%S")

    best = min(candidates, key=lambda obj: abs((_scan_time(obj) - time).total_seconds()))
    best_time = _scan_time(best)
    if abs((best_time - time).total_seconds()) > tolerance.total_seconds():
        raise FileNotFoundError(
            f"Closest {PRODUCT} scan to {time} is {best_time} "
            f"({abs((best_time - time).total_seconds()) / 60:.1f} min away, "
            f"outside the {tolerance_minutes}-min tolerance)"
        )

    return bucket, best["Key"], best_time


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
        s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
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
                values = np.where(values < -900, np.nan, values)
                valid_time = dt.datetime.strptime(
                    f"{eccodes.codes_get(gid, 'validityDate')}"
                    f"{eccodes.codes_get(gid, 'validityTime'):04d}",
                    "%Y%m%d%H%M",
                )
            finally:
                eccodes.codes_release(gid)

    return lon, lat, values, valid_time
