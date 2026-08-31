"""Finding, downloading, and reading GOES-R series ABI Cloud and Moisture
Imagery Product (CMIP) brightness temperature, for comparison against WRF
cloud-top temperature (see `wrf.read_cloud_top_temperature`).

NOTE: nothing in this module has been run yet -- see the plan this was
written from. In particular the ABI fixed-grid -> lat/lon reprojection in
`read_brightness_temperature` follows the standard, widely-documented GOES-R
convention but hasn't been exercised against a real file here.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np

# GOES-16 was the operational "GOES-East" satellite through April 2025, when
# NOAA switched that role to GOES-19 (after a GOES-16 cooling-system issue
# degraded its imager). Every WRF run date in this project so far
# (2025-05 onward) postdates that switch, so GOES-19 is the right default;
# pass satellite="G16" explicitly for older cases.
DEFAULT_SATELLITE = "G19"

# Public, unauthenticated AWS Open Data buckets -- no AWS credentials
# needed, just anonymous S3 access.
S3_BUCKET_BY_SATELLITE = {
    "G16": "noaa-goes16",
    "G17": "noaa-goes17",
    "G18": "noaa-goes18",
    "G19": "noaa-goes19",
}

# CONUS-sector Cloud and Moisture Imagery Product: 5-minute cadence,
# ~2 km resolution at nadir for the IR bands, and already converted to
# brightness temperature (no radiance -> Tb conversion needed, unlike the
# L1b radiance product).
PRODUCT = "ABI-L2-CMIPC"

# GOES-R CMIP files' "t" (scan midpoint) variable is seconds since this
# fixed epoch (the files don't vary this).
_GOES_EPOCH = dt.datetime(2000, 1, 1, 12, 0, 0)


def find_goes_file(
    time: dt.datetime,
    channel: int = 13,
    satellite: str = DEFAULT_SATELLITE,
    tolerance_minutes: float = 5.0,
) -> tuple[str, str, dt.datetime]:
    """Find the CONUS-sector GOES ABI CMIP object (for `channel`) whose scan
    start time is closest to `time`, within `tolerance_minutes`.

    Lists the relevant S3 prefixes (the bucket is public; no AWS
    credentials needed) and picks the closest-time match. Doesn't download
    anything.

    Channel 13 (10.3 micron, "Clean IR Longwave Window") is the standard
    clean-window channel most comparable to `wrf-python`'s `ctt`
    diagnostic; channel 11 (8.4 micron) is more useful for cloud-phase/dust
    discrimination than as a standalone cloud-top proxy, but is supported
    here too.

    Returns (bucket, key, scan_start_time) for the matching object.

    Raises FileNotFoundError if nothing is found within tolerance.
    """
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config

    bucket = S3_BUCKET_BY_SATELLITE[satellite]
    s3 = boto3.client(
        "s3",
        config=Config(
            signature_version=UNSIGNED,
            connect_timeout=10,
            read_timeout=30,
            retries={"max_attempts": 2},
        ),
    )

    # Objects are keyed by <product>/<year>/<day-of-year>/<hour>/<file>, so
    # list just the target hour (and the ones either side, in case the
    # closest scan crosses an hour boundary) rather than the whole day.
    candidates = []
    for hour_offset in (-1, 0, 1):
        probe_time = time + dt.timedelta(hours=hour_offset)
        prefix = (
            f"{PRODUCT}/{probe_time:%Y}/{probe_time:%j}/{probe_time:%H}/"
            f"OR_{PRODUCT}-M6C{channel:02d}_{satellite}_"
        )
        resp = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
        candidates.extend(resp.get("Contents", []))

    if not candidates:
        raise FileNotFoundError(
            f"No {PRODUCT} channel {channel} files found near {time} in s3://{bucket}"
        )

    def _scan_start(obj) -> dt.datetime:
        # Key looks like ..._sYYYYJJJHHMMSSs_eYYYYJJJHHMMSSs_cYYYYJJJHHMMSSs.nc
        s_field = obj["Key"].split("_s")[1].split("_")[0]  # "YYYYJJJHHMMSSs"
        return dt.datetime.strptime(s_field[:-1], "%Y%j%H%M%S")

    best = min(candidates, key=lambda obj: abs((_scan_start(obj) - time).total_seconds()))
    best_time = _scan_start(best)
    if abs((best_time - time).total_seconds()) > tolerance_minutes * 60:
        raise FileNotFoundError(
            f"Closest {PRODUCT} channel {channel} scan to {time} is {best_time} "
            f"({abs((best_time - time).total_seconds()) / 60:.1f} min away, "
            f"outside the {tolerance_minutes}-min tolerance)"
        )

    return bucket, best["Key"], best_time


def _goes_scan_start_from_name(name: str) -> dt.datetime:
    """Scan start time parsed from a CMIP filename's `_sYYYYJJJHHMMSSs` field."""
    s_field = name.split("_s")[1].split("_")[0]
    return dt.datetime.strptime(s_field[:-1], "%Y%j%H%M%S")


def find_local_goes_file(
    time: dt.datetime,
    goes_dir: str | Path,
    channel: int = 13,
    satellite: str = DEFAULT_SATELLITE,
    tolerance_minutes: float = 10.0,
) -> Path:
    """The already-downloaded GOES CMIP file in `goes_dir` whose scan start
    is closest to `time`, within `tolerance_minutes` -- resolved purely from
    local filenames, with no S3 call.

    Use this on hosts without outbound internet (e.g. compute nodes): the
    S3-backed `find_goes_file` blocks forever on `list_objects_v2` there.
    Pre-populate `goes_dir` with `scripts/download_obs.py` from a host that
    does have internet.

    Raises FileNotFoundError if `goes_dir` holds no matching file within
    tolerance.
    """
    goes_dir = Path(goes_dir)
    pattern = f"OR_{PRODUCT}-M6C{channel:02d}_{satellite}_s*.nc"
    best: Path | None = None
    best_delta = None
    for path in goes_dir.glob(pattern):
        try:
            scan = _goes_scan_start_from_name(path.name)
        except (IndexError, ValueError):
            continue
        delta = abs((scan - time).total_seconds())
        if best_delta is None or delta < best_delta:
            best, best_delta = path, delta
    if best is None or best_delta > tolerance_minutes * 60:
        raise FileNotFoundError(
            f"No local {PRODUCT} channel {channel} file in {goes_dir} within "
            f"{tolerance_minutes} min of {time}"
            + (f" (closest is {best_delta / 60:.1f} min away)" if best is not None else "")
        )
    return best


def download_goes_file(bucket: str, key: str, dest_dir: str | Path) -> Path:
    """Download one S3 object (as found by `find_goes_file`) into
    `dest_dir`, keeping just its filename. Skips the download if a file of
    that name is already there."""
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


def read_brightness_temperature(
    goes_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Read brightness temperature (K) and lat/lon from a GOES ABI CMIP
    NetCDF file.

    The file stores the image on the ABI fixed grid (N/S and E/W scan
    angles, in radians -- not lat/lon) plus a `goes_imager_projection`
    variable carrying the geostationary projection's parameters. This
    reconstructs lat/lon from that using cartopy's `Geostationary`
    projection (the same reprojection machinery the WRF/HRRR plots already
    use), following the standard GOES-R fixed-grid convention: projection
    x/y (meters) = scan angle (radians) * satellite height.

    Returns (lon, lat, values, scan_start_time); lon/lat/values share shape
    (y, x). Off-Earth pixels (space, beyond the visible disk) are NaN in
    all three.
    """
    import cartopy.crs as ccrs
    import netCDF4

    with netCDF4.Dataset(goes_file) as nc:
        tb = np.ma.filled(nc.variables["CMI"][:], np.nan)
        x = nc.variables["x"][:]
        y = nc.variables["y"][:]

        proj_info = nc.variables["goes_imager_projection"]
        sat_height = proj_info.perspective_point_height
        sat_lon = proj_info.longitude_of_projection_origin
        sweep_axis = proj_info.sweep_angle_axis
        semi_major = proj_info.semi_major_axis
        semi_minor = proj_info.semi_minor_axis

        X, Y = np.meshgrid(x * sat_height, y * sat_height)
        geos = ccrs.Geostationary(
            central_longitude=sat_lon,
            satellite_height=sat_height,
            sweep_axis=sweep_axis,
            globe=ccrs.Globe(semimajor_axis=semi_major, semiminor_axis=semi_minor),
        )
        lonlat = ccrs.PlateCarree().transform_points(geos, X, Y)
        lon, lat = lonlat[..., 0], lonlat[..., 1]

        # Off-Earth pixels (beyond the visible disk) transform to +/-inf;
        # mask them the same way the brightness temperatures already are.
        off_earth = ~np.isfinite(lon) | ~np.isfinite(lat)
        lon = np.where(off_earth, np.nan, lon)
        lat = np.where(off_earth, np.nan, lat)
        tb = np.where(off_earth, np.nan, tb)

        scan_start = _GOES_EPOCH + dt.timedelta(seconds=float(nc.variables["t"][:]))

    return lon, lat, tb, scan_start
