"""Gridding raw lightning observations (MRMS/NLDN, GOES GLM, NALMA) onto a
WRF domain's native grid, for comparison against WRF's Dynamic Lightning
Scheme flash counts (`wrf.read_dyn_lightning_flash_counts`).

Unlike the reflectivity/brightness-temperature comparisons in `plotting.py`
(which just reproject each source onto a shared map for a visual
side-by-side look), a lightning comparison needs actual spatial/temporal
*aggregation*: WRF's flash counts are per-grid-cell, per-output-interval
totals, so a fair comparison needs the observations binned into the same
grid cells over the same time interval -- "summing up values from MRMS
1-km boxes inside each WRF d1 grid box" and "accumulate the obs data within
the model output time interval and within each model grid box", per
`prompts/lightning_evaluation.md`. `bin_points_to_wrf_grid` is the shared
primitive behind that for every observation source (MRMS's already-gridded
~1km cells, GLM/NALMA's raw flash/event points).

See `prompts/lightning_threat_verification_notes.md` for background on
HRRR's separate (McCaul et al. 2009-based) `ltng` field, and the
conversation history for how WRF's `dyn_lightning_option` scheme itself
works.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import cartopy.crs as ccrs
import numpy as np

from . import goes as goes_reader
from . import mrms as mrms_reader
from . import nalma as nalma_reader
from . import wrf as wrf_reader


def bin_points_to_wrf_grid(
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    proj: ccrs.Projection,
    point_lon: np.ndarray,
    point_lat: np.ndarray,
    weight: np.ndarray | float = 1.0,
) -> np.ndarray:
    """Sum `weight` at each point into whichever WRF grid cell it falls in.

    `x_edges`/`y_edges` are the WRF grid's exact mass-cell edges in `proj`'s
    projected x/y coordinates (`wrf.read_cell_edges_xy` -- derived from the
    staggered U/V-point lat/lon, not extrapolated from the mass-point
    centers). `proj` is that domain's Lambert projection
    (`wrf.get_lambert_projection`). Both the grid edges and `point_lon`/
    `point_lat` get projected through this same `proj` before binning, so
    everything lives in one shared Lambert x/y plane -- this bins the
    projected points against the projected edges, equivalent to a 2D
    histogram but letting the caller supply a per-point weight (flash
    count, density value, ...) rather than counting points 1-for-1.

    `weight` broadcasts against `point_lon`/`point_lat`; pass a scalar
    (default 1.0) to just count points per cell.

    Points outside the WRF grid's extent, or with a non-finite lon/lat/
    weight, are silently dropped (not an error -- an observation source
    commonly covers a much larger area than the WRF domain, e.g. MRMS's
    CONUS-wide NLDN grid).

    Returns a `(len(y_edges)-1, len(x_edges)-1)` array (matching the WRF
    mass grid's shape), the summed weight per cell.
    """
    point_lon = np.asarray(point_lon)
    point_lat = np.asarray(point_lat)
    weight_arr = np.broadcast_to(np.asarray(weight, dtype=float), point_lon.shape)

    pt_xyz = proj.transform_points(ccrs.PlateCarree(), point_lon.ravel(), point_lat.ravel())
    px, py = pt_xyz[:, 0], pt_xyz[:, 1]
    w = weight_arr.ravel()

    col = np.searchsorted(x_edges, px) - 1
    row = np.searchsorted(y_edges, py) - 1

    ny, nx = len(y_edges) - 1, len(x_edges) - 1
    valid = (
        np.isfinite(px) & np.isfinite(py) & np.isfinite(w)
        & (col >= 0) & (col < nx) & (row >= 0) & (row < ny)
    )

    grid = np.zeros((ny, nx), dtype=float)
    np.add.at(grid, (row[valid], col[valid]), w[valid])
    return grid


def _crop_to_pad(
    lon: np.ndarray, lat: np.ndarray, values: np.ndarray, extent: tuple[float, float, float, float], pad_deg: float = 0.1
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Crop a 2D lon/lat/values grid to `extent` padded by `pad_deg`, to
    avoid projecting/binning a whole-CONUS observation grid when only a
    small area around the WRF domain is needed. Simpler than
    `plotting._crop_to_extent` (no pcolormesh-edge-artifact padding, since
    the result here is only ever binned, never drawn directly)."""
    lon_min, lon_max, lat_min, lat_max = extent
    mask = (
        (lon >= lon_min - pad_deg) & (lon <= lon_max + pad_deg)
        & (lat >= lat_min - pad_deg) & (lat <= lat_max + pad_deg)
    )
    if not mask.any():
        return lon, lat, values
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    sl = (slice(rows.min(), rows.max() + 1), slice(cols.min(), cols.max() + 1))
    return lon[sl], lat[sl], values[sl]


def mrms_cg_counts_on_wrf_grid(
    wrf_file: str | Path,
    time_start: dt.datetime,
    time_end: dt.datetime,
    mrms_dir: str | Path,
    auto_download: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dt.datetime]]:
    """Ground-based (NLDN) cloud-to-ground flash count per WRF grid cell,
    summed over `(time_start, time_end]` -- the observation-side counterpart
    of `wrf.read_dyn_lightning_flash_counts`'s `cg_pos + cg_neg`.

    Finds every `mrms.NLDN_CG_PRODUCT` (1-min) file in the window (S3
    listing + download if `auto_download`, else local-only -- see
    `mrms.find_mrms_files_in_range`/`download_mrms_file`), reads each,
    crops to the WRF domain's extent, and bins+sums every ~1km MRMS cell's
    value into whichever WRF grid cell it falls in
    (`bin_points_to_wrf_grid`), accumulating across all files in the
    window. See `mrms.read_nldn_cg_density` for the per-file value
    convention and its open caveat.

    Returns (lon, lat, counts, scan_times) on the WRF grid; `scan_times` is
    the sorted list of MRMS scan times actually summed (empty if none were
    found -- distinguish a genuine zero-flash window from a missing-data
    gap by checking whether this is empty).
    """
    wrf_lon, wrf_lat, _, _ = wrf_reader.read_field(wrf_file, "XLONG")
    proj = wrf_reader.get_lambert_projection(wrf_file)
    x_edges, y_edges = wrf_reader.read_cell_edges_xy(wrf_file, proj)
    extent = (
        float(np.nanmin(wrf_lon)), float(np.nanmax(wrf_lon)),
        float(np.nanmin(wrf_lat)), float(np.nanmax(wrf_lat)),
    )

    mrms_dir = Path(mrms_dir)
    if auto_download:
        matches = mrms_reader.find_mrms_files_in_range(
            time_start, time_end, product=mrms_reader.NLDN_CG_PRODUCT
        )
    else:
        # No S3 call at all here -- `find_mrms_files_in_range` always lists
        # S3 to discover what exists, even when every file is already
        # local, which hangs forever on an internet-less compute node.
        matches = mrms_reader.find_local_mrms_files_in_range(
            time_start, time_end, mrms_dir, product=mrms_reader.NLDN_CG_PRODUCT
        )
    mrms_dir.mkdir(parents=True, exist_ok=True)

    counts = np.zeros_like(wrf_lon, dtype=float)
    scan_times: list[dt.datetime] = []
    for bucket, key, scan_time in matches:
        dest = mrms_dir / Path(key).name
        if not dest.exists():
            if not auto_download:
                continue
            mrms_reader.download_mrms_file(bucket, key, mrms_dir)
        lon, lat, values, _ = mrms_reader.read_nldn_cg_density(dest)
        lon_c, lat_c, values_c = _crop_to_pad(lon, lat, values, extent)
        counts += bin_points_to_wrf_grid(x_edges, y_edges, proj, lon_c, lat_c, values_c)
        scan_times.append(scan_time)

    return wrf_lon, wrf_lat, counts, scan_times


def glm_total_counts_on_wrf_grid(
    wrf_file: str | Path,
    time_start: dt.datetime,
    time_end: dt.datetime,
    goes_dir: str | Path,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    auto_download: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dt.datetime]]:
    """Total-lightning (IC + CC + CG, undifferentiated) flash count per WRF
    grid cell, summed over `(time_start, time_end]` -- the observation-side
    counterpart of `wrf.read_dyn_lightning_flash_counts`'s
    `cg_pos + cg_neg + ic`.

    Finds every GLM-L2-LCFA file (one per ~20-second scan) in the window
    (S3 listing + download if `auto_download`, else local-only), reads
    each, and bins each flash centroid into whichever WRF grid cell it
    falls in (`bin_points_to_wrf_grid`, weight=1 per flash -- unlike MRMS's
    already-gridded density values, GLM flashes are discrete points),
    accumulating across all files in the window. No cropping step is
    needed first the way MRMS's whole-CONUS grid needs one: GLM files are
    already small (one scan's worth of flashes, full-disk but sparse), so
    `bin_points_to_wrf_grid`'s own out-of-domain dropping is cheap enough
    on its own.

    Returns (lon, lat, counts, scan_times) on the WRF grid; `scan_times` is
    the sorted list of GLM scan start times actually summed (empty if none
    were found).
    """
    wrf_lon, wrf_lat, _, _ = wrf_reader.read_field(wrf_file, "XLONG")
    proj = wrf_reader.get_lambert_projection(wrf_file)
    x_edges, y_edges = wrf_reader.read_cell_edges_xy(wrf_file, proj)

    goes_dir = Path(goes_dir)
    if auto_download:
        matches = goes_reader.find_glm_files_in_range(time_start, time_end, satellite=satellite)
    else:
        # No S3 call at all here -- same reasoning as
        # `mrms_cg_counts_on_wrf_grid`'s local-vs-S3 branch above.
        matches = goes_reader.find_local_glm_files_in_range(
            time_start, time_end, goes_dir, satellite=satellite
        )
    goes_dir.mkdir(parents=True, exist_ok=True)

    counts = np.zeros_like(wrf_lon, dtype=float)
    scan_times: list[dt.datetime] = []
    for bucket, key, scan_time in matches:
        dest = goes_dir / Path(key).name
        if not dest.exists():
            if not auto_download:
                continue
            goes_reader.download_goes_file(bucket, key, goes_dir)
        flash_lon, flash_lat, _energy, _start, _end = goes_reader.read_glm_flashes(dest)
        if flash_lon.size > 0:
            counts += bin_points_to_wrf_grid(x_edges, y_edges, proj, flash_lon, flash_lat)
        scan_times.append(scan_time)

    return wrf_lon, wrf_lat, counts, scan_times


def nalma_source_counts_on_wrf_grid(
    wrf_file: str | Path,
    time_start: dt.datetime,
    time_end: dt.datetime,
    nalma_dir: str | Path,
    token_file: str | Path = nalma_reader.DEFAULT_TOKEN_FILE,
    auto_download: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Raw VHF source count per WRF grid cell, summed over
    `(time_start, time_end]`.

    **Not directly comparable in magnitude** to WRF's or GLM's *flash*
    counts -- see the module docstring in `nalma.py`: a single flash
    typically produces many tens to hundreds of individual VHF sources, so
    this is a source-density field, not a flash-count field. Useful as a
    qualitative/spatial check (does the activity line up with WRF and
    GLM?), not for a like-for-like count comparison.

    Finds every NALMA granule overlapping the window (CMR search, no auth
    needed; downloads require `token_file` -- see `nalma.py`), reads each,
    filters sources to those actually falling in `(time_start, time_end]`
    (a granule's own window can extend past the requested one), and bins
    each source's location into whichever WRF grid cell it falls in
    (`bin_points_to_wrf_grid`, weight=1 per source). No chi^2/power quality
    filtering is applied yet (see `nalma.read_nalma_sources`) -- some
    fraction of raw sources are likely poorly-constrained/spurious
    solutions, though most such sources fall far outside NALMA's ~90-km
    physical network footprint and so outside the WRF domain entirely,
    naturally dropped by `bin_points_to_wrf_grid`.

    Returns (lon, lat, counts, granule_names) on the WRF grid;
    `granule_names` is the sorted list of NALMA filenames actually used
    (empty if none were found/downloaded).
    """
    wrf_lon, wrf_lat, _, _ = wrf_reader.read_field(wrf_file, "XLONG")
    proj = wrf_reader.get_lambert_projection(wrf_file)
    x_edges, y_edges = wrf_reader.read_cell_edges_xy(wrf_file, proj)

    matches = nalma_reader.find_nalma_files_in_range(time_start, time_end)
    nalma_dir = Path(nalma_dir)
    nalma_dir.mkdir(parents=True, exist_ok=True)

    counts = np.zeros_like(wrf_lon, dtype=float)
    granule_names: list[str] = []
    t_start64 = np.datetime64(time_start)
    t_end64 = np.datetime64(time_end)
    for url, name, _gstart, _gend in matches:
        dest = nalma_dir / name
        if not dest.exists():
            if not auto_download:
                continue
            nalma_reader.download_nalma_file(url, nalma_dir, token_file=token_file)
        lon, lat, times = nalma_reader.read_nalma_sources(dest)
        in_window = (times > t_start64) & (times <= t_end64)
        if np.any(in_window):
            counts += bin_points_to_wrf_grid(x_edges, y_edges, proj, lon[in_window], lat[in_window])
        granule_names.append(name)

    return wrf_lon, wrf_lat, counts, granule_names
