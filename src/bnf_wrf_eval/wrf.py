"""Reading fields and map projection info from LASSO-BNF wrfout files."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import cartopy.crs as ccrs
import netCDF4
import numpy as np


def get_lambert_projection(wrf_file: str | Path) -> ccrs.LambertConformal:
    """Cartopy Lambert Conformal projection matching a WRF domain's grid.

    Mirrors the convention wrf-python uses: the cone is defined by the
    domain's TRUELAT1/TRUELAT2 standard parallels, centered on STAND_LON /
    MOAD_CEN_LAT, on the spherical earth (radius 6,370,000 m) that WRF
    itself assumes rather than cartopy's default WGS84 ellipsoid.
    """
    with netCDF4.Dataset(wrf_file) as nc:
        globe = ccrs.Globe(ellipse=None, semimajor_axis=6370000, semiminor_axis=6370000)
        return ccrs.LambertConformal(
            central_longitude=float(nc.STAND_LON),
            central_latitude=float(nc.MOAD_CEN_LAT),
            standard_parallels=(float(nc.TRUELAT1), float(nc.TRUELAT2)),
            globe=globe,
        )


def read_field(
    wrf_file: str | Path, varname: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Read a single-time 2D field plus its lat/lon grid from a wrfout file.

    Returns (lon, lat, values, valid_time); lon/lat/values all share shape
    (south_north, west_east).
    """
    with netCDF4.Dataset(wrf_file) as nc:
        values = np.ma.filled(nc.variables[varname][0, ...], np.nan)
        lon = np.ma.filled(nc.variables["XLONG"][0, ...], np.nan)
        lat = np.ma.filled(nc.variables["XLAT"][0, ...], np.nan)
        time_str = b"".join(nc.variables["Times"][0, :]).decode("utf-8")
        valid_time = dt.datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")
    return lon, lat, values, valid_time


def _get_wrf_getvar():
    """Import wrf-python's `getvar` lazily, working around two issues from
    it being last released in 2021 (see pyproject.toml for the packaging
    issues just getting it installed):

    - It needs `numpy.float_`, removed in NumPy 2.0.
    - Its `is_multi_file()` helper decides whether an input is a single
      NetCDF file vs. a sequence of them via `isinstance(wrfin, Iterable)`.
      Modern `netCDF4.Dataset` now implements `__iter__` (dict-like, over
      variable names) -- something wrf-python never anticipated -- so a
      plain open Dataset now gets misdetected as a multi-file sequence,
      and wrf-python then tries to `copy()` it (to get a "resettable"
      iterator over the sequence), which `Dataset` doesn't support,
      crashing with "Dataset is not picklable". Patched to still treat a
      plain Dataset as a single file; a real multi-file sequence (a list
      of Datasets, e.g.) is unaffected since it isn't a Dataset itself.
    """
    if not hasattr(np, "float_"):
        np.float_ = np.float64

    from wrf import util as wrf_util

    if not isinstance(wrf_util.is_multi_file, _PatchedIsMultiFile):
        wrf_util.is_multi_file = _PatchedIsMultiFile(wrf_util.is_multi_file)

    from wrf import getvar

    return getvar


class _PatchedIsMultiFile:
    def __init__(self, orig):
        self._orig = orig

    def __call__(self, wrfin):
        if isinstance(wrfin, netCDF4.Dataset):
            return False
        return self._orig(wrfin)


def read_cloud_top_temperature(
    wrf_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Cloud-top temperature (K), via `wrf-python`'s `ctt` diagnostic.

    `ctt` derives a physically-based cloud-top temperature from the model
    state (an optical-depth-threshold method over the hydrometeor
    profiles), designed specifically as the model-side analog of what a
    clean IR window channel -- e.g. GOES ABI channel 13, 10.3 micron --
    observes, for comparison against satellite brightness temperature (see
    `goes.read_brightness_temperature`). Falls back to clear-sky
    surface/skin temperature wherever `ctt` judges there's no cloud, same
    as it always does.
    """
    getvar = _get_wrf_getvar()

    with netCDF4.Dataset(wrf_file) as nc:
        ctt = getvar(nc, "ctt", units="K")
        lon = np.asarray(ctt.coords["XLONG"].values)
        lat = np.asarray(ctt.coords["XLAT"].values)
        values = np.asarray(ctt.values)
        time_str = b"".join(nc.variables["Times"][0, :]).decode("utf-8")
        valid_time = dt.datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")
    return lon, lat, values, valid_time


def read_column_max_reflectivity(
    wrf_file: str | Path, varname: str = "REFL_10CM"
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Column-max radar reflectivity (dBZ), computed from the 3D field.

    wrfout also carries a precomputed `REFD_MAX` diagnostic, but at least
    in the 20250917lassobnfwrfhrrr3 test run it is identically zero at
    every time (likely never accumulated/reset correctly for this run),
    so this takes the max over height of the instantaneous 3D
    `REFL_10CM` field instead, which does hold real values.
    """
    with netCDF4.Dataset(wrf_file) as nc:
        values = np.nanmax(np.ma.filled(nc.variables[varname][0, ...], np.nan), axis=0)
        lon = np.ma.filled(nc.variables["XLONG"][0, ...], np.nan)
        lat = np.ma.filled(nc.variables["XLAT"][0, ...], np.nan)
        time_str = b"".join(nc.variables["Times"][0, :]).decode("utf-8")
        valid_time = dt.datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")
    return lon, lat, values, valid_time
