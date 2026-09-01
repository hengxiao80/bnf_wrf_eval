"""Reading fields from HRRR native-level analysis GRIB2 files."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np

from ._eccodes_setup import ensure_eccodes_loadable
from .wrf import brightness_temperature_from_olr


def _get_eccodes():
    """Import eccodes lazily, on first actual use rather than at module
    import time.

    `ensure_eccodes_loadable()` does real filesystem work -- dlopen'ing
    several shared libraries out of the 85 MB `ecmwflibs` package, plus
    writing a symlink under /tmp -- which used to run unconditionally the
    moment `bnf_wrf_eval.plotting`/`bnf_wrf_eval.hrrr` got imported, even
    in a cell that never ends up reading a GRIB file. Deferring it here
    means a plain `import` cell stays cheap; the cost (still cached after
    the first call, via `ensure_eccodes_loadable`'s `functools.cache`)
    only shows up once something actually calls `read_field`.
    """
    ensure_eccodes_loadable()
    import eccodes

    return eccodes


def read_field(
    grib_file: str | Path,
    *,
    type_of_level: str,
    short_name: str | None = None,
    parameter_category: int | None = None,
    parameter_number: int | None = None,
    level: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Read one field from a HRRR GRIB2 file by matching GRIB2 metadata.

    Some HRRR fields (e.g. top-of-atmosphere upward long-wave radiation,
    i.e. OLR) have no `shortName` in the local ecCodes parameter tables and
    decode as "unknown"; for those, select by
    (parameter_category, parameter_number) instead of short_name.

    Returns (lon, lat, values, valid_time); lon is wrapped to [-180, 180].
    """
    eccodes = _get_eccodes()
    with open(grib_file, "rb") as fh:
        while True:
            gid = eccodes.codes_grib_new_from_file(fh)
            if gid is None:
                raise ValueError(
                    f"No message matching type_of_level={type_of_level!r}, "
                    f"short_name={short_name!r}, "
                    f"parameter=({parameter_category}, {parameter_number}) "
                    f"found in {grib_file}"
                )
            try:
                if eccodes.codes_get(gid, "typeOfLevel") != type_of_level:
                    continue
                if level is not None and eccodes.codes_get(gid, "level") != level:
                    continue
                if short_name is not None and eccodes.codes_get(gid, "shortName") != short_name:
                    continue
                if (
                    parameter_category is not None
                    and eccodes.codes_get(gid, "parameterCategory") != parameter_category
                ):
                    continue
                if (
                    parameter_number is not None
                    and eccodes.codes_get(gid, "parameterNumber") != parameter_number
                ):
                    continue

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
                return lon, lat, values, valid_time
            finally:
                eccodes.codes_release(gid)


def read_composite_reflectivity(
    grib_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Composite (column-max) radar reflectivity, dBZ."""
    return read_field(grib_file, type_of_level="atmosphere", short_name="refc")


def read_olr(grib_file: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Outgoing long-wave radiation (upward LW flux at top of atmosphere), W/m^2."""
    return read_field(
        grib_file, type_of_level="nominalTop", parameter_category=5, parameter_number=4
    )


def read_brightness_temperature_from_olr(
    grib_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """IR-window brightness temperature (K), converted from HRRR's TOA OLR.

    Reads HRRR's top-of-atmosphere upward long-wave flux (`read_olr`) and
    passes it through `wrf.brightness_temperature_from_olr` -- see that
    function for the formula and caveats. This is the cheap,
    same-recipe-as-WRF counterpart to `read_simulated_brightness_temperature`
    (HRRR's own CRTM `SBT114`), so a WRF-vs-HRRR OLR-derived Tb comparison
    stays apples-to-apples.
    """
    lon, lat, olr, valid_time = read_olr(grib_file)
    return lon, lat, brightness_temperature_from_olr(olr), valid_time


def read_simulated_brightness_temperature(
    grib_file: str | Path, short_name: str = "SBT114"
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Simulated IR-window-channel brightness temperature (K), from HRRR's
    CRTM-derived simulated-satellite fields -- the observation-comparable
    analog of WRF's `ctt` (see `wrf.read_cloud_top_temperature`), letting a
    WRF-vs-HRRR-vs-GOES three-way brightness-temperature comparison use the
    same physical quantity throughout.

    HRRR/UPP's output carries several `SBT1{sat}{channel}`-named fields
    (`sat` 1/2 = a legacy "GOES 11"/"GOES 12" designation baked into the
    GRIB shortName regardless of which real satellite is actually being
    emulated; `channel` 3 ~= water vapor ~6.7 micron, `channel` 4 ~= the
    ~10.7 micron IR window). `SBT114` (the default) is the IR window
    channel, comparable to GOES-R ABI channel 13 (10.3 micron) and to
    WRF's `ctt`. A near-duplicate `SBT124` field also exists in the files;
    which of the two (if they differ meaningfully) is the better match
    hasn't been checked against real data yet -- pass `short_name="SBT124"`
    to compare.

    Checked against real data: the field uses 9999 K as a fill/missing
    value at some pixels (~5% of the grid in a CONUS test case) -- this
    masks that out as NaN rather than letting it dominate the color scale.
    """
    lon, lat, values, valid_time = read_field(
        grib_file, type_of_level="nominalTop", short_name=short_name
    )
    values = np.where(values == 9999.0, np.nan, values)
    return lon, lat, values, valid_time
