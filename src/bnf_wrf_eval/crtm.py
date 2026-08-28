"""Compute simulated GOES-East ABI channel 13 brightness temperature
directly from WRF model output via CRTM (through JCSDA's `pyCRTM`) -- the
same forward-model *approach* NCEP's UPP uses to produce HRRR's `SBT114`
field (see `hrrr.read_simulated_brightness_temperature`), applied directly
to our own WRF state instead of running the full UPP pipeline on it.

Uses GOES-16's ABI coefficients (`abi_g16`) rather than GOES-19's, even
though GOES-19 is the satellite actually in the GOES-East slot for our
case days: the installed CRTM release (v2.4.0, Feb 2023) predates GOES-19
(launched June 2024, operational as GOES-East April 2025) entirely, so no
`abi_g19` coefficient file exists yet. GOES-16 and GOES-19 carry the same
ABI instrument design, so its spectral response is the right available
proxy -- viewing geometry (satellite position) is still computed for
GOES-19's actual orbital slot via `sat_lon`/`sat_height`, only the
sensor's spectral characteristics come from the GOES-16 coefficients.

This deliberately does not aim for bit-for-bit parity with UPP's internal
choices (its exact vertical grid, cloud-overlap assumptions, etc.) -- the
goal is a scientifically fair comparison using the same forward model and
the right input variables, not a clone of UPP's software. See the project
chat history for that reasoning.

STATUS: working, verified against real WRF output (2026-08-28). Cross-checked
against `wrf.read_cloud_top_temperature` (`ctt`) for the same file/time: r=0.80,
similar value range (both ~215-310K), mean abs difference ~5.4K -- a strong,
physically sane agreement for two independently-implemented methods (this
module's own CRTM forward-model call vs. wrf-python's own diagnostic), not
just "it ran without crashing." Several of the WRF-state -> CRTM-input mapping
choices remain principled-but-unverified (each marked individually with what
would need checking): the geostationary zenith-angle formula (though that
formula itself is independently verified, see its docstring), the
layer-interface pressure approximation, the fixed per-hydrometeor effective
radii, and the ozone climatology placeholder.

`read_simulated_brightness_temperature` now runs at full WRF resolution by
default (memory-bounded via row-chunking, not spatial subsampling -- see its
own docstring and `_run_crtm_chunk`/`_slice_state`), verified (2026-08-28)
two ways: bit-identical to a non-chunked call on a small case (CRTM's
forward model is column-independent, so chunking doesn't change the answer
-- only how much memory a call needs at once), and a real full-resolution
run on the 149-level, 389x549-column domain that OOM'd unchunked -- 56
chunks, 213,561 profiles, ~145 seconds, zero NaN, 215-309K range (mean
296K), consistent with the earlier strided/cross-checked result. Results
are cached to disk (`.npz`, see `cache_dir`) since a full-resolution run is
still slow (tens of chunked CRTM Fortran calls, each paying its own
coefficient-loading overhead) -- confirmed a cached full-resolution result
reloads in ~0.02s vs. the original ~145s compute.

Setup (done once, 2026-08-28, on OLCF Cirrus -- see below for what it took)
----------------------------------------------------------------------------
The short version: `uv run python -m bnf_wrf_eval` (or any code importing
this module) just works once pyCRTM is installed in this project's `.venv`
(see step 3) -- nothing here needs redoing on a machine where that's already
true. The steps below are for rebuilding from scratch (a fresh clone, a
different machine, or after `uv sync` if that ever removes the manually
`uv pip install`-ed `pycrtm-jcsda` package -- it isn't declared in
pyproject.toml since it can only be built against a specific local CRTM
install path, so `uv sync` doesn't know to keep it. Check with
`uv run python -c "import pyCRTM"` after any `uv sync` and redo step 3 below
if it fails).

1. Build CRTM core (already done; lives outside this repo at
   `/ccsopen/home/hengxiao80/crtm`, its own git checkout with a `CLAUDE.md`
   documenting the exact build procedure used -- read that first if
   rebuilding). Two non-obvious requirements, both discovered by real build
   failures, not anticipated in advance:
   - **`-fPIC` is required** in `FCFLAGS` even though CRTM itself is only a
     static archive (`libcrtm.a`) -- without it, linking `libcrtm.a` into
     pyCRTM's shared-object Python extension fails with "relocation
     R_X86_64_32 ... can not be used when making a shared object; recompile
     with -fPIC". Harmless for ordinary static linking, so just always
     include it.
   - **CRTM's I/O endianness must match pyCRTM's assumptions**: CRTM's
     `configure` defaults to big-endian I/O for `ifort`
     (`-convert big_endian`, matching the standard `fix/` coefficient files'
     `Big_Endian` variants) but pyCRTM's own coefficient-linking script
     (`setup.py`'s `linkCoef()`) and its Python-side coefficient-header
     parser (`crtm_io.py`) both hardcode an assumption of **little-endian**
     files -- configure CRTM with `--disable-big-endian` instead, or pyCRTM
     fails confusingly (first with "Data file needs to be byte-swapped",
     then even after manually re-linking to `Big_Endian` coefficient files,
     with "Wrong Endian Binary Found... No Coefficients Valid (litte
     endian)" from pyCRTM's own Python-side parser, which has no
     big-endian code path at all).
   Also needs the exact same Fortran compiler CRTM itself was built with
   (`ifort`, not `gfortran` -- `.mod` files are compiler-specific binary
   formats; mixing them fails with "Fatal Error: Reading module
   'crtm_module': Unexpected EOF"), set via `FC`/`F77`/`F90` env vars
   before building pyCRTM (step 3), since pyCRTM's own CMake-based build
   otherwise picks whatever Fortran compiler `cmake` finds first on `PATH`.
2. CRTM's coefficient files: already downloaded, at
   `/ccsopen/home/hengxiao80/crtm/fix` (see that repo's `CLAUDE.md`). Point
   `pyCRTM`'s `setup.cfg` (`[Coefficients] source_path`) at it; pyCRTM's own
   `setup.py` flattens/symlinks the needed files into `path_used` (currently
   `/ccsopen/home/hengxiao80/crtm/coefficients_flat`) at build time.
3. Build pyCRTM (source at `/ccsopen/home/hengxiao80/crtm/src/pycrtm`,
   cloned from https://github.com/JCSDA/pycrtm) and install it into this
   project's `.venv`:
   ```
   module load intel/20.0.4 openmpi/4.1.4 netcdf-c/4.8.1 netcdf-fortran/4.5.3
   cd /ccsopen/home/hengxiao80/crtm/src/pycrtm
   export SETUPTOOLS_USE_DISTUTILS=stdlib  # same numpy.distutils workaround wrf-python needed
   export FC=$(which ifort) F77=$(which ifort) F90=$(which ifort)
   uv pip install --python <this-repo>/.venv/bin/python3.11 --no-build-isolation .
   ```
   Needs `scikit-build`, `cmake`, and `ninja` (installed here as project dev
   dependencies, since pyCRTM's own build uses scikit-build/CMake/f2py) and
   `h5py` (a runtime import of pyCRTM's own `pyCRTM.py`, not declared as its
   own dependency anywhere -- added as a regular project dependency here).
   `--no-build-isolation` matters: pyCRTM's `setup.py` needs the CRTM/coeff
   paths from its own `setup.cfg`, not a clean isolated build env.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import netCDF4
import numpy as np

# GOES-16's ABI coefficients -- checked against the actual installed CRTM
# fix files (2026-08-28): only abi_g16/abi_g17/abi_gr exist, no abi_g19,
# since this CRTM release (v2.4.0) predates GOES-19 entirely. See the
# module docstring for why abi_g16 is the right stand-in.
SENSOR_ID = "abi_g16"

# ABI channel 13, 10.3 micron clean IR window -- same channel goes.py reads
# observed Tb for, and the same one hrrr.read_simulated_brightness_temperature
# targets via SBT114. pyCRTM's `runDirect()` computes all of a sensor's
# channels at once; this is just which output column to keep.
CHANNEL = 13

# ABI channel number -> approximate central wavelength (microns), for the
# channels this coefficient set actually carries. CONFIRMED via real
# testing that abi_g16's SpcCoeff only has 10 channels, not ABI's full 16
# -- its own `loadInst()` numbers them 1-10 in coefficient-file order,
# which is NOT the same as the true ABI channel number (index != channel
# - 1, a bug caught by testing: reading index `channel - 1` = 12 raised
# IndexError since there are only 10 columns). The coefficient file's
# actual per-channel wavelengths (`pyCRTM.wavelengthMicrons`, read at
# runtime) turned out to be exactly the 10 IR channels (7-16) in order, so
# the right output column is found by nearest-wavelength match against
# this table rather than assumed from the channel number directly.
_ABI_CHANNEL_WAVELENGTH_UM = {
    7: 3.9, 8: 6.2, 9: 6.9, 10: 7.3, 11: 8.4,
    12: 9.6, 13: 10.3, 14: 11.2, 15: 12.3, 16: 13.3,
}

# GOES-19 East operational sub-satellite longitude and altitude. UNVERIFIED
# here as a hardcoded default -- prefer passing the actual values read from
# a downloaded GOES CMIP file's own `goes_imager_projection` metadata (see
# `goes.read_brightness_temperature`, which already extracts
# `longitude_of_projection_origin`/`perspective_point_height`) via the
# `sat_lon`/`sat_height` parameters below, since station-keeping can shift
# these slightly over a satellite's life.
DEFAULT_SAT_LON = -75.2
DEFAULT_SAT_HEIGHT = 35_786_023.0  # meters above the reference ellipsoid
EARTH_RADIUS = 6_371_000.0  # meters, spherical approximation

# Fixed per-hydrometeor-category effective radii (microns), following the
# documented, precedented simplification used in published all-sky
# assimilation work (rather than a microphysics-scheme-specific formula,
# which is closer to what UPP itself does but isn't required for a fair
# comparison here -- see project chat history). CONFIRMED as a legitimate
# approach via research; these specific numbers should be treated as a
# reasonable starting point, not something re-derived from a primary
# source here.
EFFECTIVE_RADIUS_UM = {
    "liquid": 20.0,
    "ice": 40.0,
    "rain": 400.0,
    "snow": 600.0,
    "graupel": 800.0,
}

# WRF hydrometeor mixing-ratio variable -> CRTM cloud category.
_HYDROMETEOR_VARS = ("QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP")
_HYDROMETEOR_CATEGORY = {
    "QCLOUD": "liquid",
    "QRAIN": "rain",
    "QICE": "ice",
    "QSNOW": "snow",
    "QGRAUP": "graupel",
}

# CRTM's cloud category integers -- confirmed directly from the installed
# CRTM source (src/Atmosphere/Cloud/CRTM_Cloud_Define.f90:107-114:
# WATER_CLOUD=1, ICE_CLOUD=2, RAIN_CLOUD=3, SNOW_CLOUD=4, GRAUPEL_CLOUD=5,
# HAIL_CLOUD=6), not a guess.
CLOUD_TYPE_ID = {
    "liquid": 1,
    "ice": 2,
    "rain": 3,
    "snow": 4,
    "graupel": 5,
}

# CRTM's own documented default values for the 6 categorical surface-type
# slots (`profiles.surfaceTypes[:, 0:6]` = land/soil/vegetation/water/
# snow/ice type) -- confirmed directly from the installed CRTM source
# (src/Surface/CRTM_Surface_Define.f90:165-187, DEFAULT_LAND_TYPE etc.,
# all defaulting to 1, "first item in list"). We don't have a real
# land-use classification to map WRF's own vegetation category onto
# CRTM's, so this uses CRTM's own defaults everywhere rather than
# guessing a mapping (real UPP source for HRRR *does* do this mapping --
# WRF's own IVGTYP field, using WRF's IGBP-based MODIS-NOAH categories,
# through a lookup table to CRTM's named IR land-type categories -- a
# known, scoped improvement, not yet done here).
#
# For Water_Type/Snow_Type/Ice_Type specifically, "1" isn't a placeholder
# guess -- confirmed from the CRTM v2.4.0 User Guide (Table 4.18, "Water,
# snow, and ice surface subtypes for infrared and visible sensors", which
# is what applies to our ABI IR channel; microwave sensors don't consult
# these fields at all): for IR/VIS, CRTM ships exactly
#   Water_Type: 1 = "sea water" (the *only* IR/VIS water category -- there
#     is no separate fresh-water IR emissivity dataset to pick; the
#     fresh-vs-salt distinction is `Salinity`, a continuous physical input
#     we already set to 0.0 for this inland domain, not a Water_Type index)
#   Snow_Type: 1 = "old snow", 2 = "new snow"
#   Ice_Type: 1 = "new ice" (again the *only* IR/VIS ice category -- no
#     sea-ice-vs-freshwater-ice split exists here either)
# so `1` for all three is simply CRTM's one-or-only-default IR/VIS choice,
# not an unverified stand-in.
DEFAULT_SURFACE_TYPE = 1

_GAS_CONSTANT_DRY_AIR = 287.05  # J/(kg K)
_GRAVITY = 9.81  # m/s^2, matches wrf-python's convention closely enough here


def _get_pycrtm():
    """Import pyCRTM lazily -- see the module docstring's "Setup" section
    for what has to be built/installed first. Not yet attempted."""
    from pyCRTM import profilesCreate, pyCRTM  # noqa: N813

    return pyCRTM, profilesCreate


def _geostationary_zenith_angle(
    lon: np.ndarray, lat: np.ndarray, sat_lon: float, sat_height: float
) -> np.ndarray:
    """Satellite (sensor) zenith angle, in degrees, seen from a
    geostationary satellite at `sat_lon` (degrees east) and `sat_height`
    (meters above the surface) looking at ground points at `lon`/`lat`.

    Standard geostationary-viewing-geometry formula (spherical Earth, no
    oblateness correction), derived from the law of cosines/sines in the
    triangle formed by Earth's center, the ground point, and the
    satellite. Checked against two independent references for the BNF
    site (34.34N, 87.34W) viewed from GOES-19 (-75.2E): the closed-form
    result here (~42.0 deg) matches a direct arccos cross-check computed
    a different way, and the formula's horizon limit (where the
    denominator crosses zero) lands at an angular distance of ~81.3 deg
    from the sub-satellite point, matching the well-known GEO visibility
    limit -- so the derivation itself is solid. What's still unverified
    is the input side: this hasn't been cross-checked against a real
    per-pixel zenith angle from an actual GOES product file. Nothing in
    pyCRTM computes this itself (confirmed via research); its test cases
    just read a precomputed value from a data file.
    """
    lon_rad = np.deg2rad(lon - sat_lon)
    lat_rad = np.deg2rad(lat)

    # Angular distance from the sub-satellite point (spherical Earth).
    cos_psi = np.cos(lat_rad) * np.cos(lon_rad)
    sin_psi = np.sqrt(np.clip(1.0 - cos_psi**2, 0.0, 1.0))

    orbital_radius = EARTH_RADIUS + sat_height
    zenith_rad = np.arctan2(
        orbital_radius * sin_psi, orbital_radius * cos_psi - EARTH_RADIUS
    )
    return np.rad2deg(zenith_rad)


def _layer_pressures(p_mass: np.ndarray, psfc: np.ndarray, p_top: float) -> np.ndarray:
    """Interface (level) pressures (Pa) from WRF's layer (mass-point)
    pressures, for CRTM's `Pi` (level_pressure) input.

    UNVERIFIED approximation: WRF's wrfout carries pressure at layer
    midpoints only (`P` + `PB`, on the unstaggered `bottom_top` grid), not
    at the staggered layer interfaces CRTM also wants. This takes the
    geometric mean of vertically-adjacent layer pressures as the interface
    value between them, with the surface interface set to `PSFC` and the
    model-top interface set to the run's `P_TOP` global attribute -- a
    common, reasonable approximation, but not verified against a more
    careful hydrostatic integration.

    `p_mass` has shape (nz, ny, nx), top-of-atmosphere-first (see
    `read_simulated_brightness_temperature`, which reverses WRF's native
    bottom-up ordering before calling this). Returns shape (nz+1, ny, nx).
    """
    interior = np.sqrt(p_mass[:-1] * p_mass[1:])
    top = np.full((1, *p_mass.shape[1:]), p_top)
    return np.concatenate([top, interior, psfc[np.newaxis, ...]], axis=0)


def _ozone_climatology_ppmv(pressure_hpa: np.ndarray) -> np.ndarray:
    """Rough ozone volume mixing ratio (ppmv) as a function of pressure,
    for CRTM's `O3` input.

    This is required, not optional: `O3` is unconditionally one of
    pyCRTM's core profile fields (`profilesCreate`'s default `keys =
    ['P','T','Q','O3']`, confirmed by reading its source), always passed
    through as a real trace-gas profile -- there's no automatic
    climatology fallback for it the way `profiles.climatology` (left at
    pyCRTM's own default, US Standard Atmosphere) provides for gases
    *not* included in the profile at all. WRF doesn't predict ozone (no
    chemistry option in these runs), so this is a placeholder climatology,
    not a real profile: a low, roughly constant tropospheric value with a
    simple stratospheric peak. Deliberately rough, not a verified US
    Standard Atmosphere fit -- worth replacing with a proper reference
    profile (CRTM ships some as part of its own test/reference data) if
    channel-13 Tb turns out to be more ozone-sensitive than expected for a
    clean IR window channel (it usually isn't, much more relevant to
    channels nearer the 9.6 micron ozone band).
    """
    tropospheric = 0.03
    stratospheric_peak = 8.0
    # Roughly locate the peak near 25 hPa-ish (crude, not a real profile).
    peak = stratospheric_peak * np.exp(-((np.log(pressure_hpa) - np.log(25.0)) ** 2) / 2.0)
    return np.maximum(tropospheric, peak)


def _read_wrf_state(wrf_file: str | Path) -> dict:
    """Read the raw 3D/2D WRF fields the CRTM mapping needs, without going
    through `wrf.py`'s higher-level (already-single-level) readers."""
    with netCDF4.Dataset(wrf_file) as nc:
        p = np.ma.filled(nc.variables["P"][0], np.nan) + np.ma.filled(
            nc.variables["PB"][0], np.nan
        )
        t_pot = np.ma.filled(nc.variables["T"][0], np.nan) + 300.0  # perturbation potential temp
        t = t_pot * (p / 100000.0) ** (_GAS_CONSTANT_DRY_AIR / 1004.0)  # -> actual temperature (K)
        qvapor = np.ma.filled(nc.variables["QVAPOR"][0], np.nan)
        hydrometeors = {
            var: np.ma.filled(nc.variables[var][0], 0.0) if var in nc.variables else None
            for var in _HYDROMETEOR_VARS
        }
        cldfra = np.ma.filled(nc.variables["CLDFRA"][0], 0.0) if "CLDFRA" in nc.variables else None
        ph = np.ma.filled(nc.variables["PH"][0], np.nan)
        phb = np.ma.filled(nc.variables["PHB"][0], np.nan)
        height = (ph + phb) / _GRAVITY  # geopotential height at layer interfaces (m)
        psfc = np.ma.filled(nc.variables["PSFC"][0], np.nan)
        tsk = np.ma.filled(nc.variables["TSK"][0], np.nan)
        landmask = np.ma.filled(nc.variables["LANDMASK"][0], np.nan)
        hgt = np.ma.filled(nc.variables["HGT"][0], 0.0)  # terrain elevation, m
        u10 = np.ma.filled(nc.variables["U10"][0], 0.0)
        v10 = np.ma.filled(nc.variables["V10"][0], 0.0)
        lon = np.ma.filled(nc.variables["XLONG"][0], np.nan)
        lat = np.ma.filled(nc.variables["XLAT"][0], np.nan)
        p_top = float(nc.variables["P_TOP"][0]) if "P_TOP" in nc.variables else float(nc.P_TOP)
        time_str = b"".join(nc.variables["Times"][0, :]).decode("utf-8")
        valid_time = dt.datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")

    return dict(
        p=p, t=t, qvapor=qvapor, hydrometeors=hydrometeors, cldfra=cldfra, height=height,
        psfc=psfc, tsk=tsk, landmask=landmask, hgt=hgt, u10=u10, v10=v10,
        lon=lon, lat=lat, p_top=p_top, valid_time=valid_time,
    )


def _slice_state(state: dict, row_sel: slice, col_sel: slice) -> dict:
    """Return a copy of a `_read_wrf_state()` dict with every ny/nx-shaped
    field sliced by `[..., row_sel, col_sel]`.

    Used both for `stride` subsampling (a strided `row_sel`/`col_sel`) and
    for splitting the domain into row-bands for memory-bounded CRTM calls
    (a contiguous `row_sel`, full `col_sel`) -- CRTM is a column-independent
    forward model (no cross-profile coupling), so slicing the domain either
    way and stitching/subsampling the per-slice results back together gives
    the same output a single whole-domain call would, just computed in
    smaller (or fewer) pieces.
    """
    sliced = dict(state)
    for key in ("p", "t", "qvapor", "height", "psfc", "tsk", "landmask",
                "hgt", "u10", "v10", "lon", "lat"):
        sliced[key] = state[key][..., row_sel, col_sel]
    sliced["hydrometeors"] = {
        var: (data[..., row_sel, col_sel] if data is not None else None)
        for var, data in state["hydrometeors"].items()
    }
    if state["cldfra"] is not None:
        sliced["cldfra"] = state["cldfra"][..., row_sel, col_sel]
    return sliced


def _cache_path(
    wrf_file: str | Path, cache_dir: str | Path, sensor_id: str, channel: int, stride: int
) -> Path:
    """Deterministic cache filename for one `read_simulated_brightness_temperature`
    call, derived from `wrf_file`'s own path structure (its run directory's
    last two path components, matching `plotting._output_path`'s
    convention, plus the wrfout filename itself, which already embeds the
    domain and valid time) and the parameters that change the result's
    shape/values (`sensor_id`, `channel`, `stride`). Does NOT depend on
    `max_profiles_per_chunk` (chunking is purely how the computation is
    split up, not what's computed) or `sat_lon`/`sat_height`/
    `coefficient_path` (rarely varied per call; pass `recompute=True` if
    you do change one of those and want a fresh cache entry rather than a
    stale one)."""
    wrf_file = Path(wrf_file)
    run_dir = wrf_file.parent
    run_label = f"{run_dir.parent.name}_{run_dir.name}"
    return Path(cache_dir) / f"{run_label}_{wrf_file.name}_{sensor_id}_ch{channel}_stride{stride}.npz"


def _load_cache(cache_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    with np.load(cache_path) as npz:
        lon, lat, values = npz["lon"], npz["lat"], npz["values"]
        valid_time = dt.datetime.fromisoformat(str(npz["valid_time"]))
    return lon, lat, values, valid_time


def _save_cache(
    cache_path: Path, lon: np.ndarray, lat: np.ndarray, values: np.ndarray, valid_time: dt.datetime
) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path, lon=lon, lat=lat, values=values, valid_time=np.array(valid_time.isoformat())
    )


def _run_crtm_chunk(
    pyCRTM,  # noqa: N803
    profilesCreate,  # noqa: N803
    state: dict,
    sensor_id: str,
    channel: int,
    sat_lon: float,
    sat_height: float,
    coefficient_path: str | Path | None,
) -> np.ndarray:
    """Build CRTM profiles from one (possibly row-chunked) `state` dict and
    run the forward model, returning simulated Tb (K) with shape (ny, nx)
    matching `state["lon"]`. Every choice below is explained/justified in
    `read_simulated_brightness_temperature`'s module-level context; see the
    module docstring's STATUS note for what's been verified."""
    ny, nx = state["lon"].shape
    nz = state["p"].shape[0]
    n_profiles = ny * nx

    def flatten_top_down(field_3d: np.ndarray) -> np.ndarray:
        # WRF's bottom_top dimension is surface-first; CRTM's layers are
        # top-of-atmosphere-first (standard RT-model convention).
        return field_3d[::-1].reshape(nz, n_profiles).T

    p_layer = flatten_top_down(state["p"]) / 100.0  # Pa -> hPa
    t_layer = flatten_top_down(state["t"])
    q_layer = flatten_top_down(state["qvapor"]) * 1000.0  # kg/kg -> g/kg

    height_interfaces = state["height"][::-1].reshape(nz + 1, n_profiles).T
    layer_thickness = height_interfaces[:, :-1] - height_interfaces[:, 1:]  # (n_profiles, nz)

    p_interfaces_pa = _layer_pressures(
        state["p"][::-1], state["psfc"], state["p_top"]
    )
    p_interfaces_hpa = p_interfaces_pa.reshape(nz + 1, n_profiles).T / 100.0

    air_density = (p_layer * 100.0) / (_GAS_CONSTANT_DRY_AIR * t_layer)  # kg/m^3

    n_clouds = sum(1 for v in state["hydrometeors"].values() if v is not None)
    # nAerosols=0: cleanly excludes aerosol effects (confirmed via
    # pyCRTM's own source: with nAerosols=0, the `aerosols`/`aerosolType`
    # keys are never created, and runDirect() only passes them through if
    # present -- rather than passing default/unset values through an
    # aerosol path we haven't set up).
    profiles = profilesCreate(
        nProfiles=n_profiles, nLevels=nz, nAerosols=0, nClouds=max(n_clouds, 1)
    )

    zenith = _geostationary_zenith_angle(state["lon"], state["lat"], sat_lon, sat_height)
    zenith_flat = zenith.reshape(n_profiles)
    profiles.Angles[:, 0] = zenith_flat  # satellite/instrument zenith angle
    # Instrument scan angle (index 4) differs in principle from the
    # ground-point zenith angle above (it's the same triangle's angle at
    # the satellite, not at the ground point) -- approximated here as
    # equal to the zenith angle, a simplification, not a separately
    # derived/verified formula.
    profiles.Angles[:, 4] = zenith_flat
    # profilesCreate() defaults the remaining Angles columns (sensor
    # azimuth, sun zenith, sun azimuth) to NaN. Confirmed via real testing
    # that leaving them NaN poisons the whole computation (Bt came back
    # all-NaN with no error) -- CRTM's Fortran code evidently touches them
    # even for a thermal-IR-only channel that shouldn't care about solar
    # geometry. pyCRTM's own reference test cases
    # (testCases/test_atms_no_clouds.py) set sun zenith to 100 degrees
    # (below-horizon, i.e. "night", so no solar contribution) as their
    # standard placeholder when solar geometry isn't being computed
    # properly for the case; sensor/sun azimuth are set to arbitrary
    # placeholders since they don't matter without solar contribution.
    profiles.Angles[:, 1] = 0.0  # sensor azimuth -- placeholder, unused here
    profiles.Angles[:, 2] = 100.0  # sun zenith -- below horizon, "night"
    profiles.Angles[:, 3] = 0.0  # sun azimuth -- placeholder, unused here

    profiles.SurfGeom[:, 0] = state["lat"].reshape(n_profiles)
    # CRTM_Geometry_IsValid requires longitude in [0, 360) East (confirmed
    # via real testing: it rejects WRF's native signed [-180, 180]
    # convention with "Invalid longitude"), so convert.
    profiles.SurfGeom[:, 1] = state["lon"].reshape(n_profiles) % 360.0
    profiles.SurfGeom[:, 2] = state["hgt"].reshape(n_profiles)

    profiles.DateTimes[:, 0] = state["valid_time"].year
    profiles.DateTimes[:, 1] = state["valid_time"].month
    profiles.DateTimes[:, 2] = state["valid_time"].day

    # `profiles` is a namedtuple (confirmed from pyCRTM's own source), so
    # its attributes can't be rebound (`profiles.P = ...` would raise
    # AttributeError) -- every field must be mutated in place via `[...]`
    # into the array profilesCreate already allocated.
    profiles.P[...] = p_layer
    profiles.Pi[...] = p_interfaces_hpa
    profiles.T[...] = t_layer
    profiles.Q[...] = q_layer
    profiles.O3[...] = _ozone_climatology_ppmv(p_layer)

    # surfaceFractions/surfaceTemperatures are both (n_profiles, 4):
    # [land, water, snow, ice] -- confirmed from pyCRTM's own source
    # (profilesCreate + the runDirect() call signature). WRF's LANDMASK
    # is a binary land/water split with no snow/ice fraction, so this
    # keeps it simple: 100% of whichever LANDMASK says, 0% snow/ice (a
    # simplification -- WRF does track snow cover separately via SNOWC,
    # not incorporated here). TSK fills both the land and water slots
    # since only the one with nonzero fraction actually matters.
    is_land = state["landmask"].reshape(n_profiles) > 0.5
    fractions = np.zeros((n_profiles, 4))
    fractions[is_land, 0] = 1.0
    fractions[~is_land, 1] = 1.0
    profiles.surfaceFractions[:, :] = fractions
    profiles.surfaceTemperatures[:, 0] = state["tsk"].reshape(n_profiles)
    profiles.surfaceTemperatures[:, 1] = state["tsk"].reshape(n_profiles)

    # surfaceTypes is (n_profiles, 6): land/soil/vegetation/water/snow/ice
    # classification indices. We don't have a real WRF-vegetation-category
    # -> CRTM-category mapping, so this uses CRTM's own documented
    # defaults throughout (see DEFAULT_SURFACE_TYPE) rather than inventing
    # one.
    profiles.surfaceTypes[:, :] = DEFAULT_SURFACE_TYPE

    # profilesCreate() defaults Salinity to NaN (`np.nan * np.zeros(...)`,
    # confirmed from its source), which would feed NaN into CRTM's ocean
    # surface-emissivity (FASTEM) Fortran code for any water-covered
    # profile. This domain (BNF, inland Alabama) has no true ocean --
    # any water pixels are rivers/lakes -- so 0 PSU (fresh water) is the
    # physically right value, not CRTM's usual ~33 PSU seawater assumption.
    profiles.Salinity[:] = 0.0

    wind_speed = np.sqrt(state["u10"] ** 2 + state["v10"] ** 2)
    profiles.windSpeed10m[:] = wind_speed.reshape(n_profiles)
    # windDirection10m left at its default (0) -- primarily affects ocean
    # surface roughness/foam in the microwave, not this LW window channel,
    # and this domain is mostly land.

    cloud_idx = 0
    total_water_content = np.zeros((n_profiles, nz))
    for var, data in state["hydrometeors"].items():
        if data is None:
            continue
        category = _HYDROMETEOR_CATEGORY[var]
        mixing_ratio = flatten_top_down(data)  # kg/kg
        water_content = mixing_ratio * air_density * layer_thickness  # kg/m^2 per layer
        profiles.clouds[:, :, cloud_idx, 0] = water_content
        profiles.clouds[:, :, cloud_idx, 1] = EFFECTIVE_RADIUS_UM[category]
        profiles.cloudType[:, cloud_idx] = CLOUD_TYPE_ID[category]
        total_water_content += water_content
        cloud_idx += 1

    # profilesCreate() defaults cloudFraction to all-zero. Left at zero
    # alongside nonzero cloud water content, CRTM's cloud optics produced
    # all-NaN Bt with no error message (confirmed via real testing --
    # in-cloud water content is presumably grid-mean-water /
    # cloud-fraction internally, so a zero fraction is a divide-by-zero).
    # Use WRF's own CLDFRA diagnostic where available; fall back to a
    # binary "any hydrometeor present at this level" fraction otherwise.
    if state["cldfra"] is not None:
        profiles.cloudFraction[:, :] = flatten_top_down(state["cldfra"])
    else:
        profiles.cloudFraction[:, :] = (total_water_content > 0.0).astype(float)

    crtm_ob = pyCRTM()
    crtm_ob.profiles = profiles
    crtm_ob.sensor_id = sensor_id
    if coefficient_path is not None:
        # pyCRTM.__init__ hardcodes self.coefficientPath from the
        # `path_used` baked into pycrtm_setup.txt at build time (confirmed
        # via its source -- there's no env var it reads instead); override
        # the attribute directly if a different location is wanted.
        crtm_ob.coefficientPath = str(coefficient_path) + "/"
    crtm_ob.loadInst()

    # Find which output column is `channel`, by nearest-wavelength match
    # rather than assuming column index == channel - 1 (verified via real
    # testing to be wrong for abi_g16 -- see _ABI_CHANNEL_WAVELENGTH_UM).
    target_um = _ABI_CHANNEL_WAVELENGTH_UM[channel]
    channel_index = int(np.argmin(np.abs(crtm_ob.wavelengthMicrons - target_um)))

    crtm_ob.runDirect()

    tb_flat = crtm_ob.Bt[:, channel_index]
    return tb_flat.reshape(ny, nx)


def read_simulated_brightness_temperature(
    wrf_file: str | Path,
    sensor_id: str = SENSOR_ID,
    channel: int = CHANNEL,
    sat_lon: float = DEFAULT_SAT_LON,
    sat_height: float = DEFAULT_SAT_HEIGHT,
    coefficient_path: str | Path | None = None,
    stride: int = 1,
    max_profiles_per_chunk: int = 4000,
    cache_dir: str | Path | None = "crtm_cache",
    use_cache: bool = True,
    recompute: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """CRTM-simulated brightness temperature (K) computed directly from
    WRF's own model state, for `sensor_id`/`channel` (GOES-19 ABI channel
    13 by default -- see module docstring).

    Prefer passing `sat_lon`/`sat_height` from a real downloaded GOES CMIP
    file's projection metadata (`goes.read_brightness_temperature` already
    extracts these) rather than relying on the hardcoded defaults, since
    those are unverified against the specific date being compared.

    Runs at full WRF resolution by default, in memory-bounded row-chunks
    (see `max_profiles_per_chunk`) rather than one whole-domain CRTM call --
    CRTM runs one profile per WRF column with no batching, so a
    whole-domain call was confirmed (via real testing) to exceed a
    login-node's memory budget on this project's real domains (e.g. 149
    vertical levels x 389x549 columns, ~213k profiles, killed by the OS's
    memory cgroup at ~27GB resident). Chunking is numerically exact --
    CRTM's forward model is column-independent, so splitting the domain
    into row-bands and stitching the results back together gives the same
    answer a single whole-domain call would, just slower (each chunk pays
    its own CRTM coefficient-loading overhead -- expect several dozen
    "Initializing/Destroying the CRTM" banners on a large domain) and with
    bounded memory instead of unbounded.

    `stride`: also subsample the WRF horizontal grid every `stride` points
    in both directions before running CRTM (default 1: every point, no
    subsampling) -- independent of and composable with the chunking above,
    e.g. for a quick coarse preview. Memory/time scale with
    `(ny/stride) * (nx/stride) * nz`.

    `max_profiles_per_chunk`: row-chunk size, chosen so
    `chunk_rows * nx <= max_profiles_per_chunk` (at least 1 row). The
    default (4000) is calibrated to the one domain size verified safe by
    real testing (with a safety margin below the ~213k-profile figure that
    actually OOM'd); lower it if this still runs out of memory on a
    particular machine/domain, raise it (fewer, larger chunks) to trade
    memory headroom for fewer chunks and less repeated coefficient-loading
    overhead.

    `cache_dir`: if not None (the default, `"crtm_cache"` relative to
    wherever this is run from), save the result to a `.npz` file there
    (see `_cache_path` for the naming convention) and reuse it on a later
    call with the same `wrf_file`/`sensor_id`/`channel`/`stride` instead of
    recomputing -- full-resolution CRTM runs are expensive (many chunked
    Fortran calls), so this avoids repeating that work across repeated
    plotting/notebook runs. Set `use_cache=False` to ignore any existing
    cache file without overwriting it, or `recompute=True` to force
    recomputing and overwrite it. The cache key does NOT include
    `sat_lon`/`sat_height`/`coefficient_path` -- pass `recompute=True` if
    you change one of those and want a fresh entry rather than a stale one.

    Returns (lon, lat, values, valid_time); lon/lat/values share the
    (possibly strided) WRF mass grid's shape (south_north, west_east).
    """
    cache_path = None
    if cache_dir is not None:
        cache_path = _cache_path(wrf_file, cache_dir, sensor_id, channel, stride)
        if use_cache and not recompute and cache_path.exists():
            return _load_cache(cache_path)

    pyCRTM, profilesCreate = _get_pycrtm()  # noqa: N806

    state = _read_wrf_state(wrf_file)
    if stride > 1:
        state = _slice_state(state, slice(None, None, stride), slice(None, None, stride))
    ny, nx = state["lon"].shape

    chunk_rows = max(1, max_profiles_per_chunk // max(nx, 1))
    if chunk_rows >= ny:
        values = _run_crtm_chunk(
            pyCRTM, profilesCreate, state, sensor_id, channel, sat_lon, sat_height, coefficient_path
        )
    else:
        values = np.empty((ny, nx), dtype=float)
        for start in range(0, ny, chunk_rows):
            end = min(start + chunk_rows, ny)
            chunk_state = _slice_state(state, slice(start, end), slice(None))
            values[start:end, :] = _run_crtm_chunk(
                pyCRTM, profilesCreate, chunk_state, sensor_id, channel,
                sat_lon, sat_height, coefficient_path,
            )

    result = (state["lon"], state["lat"], values, state["valid_time"])
    if cache_path is not None:
        _save_cache(cache_path, *result)
    return result
