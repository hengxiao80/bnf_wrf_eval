# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project state

No longer a bare scaffold: `src/bnf_wrf_eval/` has a working WRF-vs-HRRR-vs-observations comparison pipeline
(reflectivity and brightness temperature, described below), used from the notebooks under `notebooks/`. Treat
existing modules as the conventions to follow (docstring style, the `(lon, lat, values, valid_time)` return
convention every reader function uses, lazy imports of heavy/optional dependencies, etc.) rather than starting
from scratch.

## Purpose

This repo evaluates WRF (Weather Research and Forecasting model) simulations at the BNF site (DOE ARM Bankhead
National Forest site) that were performed by a colleague. The run outputs are linked at `satoshi_testruns`.
Evaluation is done two ways:

1. Against the analyses/forecasts used to initialize and drive the runs (HRRR, ERA5, EDA, ...), linked at
   `satoshi_forcing_data`.
2. Against independent observations (satellite, radar, sounding, ...) — currently GOES-R ABI brightness
   temperature and MRMS composite reflectivity, downloaded on demand (see below).

## Evaluation approach

Two comparisons -- reflectivity (3 columns: WRF | HRRR | independent observation) and brightness temperature (3
columns, or a 4-panel 2x2 grid -- see below) -- all reprojected onto the same Lambert Conformal projection (taken
from the WRF domain) and the same map extent, with the WRF domain's footprint outlined in red and the BNF site
marked with a red star on every panel:

- **Reflectivity** (`plotting.plot_refl_comparison` / `plot_run_refl_comparison`): WRF's column-max simulated
  `REFL_10CM` vs. HRRR's composite reflectivity (`refc`) vs. MRMS's composite reflectivity
  (`MergedReflectivityQCComposite`). All three are the same physical quantity (dBZ), so this is a direct,
  apples-to-apples comparison.
- **Brightness/cloud-top temperature** (`plotting.plot_tb_comparison` / `plot_run_tb_comparison`, 3 columns; or
  `plot_tb_comparison_4panel` / `plot_run_tb_comparison_4panel`, a 2x2 grid -- see "Structure" below): WRF's own
  CRTM-derived simulated brightness temperature (`crtm.read_simulated_brightness_temperature`, full resolution)
  vs. HRRR's CRTM-derived simulated brightness temperature (`SBT114`) vs. GOES ABI channel 13 (10.3 micron)
  observed brightness temperature -- this is now a genuine apples-to-apples comparison (same forward model, same
  channel) rather than relying on WRF's simpler `ctt` diagnostic (`wrf.read_cloud_top_temperature`, still
  available and still WRF's own purpose-built approximation of the same quantity via `wrf-python`). The 4-panel
  layout additionally shows `ctt` alongside the CRTM panel, for a direct look at how the two WRF-side methods
  compare; the notebooks under `notebooks/` currently use this 4-panel version. Both Tb layouts share a custom
  colormap (`plotting.DEFAULT_TB_CMAP`) matching the conventional grayscale-above/rainbow-below IR enhancement
  (black-to-white for 240-315 K, then a cyan-to-magenta rainbow for the coldest convective cloud tops down to
  180 K) rather than a plain grayscale.

The `plot_run_*` wrappers take a time plus a run's name/directory and locate the matching wrfout/HRRR/MRMS/GOES
files themselves (`auto_download_mrms`/`auto_download_goes` control whether missing observation files get
downloaded automatically or raise). See the docstrings in `src/bnf_wrf_eval/plotting.py` for the full parameter
list and file-layout conventions each one expects.

## Environment and commands

Package management is via `uv` (build backend `uv_build`, package layout under `src/`).

- Install/sync dependencies: `uv sync` -- **on a fresh clone/cache, the
  first sync must be run as `SETUPTOOLS_USE_DISTUTILS=stdlib uv sync`**.
  `wrf-python` builds compiled extensions via the legacy `numpy.distutils`,
  which is incompatible with modern setuptools' vendored `_distutils`;
  this env var (only meaningful on Python <3.12, which still ships a real
  stdlib `distutils`) makes setuptools use that instead. Once built once,
  the wheel is cached and plain `uv sync` works from then on. See the
  `[tool.uv]` tables in `pyproject.toml` for the other wrf-python
  packaging workarounds this needed (a hard, unconditional `basemap`
  dependency that conflicts with this project's matplotlib/numpy
  versions; a missing `numpy` build dependency). `wrf-python` is
  unmaintained on PyPI (stuck at a ~2021 release; the real latest,
  v1.4.2, only exists on GitHub/conda-forge) and needed two more
  workarounds applied directly in `src/bnf_wrf_eval/wrf.py`
  (`_get_wrf_getvar`): a `numpy.float_` shim (removed in NumPy 2.0), and
  a monkeypatch for a version-drift bug where modern `netCDF4.Dataset`
  being iterable makes wrf-python misdetect a single open file as a
  multi-file sequence. Confirmed via diffing against v1.4.2's source: the
  actual `ctt` computation (`g_ctt.py`, `fortran/wrf_fctt.f90`,
  `fortran/wrf_constants.f90`) is unchanged, so these are packaging-only
  workarounds, not something affecting the numbers.
- Add a dependency: `uv add <package>`
- Run the CLI entry point: `uv run bnf-wrf-eval` (maps to `bnf_wrf_eval:main`)
- Run any script/module in the project's environment: `uv run python -m <module>`
- Build the package: `uv build`
- Run a notebook end-to-end from the CLI (e.g. to verify a change): `uv run jupyter nbconvert --to notebook
  --execute --inplace notebooks/<name>.ipynb`

Python version is pinned to 3.11 (`.python-version`, `requires-python = ">=3.11"` in `pyproject.toml`).

No test runner or formatter is configured; `[dependency-groups] dev` has `ipykernel`/`nbconvert` for running and
executing notebooks, and `scikit-build`/`cmake`/`ninja` for building `pyCRTM` from source (see `crtm.py`'s
docstring for the full recipe -- `pyCRTM` itself isn't a declared dependency since it can only be built against a
specific local CRTM install path outside this repo).

If VS Code's Jupyter notebook UI ever hangs on a cell that never seems to execute (kernel process itself is idle
and healthy if you check, e.g. `ps`/connecting a second `jupyter_client` to its connection file), check the
installed `ms-toolsai.jupyter`/`ms-python.*` extension versions (`code --list-extensions --show-versions`)
against the marketplace -- we hit this once with a ~10-month-stale Jupyter extension; updating
(`code --install-extension ms-toolsai.jupyter@<latest> --force`) and reloading the window fixed it.

## Structure

- `src/bnf_wrf_eval/` -- the package (`src` layout).
  - `wrf.py` -- reads wrfout NetCDF fields (`read_field`, `read_column_max_reflectivity`,
    `read_cloud_top_temperature`) and builds the domain's cartopy Lambert Conformal projection
    (`get_lambert_projection`).
  - `hrrr.py` -- reads HRRR native-level analysis GRIB2 fields by GRIB2 metadata (`read_composite_reflectivity`,
    `read_olr`, `read_simulated_brightness_temperature`). Lazily bootstraps `eccodes`/`ecmwflibs` on first actual
    use (see `_eccodes_setup.py`) rather than at import time.
  - `goes.py` -- finds/downloads/reads GOES-R ABI Cloud and Moisture Imagery Product (CMIP) brightness
    temperature from the public `noaa-goes19` AWS Open Data S3 bucket (anonymous, no credentials).
  - `mrms.py` -- finds/downloads/reads MRMS composite reflectivity from the public `noaa-mrms-pds` AWS Open Data
    S3 bucket.
  - `crtm.py` -- **working, verified against real WRF output** (2026-08-28): computes simulated GOES-East ABI
    channel 13 Tb directly from WRF's own state via CRTM (through JCSDA's `pyCRTM`), the same forward-model
    approach UPP uses for HRRR's `SBT114` -- for a fairer WRF-vs-HRRR-vs-GOES comparison than
    `wrf.read_cloud_top_temperature`'s `ctt` diagnostic. Cross-checked against `ctt` for the same file/time:
    r=0.80, similar value range, ~5.4K mean absolute difference -- strong agreement for two independently
    implemented methods. Needs a from-source CRTM + `pyCRTM` build (already done, outside this repo -- see the
    module's own docstring "Setup" section for the exact recipe and the non-obvious gotchas it took: `-fPIC`,
    matching CRTM's I/O endianness to pyCRTM's little-endian-only assumption, matching Fortran compilers). Runs
    at `sensor_id="abi_g16"` rather than GOES-19's own `abi_g19` (the installed CRTM release predates GOES-19's
    coefficient files entirely; GOES-16/19 share the same ABI instrument design, so this is a reasonable
    spectral-response proxy -- see the module docstring). `pyCRTM` and `h5py` are installed into this project's
    `.venv`; `pyCRTM` itself isn't declared in `pyproject.toml` (it can only be built against a specific local
    CRTM install path outside this repo) so a bare `uv sync` can silently drop it -- check with
    `uv run python -c "import pyCRTM"` afterward and reinstall per the module docstring if needed. CRTM runs one
    profile per WRF grid column with no batching, so a single whole-domain call is too large for a login-node's
    memory budget on the full native WRF grid (e.g. 149 levels x 389x549 columns here, confirmed via real testing:
    killed by the OS's memory cgroup at ~27GB resident) -- `read_simulated_brightness_temperature` now runs at
    full resolution by default anyway, transparently splitting the domain into memory-bounded row-chunks and
    stitching the results back together (numerically exact, since CRTM's forward model is column-independent;
    verified bit-identical to a non-chunked call on a small case, and a real full-resolution run completed in
    ~145s across 56 chunks with zero NaN). Results are cached to disk under `crtm_cache/` (gitignored, like
    `outputs/`/`goes_data/`/`mrms_data/`) since a full-resolution run is still slow -- a cached result reloads in
    ~0.02s. The `stride` parameter still exists for an even-quicker coarse preview (spatial subsampling,
    independent of and composable with the chunking). Remaining approximations, individually marked in the code
    as unverified: the layer-interface pressure approximation, the fixed per-hydrometeor effective radii, and the
    ozone climatology placeholder (the geostationary zenith-angle formula *is* independently verified -- see its
    docstring).
  - `plotting.py` -- the comparison plots described above, plus shared helpers (`_plot_panel`, `_plot_row`,
    `_crop_to_extent`, `_domain_outline`, BNF site marker, etc.). `_plot_row` takes a flat list of axes, so it
    works for both the 3-panel row layouts and the 4-panel 2x2 grid (pass `axes.flatten()`). Cartopy `GeoAxes`
    force an equal-area aspect after `set_extent`, so a `figsize` whose aspect doesn't match the domain's actual
    projected aspect ratio "letterboxes" every panel (visible as dead whitespace between panels, not fixable via
    `constrained_layout` padding) -- the 2x2 `figsize` was tuned by computing the real aspect ratio directly
    (`proj.transform_points()` on the extent's corners, not the raw lon/lat degree span) and then verifying
    against the actual rendered PNG.
  - `_eccodes_setup.py` -- works around the pip `eccodes` package not being able to find its own bundled native
    library.
- `notebooks/` -- one notebook per WRF case day, each running several run/microphysics-scheme variants at one or
  more times.
  - `wrf_hrrr_obs_comparison_{20250502,20250520,20250917}_d1.ipynb` -- current, working comparisons using
    `plot_run_refl_comparison` (3-column reflectivity) and `plot_run_tb_comparison_4panel` (2x2 Tb: WRF `ctt` |
    WRF CRTM | HRRR | GOES). Each cell passes `wrf_tb_kwargs={"cache_dir": "../crtm_cache"}` explicitly --
    `nbconvert --execute` runs with cwd set to `notebooks/`, so the bare relative default (`"crtm_cache"`) would
    otherwise create a second, stray cache directory under `notebooks/` instead of reusing the top-level one (hit
    this for real: had to consolidate `notebooks/crtm_cache/` back into `crtm_cache/` after the first full run).
  - `wrf_hrrr_comparison_{20250502,20250520,20250917}_d1.ipynb` -- superseded 2x2 (WRF/HRRR x reflectivity/OLR)
    versions, kept for reference; these call `plot_run_vs_hrrr`, which no longer exists, so they won't run as-is.
- `outputs/`, `goes_data/`, `mrms_data/`, `crtm_cache/` -- generated PNGs / downloaded observation files / cached
  CRTM-derived WRF brightness temperature (see `crtm.py`), gitignored.
- `satoshi_forcing_data/` and `satoshi_testruns/` -- symlinks into shared project storage
  (`/gpfs/wolf2/arm/cli120/proj-shared/sey/bnf/wrf/...`), **not part of the git repo** (untracked, and large).
  - `satoshi_forcing_data/` holds reference/forcing datasets (`era5/`, `era5rda/`, `hrrr/`).
  - `satoshi_testruns/` holds WRF LASSO BNF test-run output directories, each named for a config/date
    (e.g. `20250502lassobnfwrfera5ml3`, `20250917lassobnfwrfhrrr2`). These are read-only reference data owned by
    another user (`sey`) — do not modify files under these symlinked paths.
