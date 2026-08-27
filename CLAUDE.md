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

Two three-column comparisons, each laid out WRF | HRRR | independent observation, all reprojected onto the same
Lambert Conformal projection (taken from the WRF domain) and the same map extent, with the WRF domain's footprint
outlined in red and the BNF site marked with a red star on every panel:

- **Reflectivity** (`plotting.plot_refl_comparison` / `plot_run_refl_comparison`): WRF's column-max simulated
  `REFL_10CM` vs. HRRR's composite reflectivity (`refc`) vs. MRMS's composite reflectivity
  (`MergedReflectivityQCComposite`). All three are the same physical quantity (dBZ), so this is a direct,
  apples-to-apples comparison.
- **Brightness/cloud-top temperature** (`plotting.plot_tb_comparison` / `plot_run_tb_comparison`): WRF's `ctt`
  diagnostic (via `wrf-python`) vs. HRRR's CRTM-derived simulated brightness temperature (`SBT114`) vs. GOES ABI
  channel 13 (10.3 micron) observed brightness temperature. WRF only outputs OLR (a broadband, all-wavelength
  flux), which is *not* directly comparable to satellite brightness temperature (a narrowband radiance at one
  wavelength) — `ctt` is the model-side diagnostic purpose-built to approximate what a clean IR window channel
  would see, making this about as close to apples-to-apples as is practical without running a full
  radiative-transfer forward model (CRTM/RTTOV) on the WRF profiles directly.

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
executing notebooks.

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
  - `plotting.py` -- the comparison plots described above, plus shared helpers (`_plot_panel`, `_crop_to_extent`,
    `_domain_outline`, BNF site marker, etc.).
  - `_eccodes_setup.py` -- works around the pip `eccodes` package not being able to find its own bundled native
    library.
- `notebooks/` -- one notebook per WRF case day, each running several run/microphysics-scheme variants at one or
  more times.
  - `wrf_hrrr_obs_comparison_{20250502,20250520,20250917}_d1.ipynb` -- current, working three-column
    (WRF/HRRR/obs) comparisons using `plot_run_refl_comparison`/`plot_run_tb_comparison`.
  - `wrf_hrrr_comparison_{20250502,20250520,20250917}_d1.ipynb` -- superseded 2x2 (WRF/HRRR x reflectivity/OLR)
    versions, kept for reference; these call `plot_run_vs_hrrr`, which no longer exists, so they won't run as-is.
- `outputs/`, `goes_data/`, `mrms_data/` -- generated PNGs / downloaded observation files, gitignored.
- `satoshi_forcing_data/` and `satoshi_testruns/` -- symlinks into shared project storage
  (`/gpfs/wolf2/arm/cli120/proj-shared/sey/bnf/wrf/...`), **not part of the git repo** (untracked, and large).
  - `satoshi_forcing_data/` holds reference/forcing datasets (`era5/`, `era5rda/`, `hrrr/`).
  - `satoshi_testruns/` holds WRF LASSO BNF test-run output directories, each named for a config/date
    (e.g. `20250502lassobnfwrfera5ml3`, `20250917lassobnfwrfhrrr2`). These are read-only reference data owned by
    another user (`sey`) — do not modify files under these symlinked paths.
