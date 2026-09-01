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
  layout additionally shows a second WRF-side method alongside the CRTM panel, for a direct look at how it
  compares -- `wrf_panel="ctt"` (default) uses `ctt`; `wrf_panel="olr"` uses WRF's quick OLR-derived Tb
  (`wrf.read_brightness_temperature_from_olr`, the Yang & Slingo 2001 broadband `OLR = sigma*Tf**4` fit, also
  exposed for HRRR as `hrrr.read_brightness_temperature_from_olr`). The notebooks under `notebooks/` currently
  use the `ctt` 4-panel version. Both Tb layouts share a custom
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
    temperature from the public `noaa-goes19` AWS Open Data S3 bucket (anonymous, no credentials). `find_goes_file`
    hits S3 to locate the nearest scan; `find_local_goes_file` resolves an already-downloaded file purely from
    local filenames (**no network** -- see the compute-node note below). boto3 clients carry connect/read timeouts
    so a stray S3 call on an offline host fails in seconds instead of hanging.
  - `mrms.py` -- finds/downloads/reads MRMS composite reflectivity from the public `noaa-mrms-pds` AWS Open Data
    S3 bucket. `find_local_mrms_file` is the no-network local-filename resolver, same as `goes.py`'s.
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
    profile per WRF grid column with no batching, so `read_simulated_brightness_temperature` splits the domain
    into row-chunks (`max_profiles_per_chunk`) and stitches the results back together -- numerically exact, since
    CRTM's forward model is column-independent, verified bit-identical to a non-chunked call. Keep chunking on: it
    costs almost nothing (see the 2026-08-30 calibration note below -- per-chunk `CRTM_Init` overhead turned out
    to be negligible, not the "~5% of runtime" once guessed) and it bounds memory (`chunk=8000` -> ~3.5 GB peak
    vs. ~49 GB for a whole-domain call). An earlier version of this note said a whole-domain call *segfaults*; on
    the current v2.4.1 build it does not (it runs, just slower and hungrier), so chunking is now a speed/memory
    choice, not a hard stability crutch -- but see the calibration note for a small-CRTM-call segfault that does
    still occur. Results are cached to disk under `crtm_cache/` (gitignored, like
    `outputs/`/`goes_data/`/`mrms_data/`) since a full-resolution run is still slow -- a cached result reloads in
    ~0.02s. The `stride` parameter still exists for an even-quicker coarse preview (spatial subsampling,
    independent of and composable with the chunking). `n_threads` sets CRTM's OpenMP thread count (pyCRTM
    defaults it to 1, so runs were single-threaded on a 128-core machine until this was added); measured on a
    compute node at full resolution, 213,561 profiles: 1 thread 243.6s, 8 threads 73.4s (3.3x), 32 threads 59.0s,
    64 threads 56.6s (4.3x ceiling) -- **8 is the practical sweet spot**, and every thread count was verified
    bit-identical to the 1-thread reference. See `scripts/crtm_thread_benchmark.sbatch` for the reproducible
    SLURM benchmark (SLURM is the live scheduler here; `bsub`/LSF is installed but broken). Cloud effective
    radius and land surface type are now
    microphysics-scheme-aware rather than fixed defaults (`_cloud_categories`; see the module's own STATUS
    docstring for the full picture): Thompson (`mp_physics`=8/28) replicates real UPP source formula-for-formula,
    verified against UPP's actual `CALRAD_WCLOUD_newcrtm.f`. P3 (`mp_physics`=52/53) has no UPP reference to
    match (UPP's public source has no P3 support at all): its liquid categories (cloud water, rain) replicate
    the Morrison 2-moment scheme's own gamma-distribution formulas instead (from the actual WRF source used for
    this project's runs), while its ice categories still use a simpler, explicitly-flagged-UNVERIFIED
    monodisperse approximation -- P3 itself computes a real internal ice effective radius via a genuine
    multi-dimensional lookup-table interpolation (confirmed by reading `module_mp_p3.F` directly, including the
    actual lookup table data files), porting that is a substantial undertaking, and it isn't in this project's
    current wrfout output anyway (`RE_ICE`, confirmed not gated behind any namelist flag -- just not requested
    in these runs' output list). This work also caught a real bug where P3 2-ice-category runs' second ice
    species (`QICE2`) was never being read at all. `Land_Type` now maps WRF's own vegetation category through
    UPP's real IGBP-to-NPOESS lookup table when available, confirmed against both UPP's source and the CRTM
    v2.4.0 User Guide. Also confirmed directly from UPP's source: HRRR's own CRTM Tb calculation uses *no* real
    ozone profile at all (`o3=0.0` for WRF-ARW models, "for now" per UPP's own comment) since HRRR's GRIB2 output
    has no ozone field to read -- our own placeholder climatology is already more complete than UPP's HRRR
    treatment, so there's no UPP approach to adopt there. Remaining approximations, individually marked in the
    code as unverified: the ozone climatology placeholder and the P3 ice effective-radius approximation (the
    geostationary zenith-angle formula *is* independently verified -- see its docstring).

    **Layer-interface ("level") pressure (2026-08-29): `_layer_pressures` now reproduces WRF's own `p8w`**
    rather than a geometric mean of adjacent layer pressures. `p8w` is the interface pressure WRF itself passes
    to its physics/radiation packages (`phy_prep` in `dyn_em/module_big_step_utilities_em.F`:
    `p8w(k) = fzm(k)*p_phy(k) + fzp(k)*p_phy(k-1)`, `p_phy = P + PB`): a linear-in-eta interpolation of the
    layer full pressure (`P` + `PB`) onto the staggered w-levels using WRF's own vertical-stretch weights
    `FNM`/`FNP` (= `fzm`/`fzp`; 1-D in `bottom_top`, in every wrfout, index 0 unused), with the surface
    interface set to `PSFC` and the model top to `P_TOP`. This is internally consistent with the `P + PB` layer
    values already fed to CRTM. Verified against WRF's hydrostatic mass integration (`p_hyd_w`, from
    `p_hyd_w(k) = p_hyd_w(k+1) - (1+qtot)*(c1h(k)*MUT+c2h(k))*dnw(k)`) and the stored `P_HYD`: the half-level
    average of the `p8w` field matches `P_HYD` to ~0.1 Pa (mean) and the field itself matches `p_hyd_w` to
    ~0.3 Pa (mean); the old geometric-mean approximation carried a systematic ~5-15 Pa mid-tropospheric bias.
    wrf-python has no ready-made tool for this -- its `getvar` "pres"/"pressure" is just `P + PB` on mass
    levels, and its only staggered products are `geopt_stag`/`zstag`; the RIP-heritage `dpfcalc`/`wrfcttcalc`
    Fortran routines do build an internal full-level pressure but only as a crude arithmetic mean of adjacent
    half-levels, and it isn't exposed to Python. This changed CRTM inputs, so every pre-`p8w` `crtm_cache/` entry
    was stale; they have since been regenerated wholesale -- see the "CRTM Tb pre-computation" note below. (The
    committed notebook outputs are still the old ones and would need a re-execute to pick up the new physics.)

    **CRTM version (2026-08-29): the environment now runs CRTM v2.4.1-jedi, not v2.4.0.** The build lives at
    `/gpfs/wolf2/arm/cli120/scratch/hengxiao80/crtm241` (a git worktree of the CRTM repo at tag `v2.4.1-jedi`),
    with v2.4.1 coefficients under `.../scratch/hengxiao80/crtm_fix/fix_REL-2.4.1_20221109/fix`. The v2.4.0
    install in `~/crtm` is untouched, and `~/crtm/src/pycrtm/setup.cfg.v240-backup` restores the old config, so
    rollback is: restore that file and rebuild/reinstall pyCRTM. Confirm which version is live with
    `uv run python -c "import pyCRTM; ..."` -- or just watch for the `CRTM Version: v2.4.x` banner CRTM prints
    at runtime. Two gotchas found while upgrading, both worth knowing if it is ever redone:
    (1) v2.4.1's `make.dependencies` omits three real dependencies of the new `CRTM_Active_Sensor.f90`
    (`CRTM_Atmosphere_Define`, `ODPS_CoordinateMapping`, `CRTM_GeometryInfo_Define`), which races under `make
    -j` and produces a misleading cascade of "not a field name" errors whose *first* error is really
    "Error in opening the compiled module file"; patch the dependency line rather than falling back to a serial
    build. (2) The upgrade alone changes *nothing* numerically -- verified bit-identical Tb (matching SHA-256)
    across all 213,561 pixels on two cases -- because pyCRTM hardcodes the old temperature-independent IR water
    emissivity table. `crtm.py`'s `IR_WATER_COEFF_FILE` now explicitly requests v2.4.1's temperature-dependent
    `Nalli2.IRwater.EmisCoeff.bin`, which is what actually delivers the upgrade's physics: verified to change
    *only* water pixels (184/187 changed, max 0.021K, mean -0.011K) and leave all 13,337 land pixels
    bit-identical. It falls back to the classic table with a `RuntimeWarning` if the file is missing, so a
    rollback to v2.4.0 coefficients keeps working rather than hard-failing.

    **CRTM Tb pre-computation (2026-08-30): the forward-model step is now separable from plotting.**
    `bnf_wrf_eval.crtm_precompute` (run as `python -m bnf_wrf_eval.crtm_precompute`) is a batch driver: it
    enumerates a run directory's `wrfout_<domain>_*` files (valid time parsed from the name), builds a
    deterministic flat work list, slices it `[task_index::task_count]` for a SLURM job array (strided, so each
    task gets an even mix of case days/times), and calls `read_simulated_brightness_temperature` per time --
    skipping any already in `crtm_cache/` (restartable) and retrying per-file exceptions.
    `scripts/precompute_crtm_cache.sbatch` + `scripts/precompute_crtm_runs.txt` are the 128-task array that
    drives it (16 tasks/node x 8 threads, `chunk=8000`, a task-level rerun loop for the segfault below). Used
    once (job 1023463, 8 nodes, ~47 min) to fill `crtm_cache/` with **all 2601 snapshots** = the 9 d1 HRRR3 runs
    (`2025{0502,0520,0917}lassobnfwrfhrrr3/{rund1,rund1-mp52,rund1-mp53}`) x every 15-min `wrfout_d01` time
    (289 each): `ok=2598 skipped=3 failed=0`, zero NaN, all `(389, 549)`, ~6.5 GB on disk.

    Two supporting bits in `crtm.py`:
    - `CRTM_INPUT_VERSION` (currently `"2026-08-29a"`, covering the `p8w` + `Nalli2` IR-water changes) is now
      part of every `_cache_path` filename, so a WRF->CRTM mapping change orphans stale `.npz` instead of
      silently serving them. **Bump it on any future change that alters the Tb for the same wrfout file.**
    - `require_cache=` makes `read_simulated_brightness_temperature` raise `FileNotFoundError` on a cache miss
      instead of computing inline, for callers that expect `crtm_precompute` to have run already. Not yet wired
      into `plotting.py` or the notebooks.

    Calibration (`scripts/crtm_calibration.sbatch` + `crtm_calibration_probe.py`, job 1023455 on a 256 GB
    `batch_all` node) set the sbatch's numbers and corrected two earlier assumptions:
    - Full resolution (213,561 profiles), one snapshot: 1 thread 160 s, 2 thr 102 s, 4 thr 71 s, 8 thr 55 s
      (~2.9x wall). `max_profiles_per_chunk` barely affects speed (a 130-chunk and a 14-chunk run finish within
      ~4 s), so per-chunk `CRTM_Init` overhead is negligible and there is no speed reason to prefer large chunks.
    - NOCHUNK did **not** segfault on the v2.4.1 build (ran clean at ~48.6 GB, 83 s) -- the old "whole-domain
      call segfaults" claim is stale. A rare tbbmalloc SIGSEGV (`ODPS_Predictor_Define.f90`, in CRTM's OpenMP
      region) *was* seen, but only on very small CRTM calls (<~2000 profiles); every full-resolution run
      (chunks >= 2000) was stable across ~55 processes and the full 2601-snapshot batch. The sbatch reruns a
      failed task up to 3x regardless.
    - Under job-array packing the real limiter is GPFS contention on the 3.4 GB `wrfout` reads (`_read_wrf_state`
      went from ~3 s solo to ~100 s at 16-wide), not CPU or CRTM: measured throughput ~366 snapshots/hr/node,
      and 16 tasks/node x 8 threads is the chosen packing.
  - `crtm_precompute.py` -- batch driver that fills `crtm_cache/` from `crtm.py`'s per-snapshot compute, so
    plotting only ever reloads `.npz`. See the "CRTM Tb pre-computation" note under `crtm.py` above for the
    full picture; run as `python -m bnf_wrf_eval.crtm_precompute` (no console-script entry, so `uv sync` can't
    drop it).
  - `batch_plot.py` -- frame-render driver for the "plots and movies" batch (see the "Plots and movies" note
    below). Enumerates a deterministic `(plot_type, run|case, time)` work list, slices it `[idx::count]` for a
    SLURM job array (strided), skips frames already on disk (restartable), retries per-frame exceptions. Renders
    both the hourly comparison frames (`tb4`, `tb4_simple`, `refl3`) and the six single-panel series (`wrf_tb`, `wrf_refl`,
    `hrrr_tb`, `hrrr_refl`, `goes_tb`, `mrms_refl`) into `<output-root>/frames/<type>/`. Run as
    `python -m bnf_wrf_eval.batch_plot --which comparison|single|all|<types>`; `scripts/batch_plot.sbatch` is the
    64-task array. **The observation files must be local already** (`scripts/download_obs.py`): `batch_plot`
    resolves them via `find_local_*_file` and never calls S3 on the compute node.
  - `make_movies.py` -- stitches the `batch_plot` PNGs into per-run / per-case MP4s. Groups frames by everything
    in the filename before the trailing `_<YYYYMMDD_HHMM>Z.png` stamp, sorts by that stamp, pads each frame on
    white to the group's max width/height (no rescale -- `bbox_inches="tight"` makes frame sizes vary by a few
    px), and encodes H.264 (`yuv420p`, 12 fps) via `imageio` + `imageio-ffmpeg`'s bundled static ffmpeg (there is
    no system `ffmpeg` here). Run as `python -m bnf_wrf_eval.make_movies --which comparison|single|all`.
  - `plotting.py` -- the comparison plots described above, plus shared helpers (`_plot_panel`, `_plot_row`,
    `_crop_to_extent`, `_domain_outline`, BNF site marker, etc.). `_plot_row` takes a flat list of axes, so it
    works for both the 3-panel row layouts and the 4-panel 2x2 grid (pass `axes.flatten()`). Cartopy `GeoAxes`
    force an equal-area aspect after `set_extent`, so a `figsize` whose aspect doesn't match the domain's actual
    projected aspect ratio "letterboxes" every panel (visible as dead whitespace between panels, not fixable via
    `constrained_layout` padding) -- the 2x2 `figsize` was tuned by computing the real aspect ratio directly
    (`proj.transform_points()` on the extent's corners, not the raw lon/lat degree span) and then verifying
    against the actual rendered PNG.
    Also has the **single-panel** plotters used by the movies: `plot_single_field` (core) and six wrappers --
    `plot_run_wrf_tb_single` / `plot_run_wrf_refl_single` (per run, keyed on a run directory) and
    `plot_case_{hrrr_tb,hrrr_refl,goes_tb,mrms_refl}_single` (per case, keyed on a case label + a `wrf_ref_dir`
    for the shared, time-independent domain projection/outline/extent). These share the comparison plots'
    projection / extent / domain-outline / BNF-star conventions and the fixed `DEFAULT_TB_CMAP` /
    `DEFAULT_*_LEVELS` colour scales (constant across frames, so movies don't flicker). No per-panel axes title;
    the only label is a suptitle built from the field's *actual* valid/scan time (so a GOES frame reads its 20:02
    scan, not the requested 20:00). `plot_single_field` draws once, reads the cartopy-shrunk axes box from
    `ax.get_position()`, and places the colorbar at exactly that height (a plain `constrained_layout` colorbar
    can't match the aspect-shrunk map). The obs wrappers resolve the observation file **locally first**
    (`_resolve_goes_file` / `_resolve_mrms_file` -> `find_local_*_file`), only falling back to the S3
    `find_*_file` when `auto_download_*=True`.
  - **Compute nodes here have no outbound internet.** Any boto3 S3 call (`goes.find_goes_file`,
    `mrms.find_mrms_file`, `download_*`) blocks indefinitely on TCP connect there -- this silently hung every
    task of a plotting job array until the local-filename resolvers above were added. Pre-download all
    observation files from a login node (`scripts/download_obs.py`), then run `batch_plot` without
    `--auto-download`. Cartopy Natural Earth features are already cached under `~/.local/share/cartopy` (shared
    `$HOME`), so those don't fetch.
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
- `scripts/` -- standalone helper scripts not part of the importable package:
  - `precompute_crtm_cache.sbatch` + `precompute_crtm_runs.txt` -- the 128-task SLURM job array that drives
    `bnf_wrf_eval.crtm_precompute` to fill `crtm_cache/` (see the "CRTM Tb pre-computation" note under `crtm.py`).
  - `crtm_calibration.sbatch` + `crtm_calibration_probe.py` -- the single-node calibration harness whose results
    (job 1023455) set the pre-compute sbatch's thread/chunk/packing numbers (documented under `crtm.py`).
  - `crtm_thread_benchmark.sbatch` -- the earlier, narrower SLURM benchmark of CRTM's OpenMP thread scaling;
    superseded for practical numbers by `crtm_calibration.sbatch` but kept for reference.
  - `download_obs.py` -- bulk anonymous-S3 prefetch of the GOES ABI ch.13 and MRMS composite-reflectivity files
    the plots need, for the three case windows, into `goes_data/` / `mrms_data/`. `--hourly` (for the comparison
    plots) or `--which goes,mrms --cadence-mrms 15` (GOES native ~5 min + MRMS every 15 min, for the single-panel
    movies). Idempotent. **Run this from a login node before `batch_plot.sbatch`** (compute nodes can't reach S3).
  - `batch_plot.sbatch` -- 64-task SLURM array (32 concurrent) driving `bnf_wrf_eval.batch_plot`; `WHICH=`
    (comparison | single | all | type list), `OUTPUT_ROOT=`, `EXTRA=` overrides. Restartable (skips existing
    frames). ~52 s/frame for the 4-panel `tb4` under 32-wide packing (the `ctt` diagnostic dominates), a few
    seconds for the single-panel types.
  The CRTM `.sbatch`/probe scripts must use the *main checkout's* `.venv` explicitly (`pyCRTM` is installed only
  there, and `uv run` from a git worktree silently creates a fresh empty venv); the two CRTM `.sbatch` files take
  `SRC=`/`RUNS=`/`PROBE=` environment overrides to point at a worktree's code instead.
- `prompts/` -- the user's task prompts / reference notes for each significant piece of work, kept for
  provenance (e.g. `CRTM_refinement*.md`, `first_comparision_plot.md`, `p8w.md` -- the WRF `p8w` interface-
  pressure formulation used to rewrite `crtm.py`'s `_layer_pressures`; `precalc_crtm.md` -- the ask that led to
  `crtm_precompute.py` and the calibration/pre-compute sbatch scripts; `plots_and_movies.md` -- the ask that led
  to `batch_plot.py` / `make_movies.py` / `download_obs.py` and the "Plots and movies" note below). Not used by
  any code.
- `outputs/`, `goes_data/`, `mrms_data/`, `crtm_cache/` -- generated PNGs / downloaded observation files / cached
  CRTM-derived WRF brightness temperature (see `crtm.py`), gitignored. `crtm_cache/` is currently fully
  populated (~6.5 GB) with all 2601 15-min snapshots of the 9 d1 HRRR3 runs, cache-key version `2026-08-29a` --
  see the "CRTM Tb pre-computation" note under `crtm.py`. `outputs/frames/<type>/` and `outputs/movies/` hold the
  batch-rendered frames and MP4s -- see the "Plots and movies" note below.
- `satoshi_forcing_data/` and `satoshi_testruns/` -- symlinks into shared project storage
  (`/gpfs/wolf2/arm/cli120/proj-shared/sey/bnf/wrf/...`), **not part of the git repo** (untracked, and large).
  - `satoshi_forcing_data/` holds reference/forcing datasets (`era5/`, `era5rda/`, `hrrr/`).
  - `satoshi_testruns/` holds WRF LASSO BNF test-run output directories, each named for a config/date
    (e.g. `20250502lassobnfwrfera5ml3`, `20250917lassobnfwrfhrrr2`). These are read-only reference data owned by
    another user (`sey`) — do not modify files under these symlinked paths.

## Plots and movies

Batch-rendered plot frames and the movies made from them, for the 9 d1 HRRR3 runs
(`2025{0502,0520,0917}lassobnfwrfhrrr3/{rund1,rund1-mp52,rund1-mp53}`, each a 72-h run at 15-min `wrfout_d01`
output). Driven by `bnf_wrf_eval.batch_plot` (frames) and `bnf_wrf_eval.make_movies` (MP4s); see those modules'
bullets under "Structure". Task prompt: `prompts/plots_and_movies.md`.

Two families:

- **Hourly comparison** (HRRR analyses are hourly, so these use on-the-hour WRF/GOES/MRMS): `tb4` = the 2x2
  `plot_run_tb_comparison_4panel` (WRF `ctt` | WRF CRTM | HRRR | GOES); `tb4_simple` = the same 2x2 but with the
  top-left panel replaced by WRF's quick OLR-derived Tb (`wrf.read_brightness_temperature_from_olr` -- the
  Yang & Slingo 2001 `OLR = sigma*Tf**4` fit, same conversion as PyFLEXTRKR's `ftfunctions.olr_to_tb`),
  selected via `plot_run_tb_comparison_4panel(..., wrf_panel="olr")` and written under the
  `wrf_simple_tb_crtm_hrrr_goes_tb_comparison` prefix so it never overwrites the `ctt` `tb4` frames; `refl3` =
  the 3-panel `plot_run_refl_comparison` (WRF | HRRR | MRMS). 73 frames/run each -> 9 movies each (27 total).
- **Single-panel**, one field per frame at its native cadence: `wrf_tb` (WRF CRTM Tb) and `wrf_refl` (WRF
  column-max reflectivity) at 15 min per run (289 frames -> 9 movies each); `hrrr_tb` / `hrrr_refl` hourly per
  case (73 frames -> 3 movies each); `goes_tb` at GOES CONUS 5 min per case (~865 frames -> 3 movies);
  `mrms_refl` at 15 min per case (289 frames -> 3 movies). 30 single-panel movies.

57 movies total, H.264 12 fps, in `outputs/movies/<frame-prefix>.mp4`; frames in `outputs/frames/<type>/`.

Regenerate (from a login node -- the download step needs internet, see the compute-node note under `plotting.py`):

    python scripts/download_obs.py --hourly                              # comparison obs
    python scripts/download_obs.py --which goes,mrms --cadence-mrms 15    # single-panel obs
    WHICH=comparison sbatch scripts/batch_plot.sbatch                     # then WHICH=single
    python -m bnf_wrf_eval.make_movies --which all

The one-time full run (2026-08-30): 10,416 frames, 0 failures; comparison array job 1023630 (~19 min), WRF/HRRR
single-panel array 1023698, GOES/MRMS single-panel array 1023731. The `tb4_simple` series was added later
(2026-09-01, `WHICH=tb4_simple` array job 1023849, 657 frames, 0 failures, ~3 min -- the OLR-fit panel is a
one-line array op so these render faster than `tb4`, whose `ctt` diagnostic dominates) followed by its 9
movies. Single-panel plots deliberately have **no
per-panel title**, a suptitle carrying the field's *actual* valid/scan time, a colorbar drawn to the exact
height of the (aspect-shrunk) map, and colorbar labels `Brightness Temp. (K)` / `Comp. Refl (dBZ)`; the GOES
suptitle reads `GOES-<n> (#<channel>)`.
