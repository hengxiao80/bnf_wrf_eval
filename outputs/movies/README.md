# BNF WRF evaluation movies

Animated comparisons of the WRF LASSO BNF test runs (Satoshi's `hrrr3` d1 (2.5-km) runs) against the HRRR
analyses that drove them and against independent satellite / radar observations,
for the DOE ARM Bankhead National Forest (BNF) site in northern Alabama.

All 48 movies are H.264 (`yuv420p`), **12 fps, no audio**. At 12 fps an hourly
(73-frame) movie runs ~6 s, a 15-min (289-frame) movie ~24 s, and the 5-min GOES
movies (~865 frames) ~72 s.

## The runs

Three case periods, each a 72-hour WRF forecast on a single 2.5-km domain
(550 x 390 grid, Lambert Conformal), initialized from HRRR at 00 UTC and writing
output every 15 minutes:

| Case label | Simulated period (UTC) |
|------------|------------------------|
| `20250502` | 2025-05-01 00Z – 2025-05-04 00Z |
| `20250520` | 2025-05-19 00Z – 2025-05-22 00Z |
| `20250917` | 2025-09-16 00Z – 2025-09-19 00Z |

Each case was run with three microphysics schemes (everything else identical):

| Run tag       | `mp_physics` | Scheme |
|---------------|--------------|--------|
| `rund1`       | 28 | Thompson aerosol-aware |
| `rund1-mp52`  | 52 | P3 (1 ice category) |
| `rund1-mp53`  | 53 | P3 (2 ice categories) |

## Reference data

- **GOES**: GOES-19 ABI channel 13 (10.3 µm "clean window" IR) brightness
  temperature, CONUS scans (~5 min).
- **HRRR**: hourly HRRR analyses — `refc` composite reflectivity and `SBT114`
  (UPP CRTM–simulated ABI ch.13 brightness temperature).
- **MRMS**: `MergedReflectivityQCComposite` composite reflectivity (sampled every
  15 min here).

## What the WRF quantities are

- **WRF CRTM Tb**: ABI ch.13 brightness temperature computed from WRF's own state
  with the CRTM forward model (the same approach UPP uses for HRRR's `SBT114`), so
  the WRF/HRRR/GOES Tb comparison is genuinely like-for-like (same forward model,
  same channel).
- **WRF `ctt`**: WRF's simpler built-in cloud-top-temperature diagnostic
  (`wrf-python`), shown alongside the CRTM panel for reference.
- **WRF reflectivity**: column maximum of the simulated `REFL_10CM` field (dBZ),
  directly comparable to the HRRR and MRMS composite reflectivity.

## Common plotting conventions

Every panel is reprojected onto the **same** Lambert Conformal projection and map
extent (taken from the WRF domain), with the WRF domain footprint outlined in
**red** and the BNF site marked with a **red star**. Colour scales are fixed
across all frames so the animations do not flicker:

- Brightness temperature: 180–315 K, 5-K steps, conventional IR enhancement
  (grayscale for 240–315 K, then a cyan-to-magenta rainbow for the coldest
  convective cloud tops down to 180 K).
- Composite / column-max reflectivity: 5–75 dBZ, 5-dBZ steps.

## Movie inventory

### Comparison movies (hourly, 73 frames each) — 18 files

| File pattern | Layout |
|--------------|--------|
| `wrf_ctt_crtm_hrrr_goes_tb_comparison_<case>_<run>.mp4` | 2×2 brightness temperature: WRF `ctt` \| WRF CRTM \| HRRR CRTM \| GOES observed |
| `wrf_hrrr_mrms_refl_comparison_<case>_<run>.mp4` | 3-panel reflectivity: WRF \| HRRR \| MRMS |

9 of each (3 cases × 3 runs).

### Single-panel movies — 30 files

One field per frame, at its native cadence.

| File pattern | Field | Cadence / frames | Count |
|--------------|-------|------------------|-------|
| `wrf_crtm_tb_<case>_<run>.mp4` | WRF CRTM ABI ch.13 Tb | 15 min / 289 | 9 |
| `wrf_refl_<case>_<run>.mp4` | WRF column-max reflectivity | 15 min / 289 | 9 |
| `hrrr_tb_<case>.mp4` | HRRR CRTM Tb (`SBT114`) | hourly / 73 | 3 |
| `hrrr_refl_<case>.mp4` | HRRR composite reflectivity (`refc`) | hourly / 73 | 3 |
| `goes_tb_<case>.mp4` | GOES-19 ABI ch.13 Tb | ~5 min / ~865 | 3 |
| `mrms_refl_<case>.mp4` | MRMS composite reflectivity | 15 min / 289 | 3 |

The observation-only single-panel movies (`hrrr_*`, `goes_tb`, `mrms_refl`) depend
only on the case, not on the microphysics run, so their filenames have no run tag.
Their suptitle carries the field's *actual* valid/scan time (e.g. a GOES frame
labelled 20:02 is the scan nearest the requested 20:00).

## Filename key

    <type>_<case>lassobnfwrfhrrr3[_<run>].mp4

- `<case>` — `20250502`, `20250520`, or `20250917` (see table above)
- `<run>` — `rund1`, `rund1-mp52`, or `rund1-mp53` (comparison and WRF
  single-panel movies only)
