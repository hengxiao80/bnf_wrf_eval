# HRRR lightning threat field and its GLM verification — reference notes

Notes from a literature dig (2026-09-03) into what HRRR's `ltng` GRIB2 field actually is, and how
Blaylock & Horel (2020) verified it against GLM observations. Written for reference while planning
this project's own WRF-vs-observation lightning evaluation (see `prompts/lightning_evaluation.md`).

## What's in the downloaded HRRR data

Confirmed directly from the GRIB2 files under `satoshi_forcing_data/hrrr/` (both `hrrr_data/` and
`hrrrnat_data/`, all of `wrfsfc`/`wrfprs`/`wrfnat`, every forecast hour checked): there is exactly
one lightning-related field, `shortName='ltng'` / `name='Lightning'` (GRIB2 discipline 0, category 17
"Electrodynamics", parameter 192; NCEP paramId 260391), `typeOfLevel='atmosphere'` (single level),
GRIB2 units tagged `dimensionless`. No flash-rate or flash-frequency field exists in this dataset —
`ltng` is HRRR's parameterized **lightning threat index**, not an observed or model-diagnosed flash
count. It is not currently read by any code in this repo (`src/bnf_wrf_eval/`).

## The underlying algorithm: McCaul et al. (2009)

HRRR's `ltng` field implements the blended lightning-threat scheme of McCaul, Goodman, LaCasse &
Cecil (2009), *"Forecasting Lightning Threat Using Cloud-Resolving Model Simulations"*
(*Wea. Forecasting*, 24, 709–729, [doi:10.1175/2008WAF2222152.1](https://doi.org/10.1175/2008WAF2222152.1);
free pre-print at [NASA NTRS](https://ntrs.nasa.gov/citations/20080048096)). Two proxy fields, each
linearly calibrated against peak observed flash-origin density from the North Alabama Lightning
Mapping Array (LMA), are blended:

**F1 — graupel flux at −15 °C** (captures temporal variability well; areal coverage too small alone):

```
F1 = k1 * (w * q_g)_m,   k1 = 0.042
```

`w` = vertical velocity (m s⁻¹), `q_g` = graupel mixing ratio, evaluated at the −15 °C level in the
mixed-phase region ("graupel flux" in name only — no air density term, so strictly a proxy, not a
true mass flux). Calibration correlation vs. observed peak flash density: r = 0.67.

**F2 — vertically integrated ice** (captures areal coverage well; temporal variability muted):

```
F2 = k2 * ∫ ρ (q_g + q_s + q_i) dz,   k2 = 0.20
```

`ρ` = air density, `q_g`/`q_s`/`q_i` = graupel/snow/cloud-ice mixing ratios, integrated over the full
storm depth (units kg m⁻²). Calibration correlation: r = 0.83.

**F3 — blended threat** (this is what HRRR outputs as `ltng`):

```
F3 = r1*F1 + r2*F2,   r1 = 0.95, r2 = 0.05
```

Weighted 95% toward the temporally-sharp graupel-flux term, with just enough of the ice-integral term
to bring areal coverage closer to observed. Both F1 and F2 (and hence F3) are calibrated in units of
**flashes (5 min)⁻¹ per model grid-box column** — a real flash-origin-density rate, not an arbitrary
0–1 index; the GRIB2 `dimensionless` tag is just NCEP's encoding choice, since no standard GRIB2 unit
exists for this locally-defined quantity.

**Caveat from the source paper**: these specific constants (k1, k2, r1, r2) were fit from only 7 North
Alabama storm cases on a 2‑km WRF grid with WSM6 microphysics; the authors explicitly say recalibration
is needed for other model configurations/resolutions/microphysics. NOAA/GSL's operational HRRR (3‑km
grid, Thompson microphysics) likely uses its own re-tuned constants, not necessarily these exact
numbers. UPP's post-processor (`INITPOST_NETCDF.f`, per a 2023-03-02 Sam Trahan commit) reads all three
components separately as `ltg1_max`/`ltg2_max`/`ltg3_max`, but only one `ltng` message (presumably F3)
appears in the public GRIB2 output we downloaded.

**HRRR's temporal packaging**: per Blaylock & Horel (below), the hourly `ltng` value HRRR writes out is
the **maximum over the preceding hour** of the instantaneous F3 field (itself nominally a 5-min-rate
snapshot at each model step) — i.e. an hourly max-of-5-min-rate field, not an hourly mean or total.

## Blaylock & Horel (2020): verifying HRRR's lightning threat against GLM

Blaylock, B. K., and J. D. Horel, 2020: *"Comparison of Lightning Forecasts from the High-Resolution
Rapid Refresh Model to Geostationary Lightning Mapper Observations,"* *Wea. Forecasting*, 35, 401–416,
[doi:10.1175/WAF-D-19-0141.1](https://doi.org/10.1175/WAF-D-19-0141.1) — HRRRv3, May–Sept 2018 & 2019,
GOES-East GLM.

### HRRR side: binary occurrence field

The hourly `ltng` (F3) forecast is thresholded at **> 0.04 flashes km⁻² (5 min)⁻¹** to produce a binary
"lightning forecast / no forecast" field on the HRRR grid. (Lower thresholds, including >0.0, were
tested but produced too many spurious weak-convection detections; results are not sensitive to the
exact choice near 0.04.)

### GLM side: binary occurrence field, not raw flash counts

GLM Level 2 *event* data (not flashes) are pulled from the public `noaa-goes` AWS bucket and aggregated
over the same hour the HRRR forecast is valid for (180 twenty-second files per hour), using the
parallax-corrected geolocation. Rather than comparing at the single grid point nearest each event, each
GLM event is given a **footprint of the 21 HRRR grid points forming a 5×5 square around the event
location, with the 4 corner points excluded** — approximating a GLM pixel's actual ~8–14 km CCD
resolution. To suppress spurious single-pixel detections, a HRRR grid point only counts as "observed
lightning" if **at least 2 GLM events** deposited a footprint on it within that hour.

This produces an hourly binary GLM occurrence field on the HRRR grid, directly comparable to HRRR's
own thresholded binary field above.

### Skill metric: neighborhood-based, not point-by-point

Point-by-point contingency-table stats (POD, FAR) are reported for illustration but are explicitly
called out as unsuitable for convection-scale features (double-penalty problem from small
displacement errors). The primary metric is the **Fractions Skill Score (FSS)**, computed within
circular neighborhoods of radius 30/60/120/240 km around each grid point, comparing the fraction of
"lightning occurred" grid points within the neighborhood between forecast and GLM-observed binary
fields. FSS = 1 is a perfect areal match; FSS = 0 is no correspondence. A forecast is called "useful"
when FSS ≥ FSS_uniform = (1 + f0)/2, where f0 is the observed areal occurrence fraction in that
neighborhood/time — this baseline itself varies by domain and hour. Headline results: HRRR lightning
placement is not skillful at 30-km scales beyond very short lead times, but is meaningfully skillful
at 60–120-km scales for several hours, with skill increasing toward shorter lead times and during the
climatological peak convective hours (late afternoon/evening).

## Why this matters for this project

If/when this repo adds a WRF-vs-observation lightning comparison (companion to the existing
reflectivity/Tb comparisons in `plotting.py`), the natural observation targets are GLM (already used
in `bill_lightning_code/`, gridded via the University of Wisconsin `cspp-geo-gridded-glm` tool) and/or
NALMA. The relevant HRRR forcing-data field for a WRF-vs-HRRR comparison is `ltng` (F3, per above) —
but note WRF itself has no equivalent built-in diagnostic: `phys/module_lightning_driver.F` in the
public WRF repo only implements Price & Rind (1992)-style parameterizations (`ltng_crmpr92w/z`,
`ltng_cpmpr92z`) and an LPI (Lightning Potential Index) option, none of which match McCaul et al.'s
graupel-flux + ice-integral scheme — so a WRF-side F1/F2/F3 analog would need to be computed
independently from wrfout hydrometeor/vertical-velocity fields, following the same formulas above
(with fresh calibration constants, per the paper's own caution against reusing theirs across model
configurations).
