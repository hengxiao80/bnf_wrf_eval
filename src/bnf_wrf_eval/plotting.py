"""Side-by-side WRF vs. HRRR vs. observations comparison plots."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Colormap, LinearSegmentedColormap

from . import crtm as crtm_reader
from . import goes as goes_reader
from . import hrrr as hrrr_reader
from . import mrms as mrms_reader
from . import wrf as wrf_reader

DEFAULT_REFL_LEVELS = np.arange(5, 76, 5)
# 180-315 K to match the enhancement curve below: grayscale above the 240K
# convective-cloud-top threshold, a rainbow enhancement below it.
DEFAULT_CTT_LEVELS = np.arange(180, 316, 5)


def _tb_enhanced_colormap() -> LinearSegmentedColormap:
    """Grayscale-above/rainbow-below IR brightness-temperature colormap, of
    the kind commonly used to make cold convective cloud tops stand out
    (e.g. RAMMB/CIRA's "Rainbow IR" enhancement): black (warm) fading to
    white at 240 K, then a sharp jump into a cyan-to-magenta rainbow for
    the coldest (highest, most vigorous convective) cloud tops down to
    180 K. Matches a reference image the user supplied
    (`notebooks/tb-colorbar.png`).

    Calibrated against `DEFAULT_CTT_LEVELS` (180-315 K): matplotlib maps a
    colormap's [0, 1] domain evenly across whatever `levels` a caller
    passes to `contourf`, so the fixed physical anchor points below (240 K
    as the gray/rainbow boundary, etc.) only land where intended if
    `ctt_levels` keeps that same 180-315 K range. Passing a very different
    `ctt_levels` still works mechanically, just without that physical
    correspondence.
    """
    t_min, t_max = 180.0, 315.0

    def frac(t: float) -> float:
        return (t - t_min) / (t_max - t_min)

    stops = [
        (frac(180), (0.98, 0.90, 0.93)),  # pale pink (coldest, <180K arrow)
        (frac(190), (0.75, 0.30, 0.60)),  # magenta
        (frac(198), (0.45, 0.05, 0.10)),  # dark maroon
        (frac(206), (0.85, 0.10, 0.10)),  # red
        (frac(214), (0.95, 0.55, 0.05)),  # orange
        (frac(222), (0.95, 0.92, 0.15)),  # yellow
        (frac(230), (0.15, 0.70, 0.20)),  # green
        (frac(236), (0.10, 0.35, 0.85)),  # blue
        (frac(240) - 1e-6, (0.40, 0.90, 0.95)),  # cyan, right at the boundary
        (frac(240), (0.97, 0.97, 0.97)),  # sharp jump to near-white
        (frac(255), (0.80, 0.80, 0.80)),
        (frac(270), (0.55, 0.55, 0.55)),
        (frac(285), (0.35, 0.35, 0.35)),
        (frac(300), (0.15, 0.15, 0.15)),
        (frac(315), (0.02, 0.02, 0.02)),  # near-black (warmest)
    ]
    return LinearSegmentedColormap.from_list("tb_enhanced", stops)


DEFAULT_TB_CMAP = _tb_enhanced_colormap()

# DOE ARM Bankhead National Forest (BNF) site.
BNF_SITE_LAT = 34.342481
BNF_SITE_LON = -87.338177


def _domain_extent(lon: np.ndarray, lat: np.ndarray, pad_deg: float) -> tuple[float, float, float, float]:
    return (
        float(np.nanmin(lon)) - pad_deg,
        float(np.nanmax(lon)) + pad_deg,
        float(np.nanmin(lat)) - pad_deg,
        float(np.nanmax(lat)) + pad_deg,
    )


def _crop_to_extent(
    lon: np.ndarray, lat: np.ndarray, values: np.ndarray, extent: tuple[float, float, float, float], pad_cells: int = 30
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Slice a curvilinear grid down to the rows/cols overlapping `extent`.

    Keeps pcolormesh fast and avoids projecting the full (e.g. CONUS-wide
    HRRR) grid when only a small sub-domain is actually shown.

    `pad_cells` needs to be generous (not just 1-2 cells): pcolormesh
    mis-renders the outermost row/col of a rotated curvilinear mesh (shows
    as a spurious blank diagonal wedge) when that edge sits close to the
    visible map extent, so the crop must extend well past what's actually
    shown to push that artifact outside the visible area.
    """
    lon_min, lon_max, lat_min, lat_max = extent
    mask = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    if not mask.any():
        return lon, lat, values
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    r0, r1 = max(rows.min() - pad_cells, 0), min(rows.max() + pad_cells + 1, lon.shape[0])
    c0, c1 = max(cols.min() - pad_cells, 0), min(cols.max() + pad_cells + 1, lon.shape[1])
    sl = (slice(r0, r1), slice(c0, c1))
    return lon[sl], lat[sl], values[sl]


def _domain_outline(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Trace the outer edge of a 2D curvilinear grid as one closed loop
    (its four sides, corner to corner), for drawing a domain boundary
    line."""
    lon_b = np.concatenate([lon[0, :], lon[:, -1], lon[-1, ::-1], lon[::-1, 0]])
    lat_b = np.concatenate([lat[0, :], lat[:, -1], lat[-1, ::-1], lat[::-1, 0]])
    return lon_b, lat_b


def _find_wrf_file(run_dir: Path, domain: str, time: dt.datetime) -> Path:
    wrf_file = run_dir / f"wrfout_{domain}_{time:%Y-%m-%d_%H_%M_%S}"
    if not wrf_file.exists():
        raise FileNotFoundError(f"No wrfout file for {time} at {wrf_file}")
    return wrf_file


def _find_hrrr_file(hrrr_base_dir: str | Path, time: dt.datetime) -> Path:
    hrrr_file = Path(hrrr_base_dir) / f"{time:%Y%m%d}" / f"hrrr.t{time:%H}z.wrfnatf00.grib2"
    if not hrrr_file.exists():
        raise FileNotFoundError(f"No HRRR native-level analysis file for {time} at {hrrr_file}")
    return hrrr_file


def _resolve_mrms_file(time: dt.datetime, mrms_dir: str | Path, auto_download: bool) -> Path:
    """Locate the MRMS file for `time`: from local filenames first (no S3
    call -- works on internet-less compute nodes), falling back to the
    S3-backed lookup+download only when `auto_download` is set."""
    mrms_dir = Path(mrms_dir)
    try:
        return mrms_reader.find_local_mrms_file(time, mrms_dir)
    except FileNotFoundError:
        if not auto_download:
            raise
    bucket, key, _ = mrms_reader.find_mrms_file(time)
    return mrms_reader.download_mrms_file(bucket, key, mrms_dir)


def _resolve_goes_file(
    time: dt.datetime, goes_dir: str | Path, channel: int, satellite: str, auto_download: bool
) -> Path:
    """GOES counterpart of `_resolve_mrms_file`."""
    goes_dir = Path(goes_dir)
    try:
        return goes_reader.find_local_goes_file(
            time, goes_dir, channel=channel, satellite=satellite
        )
    except FileNotFoundError:
        if not auto_download:
            raise
    bucket, key, _ = goes_reader.find_goes_file(time, channel=channel, satellite=satellite)
    return goes_reader.download_goes_file(bucket, key, goes_dir)


def _output_path(output_base_dir: str | Path, prefix: str, run_dir: Path, time: dt.datetime) -> Path:
    output_base_dir = Path(output_base_dir)
    output_base_dir.mkdir(parents=True, exist_ok=True)
    run_label = f"{run_dir.parent.name}_{run_dir.name}"
    return output_base_dir / f"{prefix}_{run_label}_{time:%Y%m%d_%H%MZ}.png"


def _plot_panel(
    ax,
    lon: np.ndarray,
    lat: np.ndarray,
    values: np.ndarray,
    levels: np.ndarray,
    cmap: str | Colormap,
    extent: tuple[float, float, float, float],
    title: str,
    extend: str = "both",
    domain_outline: tuple[np.ndarray, np.ndarray] | None = None,
    site: tuple[float, float] | None = (BNF_SITE_LON, BNF_SITE_LAT),
):
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    # Pre-project lon/lat to the axes' own X/Y and call contourf with no
    # `transform=` kwarg, rather than letting cartopy do it: cartopy's
    # contourf reprojection drops cells (leaves visible holes) for a
    # curvilinear grid that's rotated relative to the axes projection --
    # e.g. HRRR's own Lambert grid plotted on WRF's Lambert axes -- even
    # though the same grid renders fine with pcolormesh.
    xyz = ax.projection.transform_points(ccrs.PlateCarree(), lon, lat)
    x, y = xyz[..., 0], xyz[..., 1]
    mesh = ax.contourf(x, y, values, levels=levels, cmap=cmap, extend=extend)
    ax.coastlines(resolution="50m", linewidth=0.8)
    # Pin an explicit resolution on every feature (BORDERS defaults to
    # cartopy's AdaptiveScaler, which picks 110m/50m/10m based on each
    # plot's map extent). Left un-pinned, a domain with a different extent
    # than previously plotted can make cartopy pick a resolution that
    # hasn't been downloaded yet and try to fetch it from
    # naturalearthdata.com mid-plot -- if that network call stalls, it
    # hangs the whole kernel in a way SIGINT can't reliably interrupt.
    # Pinning means every plot uses the same, already-cached files.
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.8)
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.5, edgecolor="black")
    # Explicit dict form (rather than the top_labels/right_labels booleans)
    # so cartopy doesn't fall back to its per-label nearest-edge geometry
    # guess, which for a Lambert projection can misplace some longitude
    # labels on the side panels.
    gl = ax.gridlines(
        draw_labels={"bottom": "x", "left": "y"}, linestyle="--", color="gray", alpha=0.5
    )
    gl.x_inline = False
    gl.y_inline = False
    # Keep longitude labels horizontal and pushed out below the axis,
    # rather than cartopy's default of rotating them to follow the
    # (slightly slanted, for a Lambert projection) gridline -- which
    # otherwise lands them inside the plot near the bottom edge.
    gl.rotate_labels = False
    gl.xpadding = 8
    if domain_outline is not None:
        outline_lon, outline_lat = domain_outline
        ax.plot(
            outline_lon, outline_lat, transform=ccrs.PlateCarree(),
            color="red", linewidth=0.8, zorder=10,
        )
    if site is not None:
        site_lon, site_lat = site
        ax.plot(
            site_lon, site_lat, transform=ccrs.PlateCarree(),
            marker="*", markersize=10, color="red",
            markeredgecolor="black", markeredgewidth=0.5, linestyle="none", zorder=11,
        )
    ax.set_title(title, fontsize=10)
    return mesh


def _plot_row(
    fig,
    axes_row,
    panels: list[tuple[np.ndarray, np.ndarray, np.ndarray, str]],
    levels: np.ndarray,
    cmap: str | Colormap,
    extent: tuple[float, float, float, float],
    domain_outline: tuple[np.ndarray, np.ndarray],
    extend: str,
    colorbar_label: str,
):
    """Draw one row of side-by-side comparison panels (e.g. WRF / HRRR /
    obs for the same variable, all on the same levels/cmap), sharing one
    colorbar. `panels` is a list of (lon, lat, values, title) tuples, one
    per axes in `axes_row`."""
    mesh = None
    for ax, (lon, lat, values, title) in zip(axes_row, panels):
        mesh = _plot_panel(
            ax, lon, lat, values, levels, cmap, extent, title,
            extend=extend, domain_outline=domain_outline,
        )
    fig.colorbar(mesh, ax=axes_row[:], label=colorbar_label, shrink=0.5, pad=0.02)


def plot_refl_comparison(
    wrf_file: str | Path,
    hrrr_file: str | Path,
    mrms_file: str | Path,
    refl_levels: np.ndarray = DEFAULT_REFL_LEVELS,
    refl_cmap: str = "turbo",
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (16, 5.5),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """Compare column-max radar reflectivity across WRF, HRRR, and MRMS
    (observations), at the single time each of the given files holds.

    Unlike the OLR-vs-brightness-temperature situation (see
    `plot_tb_comparison`), MRMS composite reflectivity is a direct,
    apples-to-apples match for WRF's simulated reflectivity: same physical
    quantity (dBZ), same basis as HRRR's `refc` -- see `mrms.py`.

    All three panels are drawn on the same Lambert Conformal projection
    (taken from the WRF domain) and the same map extent -- derived from
    the WRF domain's footprint, padded by `domain_pad_deg` -- so they're
    directly comparable. Because the WRF domain is much smaller than
    HRRR's/MRMS's, the padded extent extends past the WRF grid on all
    sides, so the WRF panel shows empty corners; that's expected. The WRF
    domain's actual footprint is outlined in thin red, and the BNF site is
    marked with a red star, on every panel.

    `refl_levels` sets the contourf contour intervals (bin edges) and can
    be overridden per call.

    Parameters
    ----------
    wrf_file : path to a wrfout_d0X_* file (single time per file).
    hrrr_file : path to a HRRR native-level analysis GRIB2 file
        (hrrr.tHHz.wrfnatf00.grib2).
    mrms_file : path to a downloaded MRMS composite reflectivity GRIB2.gz
        file (see `mrms.py`) for a scan time close to `wrf_file`'s.
    suptitle : if given, drawn as a big figure-level title above all three
        panels; see `plot_run_refl_comparison` for a wrapper that fills
        this in automatically.
    out_file : if given, the figure is saved there; otherwise it's left
        for the caller to show/save.

    Returns
    -------
    matplotlib.figure.Figure
    """
    wrf_lon, wrf_lat, wrf_refl, wrf_time = wrf_reader.read_column_max_reflectivity(wrf_file)
    proj = wrf_reader.get_lambert_projection(wrf_file)

    hrrr_lon, hrrr_lat, hrrr_refl, hrrr_time = hrrr_reader.read_composite_reflectivity(hrrr_file)
    mrms_lon, mrms_lat, mrms_refl, mrms_time = mrms_reader.read_composite_reflectivity(mrms_file)

    extent = _domain_extent(wrf_lon, wrf_lat, domain_pad_deg)
    hrrr_lon_c, hrrr_lat_c, hrrr_refl_c = _crop_to_extent(hrrr_lon, hrrr_lat, hrrr_refl, extent)
    mrms_lon_c, mrms_lat_c, mrms_refl_c = _crop_to_extent(mrms_lon, mrms_lat, mrms_refl, extent)
    domain_outline = _domain_outline(wrf_lon, wrf_lat)

    fig, axes = plt.subplots(
        1, 3, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _plot_row(
        fig, axes,
        [
            (wrf_lon, wrf_lat, wrf_refl, "WRF max reflectivity (dBZ)"),
            (hrrr_lon_c, hrrr_lat_c, hrrr_refl_c, "HRRR composite reflectivity (dBZ)"),
            (
                mrms_lon_c, mrms_lat_c, mrms_refl_c,
                f"MRMS composite reflectivity (dBZ)\n{mrms_time:%Y-%m-%d %H:%M} UTC",
            ),
        ],
        refl_levels, refl_cmap, extent, domain_outline, extend="max", colorbar_label="dBZ",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")

    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_run_refl_comparison(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    mrms_dir: str | Path = "mrms_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_mrms: bool = False,
    **kwargs,
):
    """Wrapper around `plot_refl_comparison` that takes a time plus a WRF
    run's name and directory, locates the matching wrfout/HRRR/MRMS files
    itself, and labels the figure with a big suptitle showing `run_name`
    and `time`.

    Parameters
    ----------
    time : the output/analysis time to plot, e.g.
        `datetime.datetime(2025, 9, 16, 1)`, or an ISO string like
        `"2025-09-16 01:00"`.
    run_name : label for the run, used in the figure's suptitle (e.g.
        `"2.5-km (D1) Thompson"`).
    run_dir : directory holding that run's `wrfout_*` files (e.g.
        `satoshi_testruns/20250917lassobnfwrfhrrr3/rund1`).
    hrrr_base_dir : directory holding HRRR native-level GRIB2 files, laid
        out as `<hrrr_base_dir>/<YYYYMMDD>/hrrr.tHHz.wrfnatf00.grib2`.
        Defaults to `satoshi_forcing_data/hrrr/hrrrnat_data` relative to
        the repo root.
    mrms_dir : directory holding downloaded MRMS composite reflectivity
        files (as saved by `mrms.download_mrms_file`). Defaults to
        `mrms_data` relative to wherever this is run from.
    domain : WRF domain to plot, e.g. "d01" (default) or "d02".
    output_base_dir : if given (and `out_file` isn't passed explicitly via
        `**kwargs`), the figure is auto-saved to
        `<output_base_dir>/wrf_hrrr_mrms_refl_comparison_<run_dir's last
        two path components>_<time>.png`, creating the directory if
        needed -- using the directory names rather than `run_name`, since
        `run_name` is meant as a free-form label and isn't filename-safe.
    auto_download_mrms : if the matching MRMS file isn't already in
        `mrms_dir`, fetch it via `mrms.download_mrms_file` instead of
        raising `FileNotFoundError`. False by default, so plotting never
        triggers a network download unless you explicitly ask for it.
    **kwargs : forwarded to `plot_refl_comparison` (e.g. `refl_levels`,
        `out_file`).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)

    mrms_file = _resolve_mrms_file(time, mrms_dir, auto_download_mrms)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(
            output_base_dir, "wrf_hrrr_mrms_refl_comparison", run_dir, time
        )

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_refl_comparison(wrf_file, hrrr_file, mrms_file, suptitle=suptitle, **kwargs)


def plot_tb_comparison(
    wrf_file: str | Path,
    hrrr_file: str | Path,
    goes_file: str | Path,
    ctt_levels: np.ndarray = DEFAULT_CTT_LEVELS,
    cmap: str | Colormap = DEFAULT_TB_CMAP,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (16, 5.5),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
    wrf_tb_stride: int = 1,
    wrf_tb_kwargs: dict | None = None,
):
    """Compare brightness temperature across WRF, HRRR, and GOES
    (observations), at the single time each of the given files holds.

    OLR and satellite brightness temperature aren't directly comparable
    (OLR is a broadband, all-wavelength flux; Tb is a narrowband radiance
    at one wavelength converted via the Planck function), so this compares
    three things that *are* the same physical quantity throughout, all via
    the same forward-model approach (CRTM):

    - WRF: `crtm.read_simulated_brightness_temperature` -- this project's
      own CRTM forward-model calculation run directly on WRF's state, the
      same *approach* UPP uses for HRRR's `SBT114` (see below), applied
      straight to our own model output rather than relying on
      `wrf.read_cloud_top_temperature`'s simpler `ctt` diagnostic (WRF's
      own purpose-built approximation, kept in `wrf.py` for reference/
      comparison but no longer used here).
    - HRRR: `hrrr.read_simulated_brightness_temperature` -- HRRR/UPP's own
      CRTM-derived simulated brightness temperature (the same forward-model
      approach NCEP uses operationally to compare models against
      satellite).
    - GOES: `goes.read_brightness_temperature` -- the actual observed
      channel 13 (10.3 micron, the clean IR window channel) brightness
      temperature.

    Together this is about as close to an apples-to-apples three-way
    comparison as is practical: all three are either observed Tb or a
    CRTM-derived simulated Tb for the same channel.

    `cmap` (default `DEFAULT_TB_CMAP`, see `_tb_enhanced_colormap`) is a
    grayscale-above/rainbow-below IR enhancement: black (warm) fading to
    white at 240 K following the conventional satellite IR display
    convention (cold renders bright, warm renders dark), then a rainbow
    for the coldest (highest, most vigorous convective) cloud tops below
    that.

    Parameters otherwise mirror `plot_refl_comparison` -- see there for
    `domain_pad_deg`/`out_file` behavior.

    Parameters
    ----------
    wrf_file : path to a wrfout_d0X_* file (single time per file).
    hrrr_file : path to a HRRR native-level analysis GRIB2 file
        (hrrr.tHHz.wrfnatf00.grib2).
    goes_file : path to a GOES ABI CMIP NetCDF file (see `goes.py`) for a
        scan time close to `wrf_file`'s.
    suptitle : if given, drawn as a big figure-level title above all three
        panels; see `plot_run_tb_comparison` for a wrapper that fills this
        in automatically.
    wrf_tb_stride : forwarded to `crtm.read_simulated_brightness_temperature`
        as `stride` (default 1: full WRF resolution, no spatial
        subsampling). CRTM runs one profile per WRF column with no
        batching, so a naive whole-domain call can be too large for a
        login-node's memory budget -- `read_simulated_brightness_temperature`
        handles that itself by chunking the domain into memory-bounded
        row-bands and stitching the results back together (numerically
        exact -- see its docstring), and caches the result to disk (see
        `wrf_tb_kwargs`'s `cache_dir`/`use_cache`/`recompute`) since a
        full-resolution run still takes a couple of minutes the first time.
        Raise `wrf_tb_stride` for a quicker, coarser WRF panel (e.g. for a
        fast preview) if that's ever preferable to full resolution.
    wrf_tb_kwargs : extra keyword arguments forwarded to
        `crtm.read_simulated_brightness_temperature` (e.g. `sensor_id`,
        `sat_lon`/`sat_height` from a real GOES file's own projection
        metadata, `coefficient_path`, `max_profiles_per_chunk`, `cache_dir`,
        `use_cache`, `recompute`). `stride` goes through `wrf_tb_stride`
        above, not here.

    Returns
    -------
    matplotlib.figure.Figure
    """
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_file, "HGT")
    proj = wrf_reader.get_lambert_projection(wrf_file)

    wrf_lon, wrf_lat, wrf_tb, wrf_time = crtm_reader.read_simulated_brightness_temperature(
        wrf_file, stride=wrf_tb_stride, **(wrf_tb_kwargs or {})
    )
    hrrr_lon, hrrr_lat, hrrr_tb, hrrr_time = hrrr_reader.read_simulated_brightness_temperature(
        hrrr_file
    )
    goes_lon, goes_lat, goes_tb, goes_time = goes_reader.read_brightness_temperature(goes_file)

    extent = _domain_extent(domain_lon, domain_lat, domain_pad_deg)
    hrrr_lon_c, hrrr_lat_c, hrrr_tb_c = _crop_to_extent(hrrr_lon, hrrr_lat, hrrr_tb, extent)
    goes_lon_c, goes_lat_c, goes_tb_c = _crop_to_extent(goes_lon, goes_lat, goes_tb, extent)
    domain_outline = _domain_outline(domain_lon, domain_lat)

    fig, axes = plt.subplots(
        1, 3, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _plot_row(
        fig, axes,
        [
            (wrf_lon, wrf_lat, wrf_tb, "WRF simulated brightness temp. (CRTM) (K)"),
            (hrrr_lon_c, hrrr_lat_c, hrrr_tb_c, "HRRR simulated brightness temp. (K)"),
            (
                goes_lon_c, goes_lat_c, goes_tb_c,
                f"GOES brightness temperature (K)\n{goes_time:%Y-%m-%d %H:%M} UTC scan",
            ),
        ],
        ctt_levels, cmap, extent, domain_outline, extend="both", colorbar_label="K",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")

    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_run_tb_comparison(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    goes_dir: str | Path = "goes_data",
    channel: int = 13,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_goes: bool = False,
    **kwargs,
):
    """Wrapper around `plot_tb_comparison` that takes a time plus a WRF
    run's name and directory, locates the matching wrfout/HRRR/GOES files
    itself, and labels the figure with a big suptitle showing `run_name`
    and `time`.

    Parameters
    ----------
    time : the output/analysis time to plot; see `plot_run_refl_comparison`
        for accepted forms.
    run_name : label for the run, used in the figure's suptitle.
    run_dir : directory holding that run's `wrfout_*` files.
    hrrr_base_dir : directory holding HRRR native-level GRIB2 files.
        Defaults to `satoshi_forcing_data/hrrr/hrrrnat_data` relative to
        the repo root.
    goes_dir : directory holding downloaded GOES CMIP files (as saved by
        `goes.download_goes_file`, e.g. named
        `OR_ABI-L2-CMIPC-M6C13_G19_s...nc`). Defaults to `goes_data`
        relative to wherever this is run from.
    channel : GOES ABI channel to compare against (13, the clean IR window
        channel, by default; see `goes.find_goes_file`).
    satellite : which GOES satellite/bucket to look in (see
        `goes.DEFAULT_SATELLITE` for the GOES-16 -> GOES-19 operational
        switch this defaults around).
    domain : WRF domain to plot, e.g. "d01" (default) or "d02".
    output_base_dir : if given (and `out_file` isn't passed explicitly via
        `**kwargs`), the figure is auto-saved to
        `<output_base_dir>/wrf_hrrr_goes_tb_comparison_<run_dir's last two
        path components>_<time>.png`, creating the directory if needed --
        same convention as `plot_run_refl_comparison`'s `output_base_dir`.
    auto_download_goes : if the matching GOES file isn't already in
        `goes_dir`, fetch it via `goes.download_goes_file` instead of
        raising `FileNotFoundError`. False by default, so plotting never
        triggers a network download unless you explicitly ask for it.
    **kwargs : forwarded to `plot_tb_comparison` (e.g. `ctt_levels`,
        `out_file`, `wrf_tb_stride`, `wrf_tb_kwargs`).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)

    goes_file = _resolve_goes_file(time, goes_dir, channel, satellite, auto_download_goes)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(
            output_base_dir, "wrf_hrrr_goes_tb_comparison", run_dir, time
        )

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_tb_comparison(wrf_file, hrrr_file, goes_file, suptitle=suptitle, **kwargs)


def plot_tb_comparison_4panel(
    wrf_file: str | Path,
    hrrr_file: str | Path,
    goes_file: str | Path,
    ctt_levels: np.ndarray = DEFAULT_CTT_LEVELS,
    cmap: str | Colormap = DEFAULT_TB_CMAP,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (10.5, 7.5),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
    wrf_tb_stride: int = 1,
    wrf_tb_kwargs: dict | None = None,
):
    """Like `plot_tb_comparison`, but with an extra panel showing WRF's
    older `wrf.read_cloud_top_temperature` (`ctt`) diagnostic alongside the
    newer CRTM-derived WRF panel, for a direct side-by-side look at how the
    two WRF-side methods compare -- laid out as a 2x2 grid: WRF `ctt` | WRF
    simulated Tb (CRTM) on top, HRRR simulated Tb | GOES observed Tb below.

    See `plot_tb_comparison` for the rationale behind the CRTM-based panels
    and for what all the other parameters do; `wrf_tb_stride`/
    `wrf_tb_kwargs` apply only to the CRTM panel (the `ctt` panel always
    runs at WRF's native resolution -- it's a lightweight wrf-python
    diagnostic, not a per-column CRTM forward-model call, so it doesn't
    need the same memory-driven subsampling).

    Returns
    -------
    matplotlib.figure.Figure
    """
    ctt_lon, ctt_lat, wrf_ctt, _ = wrf_reader.read_cloud_top_temperature(wrf_file)
    proj = wrf_reader.get_lambert_projection(wrf_file)

    wrf_lon, wrf_lat, wrf_tb, wrf_time = crtm_reader.read_simulated_brightness_temperature(
        wrf_file, stride=wrf_tb_stride, **(wrf_tb_kwargs or {})
    )
    hrrr_lon, hrrr_lat, hrrr_tb, hrrr_time = hrrr_reader.read_simulated_brightness_temperature(
        hrrr_file
    )
    goes_lon, goes_lat, goes_tb, goes_time = goes_reader.read_brightness_temperature(goes_file)

    extent = _domain_extent(ctt_lon, ctt_lat, domain_pad_deg)
    hrrr_lon_c, hrrr_lat_c, hrrr_tb_c = _crop_to_extent(hrrr_lon, hrrr_lat, hrrr_tb, extent)
    goes_lon_c, goes_lat_c, goes_tb_c = _crop_to_extent(goes_lon, goes_lat, goes_tb, extent)
    domain_outline = _domain_outline(ctt_lon, ctt_lat)

    fig, axes = plt.subplots(
        2, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _plot_row(
        fig, axes.flatten(),
        [
            (ctt_lon, ctt_lat, wrf_ctt, "WRF cloud-top temperature (ctt) (K)"),
            (wrf_lon, wrf_lat, wrf_tb, "WRF simulated brightness temp. (CRTM) (K)"),
            (hrrr_lon_c, hrrr_lat_c, hrrr_tb_c, "HRRR simulated brightness temp. (K)"),
            (
                goes_lon_c, goes_lat_c, goes_tb_c,
                f"GOES brightness temperature (K)\n{goes_time:%Y-%m-%d %H:%M} UTC scan",
            ),
        ],
        ctt_levels, cmap, extent, domain_outline, extend="both", colorbar_label="K",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")

    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_run_tb_comparison_4panel(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    goes_dir: str | Path = "goes_data",
    channel: int = 13,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_goes: bool = False,
    **kwargs,
):
    """Wrapper around `plot_tb_comparison_4panel`, otherwise identical to
    `plot_run_tb_comparison` -- see there for what every parameter does.
    Saved output filenames use the `wrf_ctt_crtm_hrrr_goes_tb_comparison`
    prefix (vs. `plot_run_tb_comparison`'s `wrf_hrrr_goes_tb_comparison`)
    so the two don't overwrite each other when both are run for the same
    run/time.

    **kwargs : forwarded to `plot_tb_comparison_4panel` (e.g. `ctt_levels`,
        `out_file`, `wrf_tb_stride`, `wrf_tb_kwargs`).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)

    goes_file = _resolve_goes_file(time, goes_dir, channel, satellite, auto_download_goes)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(
            output_base_dir, "wrf_ctt_crtm_hrrr_goes_tb_comparison", run_dir, time
        )

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_tb_comparison_4panel(wrf_file, hrrr_file, goes_file, suptitle=suptitle, **kwargs)


# ---------------------------------------------------------------------------
# Single-panel plots
#
# One field, one map, on the same Lambert projection / map extent / WRF-domain
# outline + BNF site marker as the comparison plots above, and -- crucially for
# stitching frames into a movie -- the same fixed color scale
# (`DEFAULT_CTT_LEVELS` / `DEFAULT_REFL_LEVELS`) for every frame. Used by
# `bnf_wrf_eval.batch_plot` to render whole time series, one PNG per frame,
# which `bnf_wrf_eval.make_movies` then encodes.
#
# The WRF panels (`plot_run_wrf_tb_single`, `plot_run_wrf_refl_single`) are
# per-run, keyed on a run directory like the comparison wrappers. The HRRR /
# GOES / MRMS panels don't depend on the WRF microphysics, so they're keyed on
# a *case* (a `.../<case_dir>/` holding the per-scheme run subdirs) plus a
# `wrf_ref_dir` -- any one run subdir, used only to pull the shared domain
# projection / outline / extent (all time-independent).
# ---------------------------------------------------------------------------


def _first_wrf_file(run_dir: str | Path, domain: str) -> Path:
    """First (time-sorted) `wrfout_<domain>_*` file in `run_dir` -- for
    pulling the domain's projection/outline/extent, which don't change with
    time, without needing a file at one specific time."""
    files = sorted(Path(run_dir).glob(f"wrfout_{domain}_*"))
    if not files:
        raise FileNotFoundError(f"No wrfout_{domain}_* files in {run_dir}")
    return files[0]


def _case_output_path(
    output_base_dir: str | Path, prefix: str, case_label: str, time: dt.datetime
) -> Path:
    output_base_dir = Path(output_base_dir)
    output_base_dir.mkdir(parents=True, exist_ok=True)
    return output_base_dir / f"{prefix}_{case_label}_{time:%Y%m%d_%H%MZ}.png"


def plot_single_field(
    lon: np.ndarray,
    lat: np.ndarray,
    values: np.ndarray,
    *,
    proj,
    domain_lon: np.ndarray,
    domain_lat: np.ndarray,
    levels: np.ndarray,
    cmap: str | Colormap,
    extend: str,
    colorbar_label: str,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (7.0, 6.0),
    crop: bool = False,
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """Draw one field on one map panel, using the same projection / extent /
    domain-outline / BNF-site conventions as the comparison plots.

    There is no per-panel axes title -- the only text label is `suptitle`
    (the model/obs name plus the field's actual valid/scan time, built by
    the caller). The colorbar is drawn to exactly the height of the
    rendered map: cartopy shrinks a `set_extent` axes to the data's aspect
    ratio on draw, so the colorbar axes is placed from `ax.get_position()`
    *after* a first draw rather than via `constrained_layout` (which can't
    match that shrunk box).

    `proj` / `domain_lon` / `domain_lat` come from the WRF domain (via
    `wrf.get_lambert_projection` and `wrf.read_field(..., "HGT")` or an
    equivalent). `crop=True` slices `lon/lat/values` to the map extent first
    (do this for the wide HRRR / GOES / MRMS grids; leave it False for WRF's
    own grid, which already matches).
    """
    extent = _domain_extent(domain_lon, domain_lat, domain_pad_deg)
    if crop:
        lon, lat, values = _crop_to_extent(lon, lat, values, extent)
    domain_outline = _domain_outline(domain_lon, domain_lat)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.04, 0.04, 0.84, 0.88], projection=proj)
    mesh = _plot_panel(
        ax, lon, lat, values, levels, cmap, extent, "",
        extend=extend, domain_outline=domain_outline,
    )

    # First draw so cartopy applies the map aspect and `ax` settles to the
    # box it actually occupies; then match the colorbar to that box.
    fig.canvas.draw()
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.015, pos.y0, 0.022, pos.height])
    fig.colorbar(mesh, cax=cax, label=colorbar_label)

    if suptitle is not None:
        fig.text(
            0.5, min(pos.y1 + 0.035, 0.99), suptitle,
            ha="center", va="bottom", fontsize=13, fontweight="bold",
        )
    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_run_wrf_tb_single(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    wrf_tb_kwargs: dict | None = None,
    **kwargs,
):
    """Single-panel WRF CRTM-simulated ABI ch.13 brightness temperature for
    one run/time. Reads from `crtm_cache/` -- pass
    `wrf_tb_kwargs={"cache_dir": ..., "require_cache": True}` so a cache miss
    is a fast, loud error rather than a minutes-long inline CRTM run."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)

    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_file, "HGT")
    proj = wrf_reader.get_lambert_projection(wrf_file)
    lon, lat, tb, valid_time = crtm_reader.read_simulated_brightness_temperature(
        wrf_file, **(wrf_tb_kwargs or {})
    )

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _output_path(output_base_dir, "wrf_crtm_tb", run_dir, time)
    suptitle = kwargs.pop("suptitle", f"{run_name} -- {valid_time:%Y-%m-%d %H:%M} UTC")
    return plot_single_field(
        lon, lat, tb, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=DEFAULT_CTT_LEVELS, cmap=DEFAULT_TB_CMAP, extend="both",
        colorbar_label="Brightness Temp. (K)",
        suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_run_wrf_refl_single(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    **kwargs,
):
    """Single-panel WRF column-max `REFL_10CM` reflectivity for one run/time."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)

    proj = wrf_reader.get_lambert_projection(wrf_file)
    lon, lat, refl, valid_time = wrf_reader.read_column_max_reflectivity(wrf_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _output_path(output_base_dir, "wrf_refl", run_dir, time)
    suptitle = kwargs.pop("suptitle", f"{run_name} -- {valid_time:%Y-%m-%d %H:%M} UTC")
    levels = kwargs.pop("refl_levels", DEFAULT_REFL_LEVELS)
    return plot_single_field(
        lon, lat, refl, proj=proj, domain_lon=lon, domain_lat=lat,
        levels=levels, cmap=kwargs.pop("refl_cmap", "turbo"), extend="max",
        colorbar_label="Comp. Refl (dBZ)",
        suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_case_hrrr_tb_single(
    time: dt.datetime | str,
    case_label: str,
    wrf_ref_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    **kwargs,
):
    """Single-panel HRRR CRTM-simulated brightness temperature (`SBT114`) for
    one case/time, cropped to the WRF domain's map extent."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    wrf_ref_file = _first_wrf_file(wrf_ref_dir, domain)
    proj = wrf_reader.get_lambert_projection(wrf_ref_file)
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_ref_file, "HGT")

    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    lon, lat, tb, valid_time = hrrr_reader.read_simulated_brightness_temperature(hrrr_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _case_output_path(output_base_dir, "hrrr_tb", case_label, time)
    suptitle = kwargs.pop("suptitle", f"HRRR -- {valid_time:%Y-%m-%d %H:%M} UTC")
    return plot_single_field(
        lon, lat, tb, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=DEFAULT_CTT_LEVELS, cmap=DEFAULT_TB_CMAP, extend="both",
        colorbar_label="Brightness Temp. (K)",
        crop=True, suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_case_hrrr_refl_single(
    time: dt.datetime | str,
    case_label: str,
    wrf_ref_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    **kwargs,
):
    """Single-panel HRRR composite reflectivity (`refc`) for one case/time."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    wrf_ref_file = _first_wrf_file(wrf_ref_dir, domain)
    proj = wrf_reader.get_lambert_projection(wrf_ref_file)
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_ref_file, "HGT")

    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    lon, lat, refl, valid_time = hrrr_reader.read_composite_reflectivity(hrrr_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _case_output_path(output_base_dir, "hrrr_refl", case_label, time)
    suptitle = kwargs.pop("suptitle", f"HRRR -- {valid_time:%Y-%m-%d %H:%M} UTC")
    levels = kwargs.pop("refl_levels", DEFAULT_REFL_LEVELS)
    return plot_single_field(
        lon, lat, refl, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=levels, cmap=kwargs.pop("refl_cmap", "turbo"), extend="max",
        colorbar_label="Comp. Refl (dBZ)",
        crop=True, suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_case_goes_tb_single(
    time: dt.datetime | str,
    case_label: str,
    wrf_ref_dir: str | Path,
    goes_dir: str | Path = "goes_data",
    channel: int = 13,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_goes: bool = False,
    **kwargs,
):
    """Single-panel GOES ABI observed brightness temperature for the scan
    closest to `time`. The saved filename carries the *requested* `time`
    (a regular grid), not the scan time, so a movie's frames stay evenly
    spaced even where two grid steps snap to the same scan."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    wrf_ref_file = _first_wrf_file(wrf_ref_dir, domain)
    proj = wrf_reader.get_lambert_projection(wrf_ref_file)
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_ref_file, "HGT")

    goes_file = _resolve_goes_file(time, goes_dir, channel, satellite, auto_download_goes)
    lon, lat, tb, scan_time = goes_reader.read_brightness_temperature(goes_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _case_output_path(output_base_dir, "goes_tb", case_label, time)
    sat_num = satellite.lstrip("Gg")
    suptitle = kwargs.pop(
        "suptitle", f"GOES-{sat_num} (#{channel}) -- {scan_time:%Y-%m-%d %H:%M} UTC"
    )
    return plot_single_field(
        lon, lat, tb, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=DEFAULT_CTT_LEVELS, cmap=DEFAULT_TB_CMAP, extend="both",
        colorbar_label="Brightness Temp. (K)",
        crop=True, suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_case_mrms_refl_single(
    time: dt.datetime | str,
    case_label: str,
    wrf_ref_dir: str | Path,
    mrms_dir: str | Path = "mrms_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_mrms: bool = False,
    **kwargs,
):
    """Single-panel MRMS composite reflectivity for the scan closest to
    `time`; filename carries the requested `time` (see
    `plot_case_goes_tb_single`)."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    wrf_ref_file = _first_wrf_file(wrf_ref_dir, domain)
    proj = wrf_reader.get_lambert_projection(wrf_ref_file)
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_ref_file, "HGT")

    mrms_file = _resolve_mrms_file(time, mrms_dir, auto_download_mrms)
    lon, lat, refl, scan_time = mrms_reader.read_composite_reflectivity(mrms_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _case_output_path(output_base_dir, "mrms_refl", case_label, time)
    suptitle = kwargs.pop("suptitle", f"MRMS -- {scan_time:%Y-%m-%d %H:%M} UTC")
    levels = kwargs.pop("refl_levels", DEFAULT_REFL_LEVELS)
    return plot_single_field(
        lon, lat, refl, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=levels, cmap=kwargs.pop("refl_cmap", "turbo"), extend="max",
        colorbar_label="Comp. Refl (dBZ)",
        crop=True, suptitle=suptitle, out_file=out_file, **kwargs,
    )
