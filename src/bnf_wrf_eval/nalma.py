"""Finding, downloading, and reading North Alabama Lightning Mapping Array
(NALMA) VHF source data, from NASA's GHRC DAAC.

Unlike GLM (`goes.py`) and MRMS/NLDN (`mrms.py`), NALMA's raw data are
**VHF radiation sources**, not flash-clustered totals -- confirmed against a
real downloaded file: a single 10-minute granule held 399,149 individual
sources. A lightning flash typically produces many tens to hundreds of VHF
sources as its channel propagates, so NALMA source counts are not directly
comparable in magnitude to WRF's or GLM's *flash* counts -- see
`lightning.nalma_source_counts_on_wrf_grid`. Total lightning (IC + CC + CG,
like GLM): NALMA's VHF sources occur throughout the whole 3D discharge
channel, including the portion that reaches the ground, so it detects both
components without distinguishing them, same caveat as GLM (see
`goes.py`'s `GLM_PRODUCT` note).

**Requires a NASA Earthdata Login token** (unlike GLM/MRMS's anonymous
public S3 buckets): GHRC DAAC's actual file download endpoint
(`data.ghrc.earthdata.nasa.gov/ghrcw-protected/...`) returns 401 without
one, confirmed directly. The CMR granule *search* API is anonymous, so
`find_nalma_files_in_range` needs no token -- only `download_nalma_file`
does. Reads the token from a local file (default `~/.edl_token`, a bare
token string with no other formatting) rather than accepting it as a
function argument anywhere it could end up logged.
"""

from __future__ import annotations

import datetime as dt
import gzip
import re
from pathlib import Path

import numpy as np

# NASA CMR (Common Metadata Repository) granule search -- anonymous,
# no Earthdata login needed for search itself (only for the actual file
# download, see `download_nalma_file`).
CMR_GRANULES_URL = "https://cmr.earthdata.nasa.gov/search/granules.json"

# "North Alabama Lightning Mapping Array (LMA) V1" collection, GHRC DAAC.
COLLECTION_CONCEPT_ID = "C2683433889-GHRC_DAAC"

DEFAULT_TOKEN_FILE = "~/.edl_token"


def find_nalma_files_in_range(
    time_start: dt.datetime, time_end: dt.datetime
) -> list[tuple[str, str, dt.datetime, dt.datetime]]:
    """List every NALMA granule overlapping `[time_start, time_end]`, via
    an anonymous CMR granule search (confirmed working with no
    authentication -- only the actual file download needs a token).

    NALMA granules are ~10-minute chunks (`NALMA_YYMMDD_HHMMSS_SSSS.dat.gz`,
    the trailing field being the number of seconds analyzed), so a
    multi-hour comparison window will span several granules.

    Returns a list of (download_url, filename, granule_start, granule_end)
    tuples, sorted by start time. Empty if nothing overlaps the window
    (e.g. outside NALMA's operational period, or a gap in processing).
    """
    import requests

    params = {
        "collection_concept_id": COLLECTION_CONCEPT_ID,
        "temporal": f"{time_start:%Y-%m-%dT%H:%M:%SZ},{time_end:%Y-%m-%dT%H:%M:%SZ}",
        "page_size": 200,
    }
    resp = requests.get(CMR_GRANULES_URL, params=params, timeout=30)
    resp.raise_for_status()
    entries = resp.json()["feed"]["entry"]

    matches = []
    for entry in entries:
        data_link = next(
            (
                link["href"] for link in entry["links"]
                if link.get("rel", "").endswith("/data#") and link["href"].startswith("https://")
            ),
            None,
        )
        if data_link is None:
            continue
        granule_start = dt.datetime.strptime(entry["time_start"], "%Y-%m-%dT%H:%M:%S.%fZ")
        granule_end = dt.datetime.strptime(entry["time_end"], "%Y-%m-%dT%H:%M:%S.%fZ")
        matches.append((data_link, Path(data_link).name, granule_start, granule_end))

    matches.sort(key=lambda item: item[2])
    return matches


def _read_token(token_file: str | Path) -> str:
    path = Path(token_file).expanduser()
    if not path.exists():
        raise FileNotFoundError(
            f"No Earthdata token file at {path} -- see nalma.py's module docstring "
            "(generate one at https://urs.earthdata.nasa.gov/profile)"
        )
    return path.read_text().strip()


def download_nalma_file(
    url: str, dest_dir: str | Path, token_file: str | Path = DEFAULT_TOKEN_FILE
) -> Path:
    """Download one NALMA granule (as found by `find_nalma_files_in_range`)
    into `dest_dir`, keeping just its filename. Skips the download if a
    file of that name is already there.

    Requires a NASA Earthdata Login token (see the module docstring) --
    sent as an `Authorization: Bearer` header, confirmed to work against a
    real file (plain anonymous access 401s).
    """
    import requests

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / Path(url).name
    if not dest.exists():
        token = _read_token(token_file)
        resp = requests.get(
            url, headers={"Authorization": f"Bearer {token}"}, timeout=120, stream=True
        )
        resp.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                f.write(chunk)
        tmp.rename(dest)
    return dest


def read_nalma_sources(
    nalma_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read VHF source lon/lat/time from a downloaded NALMA `.dat.gz` file.

    File format (confirmed against a real granule): a text header ending in
    a `*** data ***` marker, a `Data start time: MM/DD/YY HH:MM:SS` line
    giving the file's reference day, then one row per source: `time (UT sec
    of day), lat, lon, alt(m), reduced chi^2, P(dBW), mask` -- fixed-width
    fields per the header's `Data format` line, parsed here by whitespace
    splitting instead (equivalent for this format, simpler).

    No quality filtering is applied (e.g. by reduced chi^2 or power) --
    every parsed source is returned as-is. Time-of-day values are not
    unwrapped across a UTC midnight rollover within a single file, since
    NALMA granules are only ~10 minutes long and don't cross it in
    practice.

    Returns (lon, lat, times) as plain 1D arrays; `times` are
    `datetime64[ns]`, combining the row's UT-seconds-of-day with the
    file's `Data start time` day.
    """
    opener = gzip.open if str(nalma_file).endswith(".gz") else open
    with opener(nalma_file, "rt", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    file_day = None
    data_start_idx = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        match = re.match(r"^Data start time:\s*(\d{2})/(\d{2})/(\d{2})", stripped)
        if match:
            mm, dd, yy = (int(x) for x in match.groups())
            file_day = dt.datetime(2000 + yy, mm, dd)
        if re.match(r"^\*+\s*data\s*\*+$", stripped, flags=re.IGNORECASE):
            data_start_idx = idx + 1
            break

    if data_start_idx is None or file_day is None:
        raise ValueError(f"{nalma_file}: couldn't find '*** data ***' marker/start time")

    secs, lats, lons = [], [], []
    for line in lines[data_start_idx:]:
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            secs.append(float(parts[0]))
            lats.append(float(parts[1]))
            lons.append(float(parts[2]))
        except ValueError:
            continue

    secs = np.asarray(secs)
    lats = np.asarray(lats)
    lons = np.asarray(lons)
    times = np.datetime64(file_day) + (secs * 1e9).astype("timedelta64[ns]")
    return lons, lats, times
