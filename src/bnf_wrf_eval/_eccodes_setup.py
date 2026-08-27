"""Make the pip-installed `eccodes` package able to find its C library.

The PyPI `eccodes` wheel ships only the Python bindings; the actual ecCodes
shared library and its dependencies (jasper, openjpeg, ...) come from the
`ecmwflibs` package instead. `findlibs` (used internally by `eccodes`)
doesn't know how to locate libraries inside `ecmwflibs`, so a plain
`import eccodes` fails with "Cannot find the ecCodes library". Call
`ensure_eccodes_loadable()` before importing `eccodes`.
"""

from __future__ import annotations

import ctypes
import functools
import os
import re
import tempfile
from pathlib import Path

# Matches the dynamic linker's "cannot open shared object file" message,
# which names only the missing dependency, e.g.
# "libjasper-1f616703.so.7.0.0: cannot open shared object file: ...".
_MISSING_LIB_RE = re.compile(r"^(lib[\w.+-]+\.so[\w.]*): cannot open shared object file")


@functools.cache
def ensure_eccodes_loadable() -> None:
    import ecmwflibs

    target = Path(ecmwflibs.find("eccodes"))
    lib_dir = target.parent

    # `ecmwflibs` bundles ecCodes' actual dependencies (jasper, jpeg,
    # png, ...) alongside unrelated libraries for other ECMWF tools --
    # including its own libnetcdf/libhdf5, built with a different ABI
    # than the `netcdf4` package's own bundled copies of the same libs.
    # Preloading everything in the directory (RTLD_GLOBAL) would load
    # both copies into one process and segfault, so instead we only
    # preload what libeccodes.so actually fails to resolve, discovered
    # by retrying and reading the missing library's name off the error.
    for _ in range(20):
        try:
            ctypes.CDLL(str(target), mode=ctypes.RTLD_GLOBAL)
            break
        except OSError as e:
            match = _MISSING_LIB_RE.match(str(e))
            if match is None:
                raise
            dep_path = lib_dir / match.group(1)
            if not dep_path.exists():
                raise
            ctypes.CDLL(str(dep_path), mode=ctypes.RTLD_GLOBAL)
    else:
        raise OSError(f"Could not resolve all dependencies of {target}")

    # findlibs looks for a file literally named "libeccodes.so" under
    # "$ECCODES_DIR/lib"; ecmwflibs ships it with a hash suffix in the
    # name, so point ECCODES_DIR at a cache directory holding a symlink.
    cache_dir = Path(tempfile.gettempdir()) / "bnf_wrf_eval_eccodes"
    lib_link_dir = cache_dir / "lib"
    lib_link_dir.mkdir(parents=True, exist_ok=True)
    link = lib_link_dir / "libeccodes.so"
    if not link.is_symlink() or link.resolve() != target.resolve():
        link.unlink(missing_ok=True)
        link.symlink_to(target)
    os.environ.setdefault("ECCODES_DIR", str(cache_dir))
