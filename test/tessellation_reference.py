"""Load the Python implementation and build a separate Fortran oracle.

Only the tessellation sources are compiled; an installed PyBDSF is not needed.
Builds live in a content-addressed temporary directory, never in the package.
"""

import hashlib
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def _load(name, path):
    """Load a Python file or extension directly, without importing PyBDSF.

    This lets the scripts run without building the unrelated fitting and
    natgrid extensions. Native-extension names must match their init symbol.
    """
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def implementations(checked=False):
    """Return (Python module, independently compiled Fortran reference module).

    Reuse a temporary cached extension when its sources, compiler version,
    flags, Python version and NumPy version match. On a cache miss, invoke
    F2PY with the current interpreter and retain compiler output in build.log.
    A failed build raises RuntimeError pointing to that log.

    ``checked=True`` enables Fortran bounds checking for diagnostic use; keep
    it false for benchmarks. Runtime bounds violations can terminate the
    process, so probes of unsafe legacy inputs belong in a separate process.
    """
    import numpy as np

    python = _load("tessellation_python", ROOT / "bdsf" / "_tessellation.py")
    sources = [ROOT / "src" / "fortran" / name for name in
               ("pytess_simple.f", "pytess_roundness.f", "constants.h")]
    flags = "-O3 -fno-fast-math" + (" -fcheck=bounds" if checked else "")
    compiler = subprocess.check_output(["gfortran", "--version"])
    # Include the constants header as well as the compiled source files: edits
    # to any reference input must invalidate a previously built oracle.
    key = hashlib.sha256(b"".join(p.read_bytes() for p in sources)
                         + compiler + flags.encode() + sys.version.encode()
                         + np.__version__.encode()).hexdigest()[:16]
    directory = Path(tempfile.gettempdir()) / ("pybdsf-tessellation-" + key)
    directory.mkdir(exist_ok=True)
    name = "_tessellation_reference"
    matches = list(directory.glob(name + "*.so")) + list(directory.glob(name + "*.pyd"))
    if not matches:
        command = [sys.executable, "-m", "numpy.f2py", "-c", "-m", name,
                   *(str(p) for p in sources[:2]), "-I" + str(sources[0].parent),
                   "--f77flags=" + flags, "--f90flags=" + flags]
        if np.lib.NumpyVersion(np.__version__) >= "2.0.0":
            # NumPy 2 uses Meson; older F2PY selects its supported backend.
            command.extend(["--backend", "meson"])
        env = os.environ.copy()
        env.update(FFLAGS=flags, FCFLAGS=flags)
        # Preserve floating-point semantics: fast-math could change NaN handling
        # and reductions, invalidating the reference used for exact parity.
        completed = subprocess.run(command, cwd=directory, env=env, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (directory / "build.log").write_text(completed.stdout)
        if completed.returncode:
            raise RuntimeError("Fortran oracle build failed; see " + str(directory / "build.log"))
        matches = list(directory.glob(name + "*.so")) + list(directory.glob(name + "*.pyd"))
    return python, _load(name, matches[0])
