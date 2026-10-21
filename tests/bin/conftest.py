"""Shared helpers for the tests of the Python scripts in bin/.

Each script is run as a subprocess with the interpreter running pytest, the way
the pipeline runs it: the script's own directory is then first on sys.path, so
its sibling imports (sahmi_stats, filter_abundance, ...) resolve exactly as they
do in production.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / "bin"


def _run(script, *args, cwd=None, timeout=120):
    """Run bin/<script> with `args`; return the CompletedProcess (never raises on exit code)."""
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1")
    return subprocess.run(
        [sys.executable, str(BIN / script), *[str(arg) for arg in args]],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


@pytest.fixture
def run_script():
    """`run_script("name.py", *args, cwd=...)` -> subprocess.CompletedProcess."""
    return _run
