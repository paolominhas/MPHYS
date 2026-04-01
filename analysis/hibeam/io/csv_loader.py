"""
hibeam.io.csv_loader — Load auxiliary CSV and text data
========================================================
Thin wrappers around pandas/numpy for supplementary data files
such as npoints_dump.csv, energies.csv, fit_results_all_runs.txt,
and Bethe-Bloch reference tables.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def load_csv(
    filepath: str | Path,
    **pandas_kwargs: Any,
) -> "pd.DataFrame":
    """Load a CSV file into a pandas DataFrame.

    Any keyword arguments are forwarded to ``pd.read_csv()``.
    Common options: ``sep``, ``header``, ``names``, ``comment``.
    """
    import pandas as pd

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"CSV file not found: {filepath}")

    df = pd.read_csv(filepath, **pandas_kwargs)
    print(f"  Loaded {len(df):,} rows from {filepath.name}")
    return df


def load_npoints(filepath: str | Path) -> "pd.DataFrame":
    """Load npoints_dump.csv (produced by extract_npoints.C).

    Expected columns: ``file``, ``nPoints``.
    """
    return load_csv(filepath)


def load_dedx_table(filepath: str | Path) -> np.ndarray:
    """Load a dE/dx reference table (e.g. dedx_p_in_CD2.txt).

    Expects whitespace-separated columns with optional comment lines
    starting with ``#``.
    """
    return np.loadtxt(filepath, comments="#")


def load_fit_results(filepath: str | Path) -> list[dict]:
    """Parse fit_results_all_runs.txt into a list of dicts.

    The file format is assumed to be tab- or space-separated with
    a header row.
    """
    import pandas as pd

    df = pd.read_csv(filepath, sep=r"\s+")
    return df.to_dict("records")
