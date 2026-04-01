"""
hibeam.config — Central configuration loader
=============================================
Reads config.yaml and returns a dot-accessible namespace so you can write
``cfg.paths.simulation.krakow`` instead of ``cfg["paths"]["simulation"]["krakow"]``.

Usage
-----
    from hibeam import config
    cfg = config.load()                    # auto-finds config.yaml
    cfg = config.load("my_config.yaml")    # explicit path
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml


# ── Dot-accessible dict ──────────────────────────────────────────────────────

class DotDict(dict):
    """Dictionary subclass that supports attribute-style access.

    Nested dicts are converted recursively so ``d.a.b.c`` works at any depth.
    Falls back to regular ``d["key"]`` access for keys containing special
    characters or reserved words.
    """

    def __getattr__(self, key: str) -> Any:
        try:
            val = self[key]
        except KeyError:
            raise AttributeError(
                f"Config has no key '{key}'.  "
                f"Available: {list(self.keys())}"
            ) from None
        return val

    def __setattr__(self, key: str, value: Any) -> None:
        self[key] = value

    def __delattr__(self, key: str) -> None:
        del self[key]

    @classmethod
    def from_dict(cls, d: dict) -> "DotDict":
        """Recursively convert a plain dict to DotDict."""
        out = cls()
        for k, v in d.items():
            if isinstance(v, dict):
                out[k] = cls.from_dict(v)
            elif isinstance(v, list):
                out[k] = [cls.from_dict(i) if isinstance(i, dict) else i for i in v]
            else:
                out[k] = v
        return out


# ── TPC geometry dataclass (mirrors the YAML section) ────────────────────────

@dataclass
class TPCGeometry:
    """Physical TPC parameters from config.yaml → ``geometry:`` section.

    All values must match those hard-coded in TPCHit.hh.
    If any field is None, track-length normalisation is skipped and
    the analysis falls back to sum(Edep) [MeV].
    """
    pad_width:  Optional[float] = None   # mm
    pad_height: Optional[float] = None   # mm
    drift_v:    Optional[float] = None   # mm/ns
    time_bin:   Optional[float] = None   # ns

    @property
    def is_complete(self) -> bool:
        return all(v is not None for v in
                   [self.pad_width, self.pad_height, self.drift_v, self.time_bin])

    @classmethod
    def from_config(cls, cfg: DotDict) -> "TPCGeometry":
        g = cfg.get("geometry", {})
        return cls(
            pad_width=g.get("pad_width"),
            pad_height=g.get("pad_height"),
            drift_v=g.get("drift_v"),
            time_bin=g.get("time_bin"),
        )


# ── Loader ────────────────────────────────────────────────────────────────────

_DEFAULT_NAMES = ["config.yaml", "config.yml"]


def _find_config() -> Path:
    """Walk upward from cwd looking for config.yaml in the repo root."""
    cwd = Path.cwd()
    for parent in [cwd, *cwd.parents]:
        for name in _DEFAULT_NAMES:
            candidate = parent / name
            if candidate.is_file():
                return candidate
    raise FileNotFoundError(
        "Cannot find config.yaml.  Either pass the path explicitly to "
        "config.load() or run from within the hibeam-tpc-analysis directory."
    )


def load(path: str | Path | None = None) -> DotDict:
    """Load and return the analysis configuration.

    Parameters
    ----------
    path : str or Path, optional
        Explicit path to a YAML config file.  If omitted, searches upward
        from the current directory for ``config.yaml``.

    Returns
    -------
    DotDict
        Dot-accessible configuration namespace.
    """
    if path is None:
        path = _find_config()
    else:
        path = Path(path)

    with open(path) as f:
        raw = yaml.safe_load(f)

    cfg = DotDict.from_dict(raw)

    # Attach the config file location so paths can be resolved relative to it.
    cfg["_config_dir"] = str(path.parent.resolve())

    return cfg


def resolve_path(cfg: DotDict, relative: str) -> Path:
    """Resolve a relative path from config.yaml against the repo root."""
    base = Path(cfg.get("_config_dir", "."))
    p = base / relative
    return p.resolve()
