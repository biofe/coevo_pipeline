"""Configuration loading and validation for the coevo pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


DEFAULT_CONFIG: dict[str, Any] = {
    "blast": {
        "protein_db": "nr",
        "nucleotide_db": "core_nt",
        "threads": 8,
        "evalue": 0.05,
        "max_target_seqs": 50000,
    },
    "filters": {
        "min_identity": 30,
        "min_alignment_coverage": 0.7,
    },
    "motif": {
        "positions": [1408, 1409, 1410],
        "residues": ["A", "C", "G"],
    },
    "alignment": {
        "method": "mafft",
    },
    "output": {
        "results_dir": "results",
    },
}


def load_config(config_path: str | Path) -> dict[str, Any]:
    """Load pipeline configuration from a YAML file.

    Parameters
    ----------
    config_path:
        Path to the YAML configuration file.

    Returns
    -------
    dict
        Merged configuration dictionary (defaults overridden by file values).
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with config_path.open() as fh:
        user_config: dict[str, Any] = yaml.safe_load(fh) or {}

    config = _deep_merge(DEFAULT_CONFIG, user_config)
    return config


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge *override* into *base* and return a new dict."""
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def get_results_dir(config: dict[str, Any]) -> Path:
    """Return the resolved results directory from config."""
    return Path(config["output"]["results_dir"])
