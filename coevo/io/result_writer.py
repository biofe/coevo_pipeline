"""Result writer utilities: persist DataFrames and dicts to disk."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger


def write_dataframe(df: pd.DataFrame, file_path: str | Path, sep: str = "\t") -> None:
    """Write a pandas DataFrame to a delimited text file.

    Parameters
    ----------
    df:
        DataFrame to write.
    file_path:
        Destination file path.
    sep:
        Column separator (default: tab).
    """
    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(str(file_path), sep=sep, index=False)
    logger.info(f"Wrote DataFrame ({len(df)} rows) to {file_path}")


def write_dict(data: dict[str, Any], file_path: str | Path, sep: str = "\t") -> None:
    """Write a flat dictionary as a two-column (key, value) TSV file.

    Parameters
    ----------
    data:
        Dictionary to serialize.
    file_path:
        Destination file path.
    sep:
        Column separator (default: tab).
    """
    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(list(data.items()), columns=["key", "value"])
    df.to_csv(str(file_path), sep=sep, index=False)
    logger.info(f"Wrote dict ({len(data)} entries) to {file_path}")


def read_dataframe(file_path: str | Path, sep: str = "\t") -> pd.DataFrame:
    """Read a delimited text file into a pandas DataFrame.

    Parameters
    ----------
    file_path:
        Source file path.
    sep:
        Column separator (default: tab).

    Returns
    -------
    pandas.DataFrame

    Raises
    ------
    FileNotFoundError
        If *file_path* does not exist.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Result file not found: {file_path}")
    df = pd.read_csv(str(file_path), sep=sep)
    logger.info(f"Read DataFrame ({len(df)} rows) from {file_path}")
    return df
