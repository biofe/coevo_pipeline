"""Structured logging utilities using loguru."""

from __future__ import annotations

import sys
from pathlib import Path

from loguru import logger


def setup_logging(
    log_level: str = "INFO",
    log_file: str | Path | None = None,
) -> None:
    """Configure loguru logging for the pipeline.

    Parameters
    ----------
    log_level:
        Minimum log level (DEBUG, INFO, WARNING, ERROR).
    log_file:
        Optional path to write log output to a file.
    """
    logger.remove()
    logger.add(
        sys.stderr,
        level=log_level.upper(),
        format=(
            "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
            "<level>{level: <8}</level> | "
            "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
            "<level>{message}</level>"
        ),
        colorize=True,
    )
    if log_file is not None:
        logger.add(
            str(log_file),
            level=log_level.upper(),
            rotation="10 MB",
            retention="30 days",
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} - {message}",
        )


def get_logger(name: str) -> "logger":  # type: ignore[valid-type]
    """Return a bound loguru logger for a named component.

    Parameters
    ----------
    name:
        Component or module name to bind to log records.
    """
    return logger.bind(component=name)
