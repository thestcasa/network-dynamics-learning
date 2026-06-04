"""Shared helper functions for Homework 2."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def ensure_parent(path):
    """Create the parent directory for a file path if needed."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def write_single_row_csv(path, row):
    """Write one dictionary row to a CSV file."""
    ensure_parent(path)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def write_rows_csv(path, rows, fieldnames=None):
    """Write a list of dictionary rows to a CSV file."""
    ensure_parent(path)
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def confidence_interval_95(samples):
    """Return mean, sample std, standard error, and normal 95 percent CI."""
    values = np.asarray(samples, dtype=float)
    mean = float(np.mean(values))
    std = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
    stderr = float(std / np.sqrt(values.size)) if values.size > 0 else float("nan")
    half_width = 1.96 * stderr
    return mean, std, stderr, float(mean - half_width), float(mean + half_width)
