"""Shared helper functions for Homework 3."""

from __future__ import annotations

import csv
from pathlib import Path


def ensure_parent(path):
    """Create the parent directory for a file path if needed."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def format_float(value):
    """Format a float compactly with up to 12 significant digits."""
    return f"{float(value):.12g}"


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


def write_weekly_csv(path, weeks, columns):
    """Write per-week time series to CSV.

    'columns' maps a column name to an iterable of values aligned with 'weeks'.
    Floats are formatted via 'format_float' for compact, reproducible output.
    """
    rows = []
    for index, week in enumerate(weeks):
        row = {"week": int(week)}
        for name, values in columns.items():
            row[name] = format_float(values[index])
        rows.append(row)
    write_rows_csv(path, rows, fieldnames=["week", *columns.keys()])
