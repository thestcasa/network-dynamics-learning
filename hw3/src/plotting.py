"""Plotting helpers for Homework 3."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # non-interactive backend, safe for headless batch runs
import matplotlib.pyplot as plt  # noqa: E402


def save_figure(fig, path):
    """Save a Matplotlib figure with consistent layout handling."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_lines(weeks, series, xlabel, ylabel, title, path, bands=None, markers=True):
    """Plot one or more per-week curves on shared axes and save to 'path'.

    'series' maps a legend label to a mean array aligned with 'weeks'. 'bands'
    optionally maps the same labels to standard-deviation arrays drawn as a light
    +/-1 std shaded region around the corresponding mean.
    """
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for label, values in series.items():
        (line,) = ax.plot(
            weeks,
            values,
            marker="o" if markers else None,
            markersize=4,
            linewidth=1.8,
            label=label,
        )
        if bands is not None and label in bands:
            std = bands[label]
            ax.fill_between(
                weeks,
                values - std,
                values + std,
                color=line.get_color(),
                alpha=0.15,
            )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    save_figure(fig, path)


def plot_histogram(values, xlabel, ylabel, title, path, bins=None):
    """Plot a histogram of 'values' and save to 'path'."""
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.hist(
        values,
        bins=bins if bins is not None else "auto",
        color="#2f6f73",
        edgecolor="white",
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25, axis="y")
    save_figure(fig, path)
