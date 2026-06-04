"""Plotting helpers for Homework 2."""

from pathlib import Path

import matplotlib.pyplot as plt


def save_figure(fig, path):
    """Save a Matplotlib figure with consistent layout handling."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
