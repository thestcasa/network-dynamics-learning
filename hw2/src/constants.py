"""Reusable constants for Homework 2.

Matrix convention: Lambda[i, j] is the rate/weight from source node i to destination node j.
"""

from pathlib import Path

import numpy as np

NODES = ["o", "a", "b", "c", "d"]
NODE_TO_INDEX = {node: index for index, node in enumerate(NODES)}

LAMBDA = np.array(
    [
        [0.0, 2.0 / 5.0, 1.0 / 5.0, 0.0, 0.0],
        [0.0, 0.0, 3.0 / 4.0, 1.0 / 4.0, 0.0],
        [1.0 / 2.0, 0.0, 0.0, 1.0 / 3.0, 0.0],
        [0.0, 0.0, 1.0 / 3.0, 0.0, 2.0 / 3.0],
        [0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0],
    ],
    dtype=float,
)

LAMBDA_OPEN = np.array(
    [
        [0.0, 1.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0 / 4.0, 1.0 / 4.0, 2.0 / 4.0],
        [0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
    ],
    dtype=float,
)

DEFAULT_SEED = 20260604

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"

