from __future__ import annotations

import numpy as np

DEFAULT_SEED = 5489


def matlab_stream(seed: int) -> np.random.RandomState:
    return np.random.RandomState(DEFAULT_SEED if int(seed) == 0 else int(seed))


def matlab_rand(stream: np.random.RandomState, rows: int, cols: int) -> np.ndarray:
    return stream.random_sample(rows * cols).reshape((rows, cols), order="F")


def matlab_randperm(stream: np.random.RandomState, count: int) -> np.ndarray:
    return np.argsort(stream.random_sample(count), kind="stable") + 1
