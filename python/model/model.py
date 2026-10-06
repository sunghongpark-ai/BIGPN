from __future__ import annotations

import itertools
import math
import numbers
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .dataset import Dataset
from .errors import InvalidDatasetError, InvalidOptionError
from .stream import matlab_randperm, matlab_stream

GRADIENTS = ("exact", "legacy")


@dataclass(frozen=True)
class Options:
    max_epoch: int
    learn_rate: float
    reg_gamma: float
    num_iter: int = 100
    num_fold: int = 5
    seed: int = 1
    edge_threshold: float = 0.15
    gradient: str = "exact"

    def __post_init__(self) -> None:
        object.__setattr__(self, "num_iter", integer_option(self.num_iter, "num_iter", 1))
        object.__setattr__(self, "num_fold", integer_option(self.num_fold, "num_fold", 3))
        object.__setattr__(self, "max_epoch", integer_option(self.max_epoch, "max_epoch", 1))
        object.__setattr__(self, "seed", integer_option(self.seed, "seed", 0))
        if self.seed + self.num_iter - 1 > 2**32 - 1:
            raise InvalidOptionError("seed + num_iter - 1 must not exceed 2**32 - 1.")
        object.__setattr__(self, "learn_rate", real_option(self.learn_rate, "learn_rate", False))
        object.__setattr__(self, "reg_gamma", real_option(self.reg_gamma, "reg_gamma", True))
        object.__setattr__(self, "edge_threshold", real_option(self.edge_threshold, "edge_threshold", True))
        if not isinstance(self.gradient, str) or self.gradient not in GRADIENTS:
            raise InvalidOptionError("gradient must be 'exact' or 'legacy'.")


@dataclass(frozen=True)
class Level:
    index: np.ndarray
    num_row: int
    num_col: int
    parent: np.ndarray
    child: np.ndarray
    row: np.ndarray
    col: np.ndarray
    indptr: np.ndarray
    indices: np.ndarray
    order: np.ndarray


@dataclass(frozen=True)
class Model:
    target: tuple[str, ...]
    num_iter: int
    num_fold: int
    max_epoch: int
    learn_rate: float
    reg_gamma: float
    edge_threshold: float
    gradient: str
    seed: np.ndarray
    subject: tuple[str, ...]
    protein: tuple[str, ...]
    x: np.ndarray
    y: np.ndarray
    laplacian: np.ndarray
    node: tuple[str, ...]
    gene_index: np.ndarray
    levels: tuple[Level, ...]
    index_u: slice
    index_w: tuple[slice, ...]
    index_b: slice
    num_param: int
    cv_data: np.ndarray = field(repr=False)
    cv_list: np.ndarray = field(repr=False)

    @property
    def num_target(self) -> int:
        return len(self.target)

    @property
    def num_gene(self) -> int:
        return self.x.shape[0]

    @property
    def num_subject(self) -> int:
        return self.x.shape[1]

    @property
    def num_path(self) -> int:
        return len(self.node)

    @property
    def num_level(self) -> int:
        return len(self.levels)

    @property
    def num_model(self) -> int:
        return len(self.cv_list)


def build_model(dataset: Dataset, options: Options) -> Model:
    laplacian = network_laplacian(dataset.network, options.edge_threshold)
    gene_index, levels = pathway_structure(dataset.node, dataset.depth, dataset.link, len(dataset.protein))
    count = [len(dataset.protein), *(len(level.row) for level in levels), len(dataset.node) * len(dataset.target)]
    edge = np.cumsum([0, *count])
    block = [slice(int(edge[index]), int(edge[index + 1])) for index in range(len(count))]
    seed = options.seed + np.arange(options.num_iter, dtype=np.int64)
    cv_data, cv_list = cross_validation(dataset.y, options.num_iter, options.num_fold, seed)
    return Model(
        target=dataset.target,
        num_iter=options.num_iter,
        num_fold=options.num_fold,
        max_epoch=options.max_epoch,
        learn_rate=options.learn_rate,
        reg_gamma=options.reg_gamma,
        edge_threshold=options.edge_threshold,
        gradient=options.gradient,
        seed=seed,
        subject=dataset.subject,
        protein=dataset.protein,
        x=dataset.x,
        y=dataset.y,
        laplacian=laplacian,
        node=dataset.node,
        gene_index=gene_index,
        levels=levels,
        index_u=block[0],
        index_w=tuple(block[1:-1]),
        index_b=block[-1],
        num_param=int(edge[-1]),
        cv_data=cv_data,
        cv_list=cv_list,
    )


def integer_option(value: Any, name: str, lower: int) -> int:
    valid = isinstance(value, numbers.Real) and not isinstance(value, bool) and math.isfinite(value)
    if not valid or value != round(value) or value < lower:
        raise InvalidOptionError(f"{name} must be an integer of at least {lower}.")
    return int(value)


def real_option(value: Any, name: str, allow_zero: bool) -> float:
    valid = isinstance(value, numbers.Real) and not isinstance(value, bool) and math.isfinite(value)
    if not valid or not (value > 0 or (allow_zero and value == 0)):
        bound = "nonnegative" if allow_zero else "positive"
        raise InvalidOptionError(f"{name} must be a finite {bound} scalar.")
    return float(value)


def network_laplacian(weight: np.ndarray, threshold: float) -> np.ndarray:
    count = weight.shape[0]
    flat = np.where(weight < threshold, 0.0, weight).ravel(order="F")
    edge = flat > 0
    if edge.any():
        value = flat[edge]
        score = np.zeros_like(value)
        spread = float(np.std(value, ddof=1)) if value.size > 1 else 0.0
        if spread > 0:
            score = (value - value.mean()) / spread
        flat[edge] = 1.0 / (1.0 + np.exp(-score))
    weight = flat.reshape((count, count), order="F")
    degree = weight.sum(axis=1)
    scale = np.zeros(count)
    scale[degree > 0] = 1.0 / np.sqrt(degree[degree > 0])
    return np.eye(count) - (scale[:, None] * weight) * scale[None, :]


def pathway_structure(
    node: tuple[str, ...], depth: np.ndarray, link: np.ndarray, num_gene: int
) -> tuple[np.ndarray, tuple[Level, ...]]:
    num_level = int(depth.max())
    node_row = np.zeros(len(node), dtype=np.int64)
    member = []
    for value in range(num_level + 1):
        index = np.flatnonzero(depth == value)
        node_row[index] = np.arange(index.size)
        member.append(index)
    pair = np.unique(link[:, ::-1], axis=0)
    parent_depth = depth[pair[:, 0]]
    levels = []
    for value in range(1, num_level + 1):
        num_row = member[value].size
        num_col = num_gene + (member[value - 1].size if value > 1 else 0)
        parent = pair[parent_depth == value, 0]
        child = pair[parent_depth == value, 1]
        order = np.argsort(node_row[parent] + child * num_row, kind="stable")
        parent = parent[order]
        child = child[order]
        row = node_row[parent]
        col = node_row[child] + num_gene * (depth[child] > 0)
        sparse = np.lexsort((col, row))
        indptr = np.concatenate(([0], np.cumsum(np.bincount(row, minlength=num_row))))
        levels.append(Level(member[value], num_row, num_col, parent, child, row, col, indptr, col[sparse], sparse))
    return member[0], tuple(levels)


def cross_validation(y: np.ndarray, num_iter: int, num_fold: int, seed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    pattern = label_pattern(y.shape[0])
    lookup = {item: index for index, item in enumerate(pattern)}
    group = np.array([lookup[tuple(int(value) for value in column)] for column in y.T])
    member = [np.flatnonzero(group == index) for index in range(len(pattern))]
    cv_data = np.zeros((num_iter, y.shape[1]), dtype=np.int64)
    for iteration in range(num_iter):
        stream = matlab_stream(int(seed[iteration]))
        for index in member:
            cv_data[iteration, index] = np.mod(matlab_randperm(stream, index.size), num_fold) + 1
    for fold in range(1, num_fold + 1):
        empty = np.flatnonzero(~(cv_data == fold).any(axis=1))
        if empty.size:
            raise InvalidDatasetError(f"Fold {fold} of iteration {empty[0] + 1} is empty; reduce num_fold.")
    cv_list = np.array(
        [
            (iteration, test, valid)
            for iteration in range(1, num_iter + 1)
            for test in range(1, num_fold + 1)
            for valid in range(1, num_fold + 1)
            if valid != test
        ],
        dtype=np.int64,
    )
    return cv_data, cv_list


def label_pattern(num_target: int) -> list[tuple[int, ...]]:
    pattern = []
    for num_positive in range(num_target + 1):
        positive_minority = 2 * num_positive <= num_target
        num_minority = min(num_positive, num_target - num_positive)
        for minority in itertools.combinations(range(num_target), num_minority):
            row = [int(not positive_minority)] * num_target
            for position in minority:
                row[position] = int(positive_minority)
            pattern.append(tuple(row))
    return pattern
