from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import lapack
from scipy.sparse import csr_matrix

from .errors import SingularPropagationError
from .model import Model

EPSILON = np.finfo(float).eps


@dataclass(frozen=True)
class Param:
    u: np.ndarray
    b: np.ndarray
    w: tuple[np.ndarray, ...]
    s: tuple[np.ndarray, ...]
    operator: tuple[csr_matrix, ...]
    inverse: np.ndarray


@dataclass(frozen=True)
class Cache:
    input: np.ndarray
    gene: np.ndarray
    source: tuple[np.ndarray, ...]
    path: tuple[np.ndarray, ...]
    logit: np.ndarray
    probability: np.ndarray


def reshape_weight(model: Model, weight: np.ndarray) -> Param:
    u = weight[model.index_u]
    b = weight[model.index_b].reshape((model.num_path, model.num_target), order="F")
    w = tuple(weight[block] for block in model.index_w)
    s = tuple(1.0 / (1.0 + np.exp(-value)) for value in w)
    operator = tuple(
        csr_matrix((value[level.order], level.indices, level.indptr), shape=(level.num_row, level.num_col))
        for level, value in zip(model.levels, s)
    )
    propagation = model.laplacian + np.diag(u)
    lu, pivot, info = lapack.dgetrf(propagation)
    condition = 0.0
    if info == 0:
        condition = lapack.dgecon(lu, float(np.abs(propagation).sum(axis=0).max()), norm="1")[0]
    if not condition >= EPSILON:
        raise SingularPropagationError(
            "The propagation matrix diag(u) + L is singular to working precision; lower learn_rate."
        )
    return Param(u, b, w, s, operator, lapack.dgetri(lu, pivot)[0])


def solve(param: Param, right: np.ndarray, transpose: bool = False) -> np.ndarray:
    return (param.inverse.T if transpose else param.inverse) @ right


def forward(model: Model, param: Param, index: np.ndarray) -> Cache:
    data = model.x[:, index]
    gene = solve(param, param.u[:, None] * data)
    logit = param.b[model.gene_index].T @ gene
    previous = np.zeros((0, gene.shape[1]))
    source = []
    path = []
    for level, operator in zip(model.levels, param.operator):
        stacked = np.vstack((gene, previous))
        current = operator @ stacked
        logit = logit + param.b[level.index].T @ current
        source.append(stacked)
        path.append(current)
        previous = current
    return Cache(data, gene, tuple(source), tuple(path), logit, 1.0 / (1.0 + np.exp(-logit)))


def loss(logit: np.ndarray, label: np.ndarray) -> float:
    return float(np.sum(np.mean(np.maximum(logit, 0.0) - label * logit + np.log1p(np.exp(-np.abs(logit))), axis=1)))


def backward(model: Model, param: Param, cache: Cache, label: np.ndarray, weight: np.ndarray) -> np.ndarray:
    delta = (cache.probability - label) / label.shape[1]
    gradient = 2.0 * model.reg_gamma * weight
    legacy = model.gradient == "legacy"
    grad_b = np.zeros((model.num_path, model.num_target))
    grad_b[model.gene_index] = cache.gene @ delta.T
    error_gene = param.b[model.gene_index] @ delta
    error_path = []
    for level, current in zip(model.levels, cache.path):
        grad_b[level.index] = current @ delta.T
        error_path.append(param.b[level.index] @ delta)
    for depth in range(model.num_level - 1, -1, -1):
        level = model.levels[depth]
        outer = error_path[depth] @ cache.source[depth].T
        s = param.s[depth]
        gradient[model.index_w[depth]] += outer[level.row, level.col] * (s * (1.0 - s))
        back = param.operator[depth].T @ error_path[depth]
        if depth == 0 or not legacy:
            error_gene = error_gene + back[: model.num_gene]
        if depth > 0:
            error_path[depth - 1] = error_path[depth - 1] + back[model.num_gene :]
    if legacy:
        adjoint = solve(param, cache.input @ error_gene.T, transpose=True)
        grad_u = np.diag(adjoint) - np.diag(solve(param, adjoint, transpose=True)) * param.u
    else:
        grad_u = np.sum(solve(param, error_gene, transpose=True) * (cache.input - cache.gene), axis=1)
    gradient[model.index_u] += grad_u
    gradient[model.index_b] += grad_b.ravel(order="F")
    return gradient
