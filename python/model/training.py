from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .errors import TrainingFailedError
from .model import Model
from .optimizer import Adam
from .propagation import backward, forward, loss, reshape_weight
from .stream import matlab_rand, matlab_stream


@dataclass(frozen=True)
class Split:
    model_index: int
    iteration: int
    test_fold: int
    valid_fold: int
    train_folds: np.ndarray
    train: np.ndarray
    valid: np.ndarray
    test: np.ndarray


@dataclass(frozen=True)
class Fit:
    loss_train: np.ndarray
    loss_valid: np.ndarray
    best_epoch: int
    weight: np.ndarray
    probability: np.ndarray


def split_indices(model: Model, model_index: int) -> Split:
    iteration, test_fold, valid_fold = (int(value) for value in model.cv_list[model_index])
    train_folds = np.setdiff1d(np.arange(1, model.num_fold + 1), [test_fold, valid_fold])
    fold = model.cv_data[iteration - 1]
    return Split(
        model_index=model_index,
        iteration=iteration,
        test_fold=test_fold,
        valid_fold=valid_fold,
        train_folds=train_folds,
        train=np.flatnonzero(np.isin(fold, train_folds)),
        valid=np.flatnonzero(fold == valid_fold),
        test=np.flatnonzero(fold == test_fold),
    )


def initialize_weight(model: Model, iteration: int) -> np.ndarray:
    stream = matlab_stream(int(model.seed[iteration - 1]))
    weight = np.zeros(model.num_param)
    weight[model.index_u] = 1.0
    for level, block in zip(model.levels, model.index_w):
        draw = matlab_rand(stream, level.num_row, model.num_path)
        weight[block] = (2.0 * draw[level.row, level.child] - 1.0) * math.sqrt(6.0 / (level.num_row + model.num_path))
    draw = matlab_rand(stream, model.num_path, model.num_target)
    weight[model.index_b] = (2.0 * draw.ravel(order="F") - 1.0) * math.sqrt(6.0 / (model.num_path + 1))
    return weight


def train(model: Model, split: Split, weight: np.ndarray) -> Fit:
    adam = Adam(model.num_param, model.learn_rate)
    label_train = model.y[:, split.train]
    label_valid = model.y[:, split.valid]
    loss_train = np.full(model.max_epoch, np.nan)
    loss_valid = np.full(model.max_epoch, np.nan)
    best_epoch = 0
    best_weight = weight
    best_loss = math.inf
    for epoch in range(1, model.max_epoch + 1):
        param = reshape_weight(model, weight)
        cache_train = forward(model, param, split.train)
        cache_valid = forward(model, param, split.valid)
        loss_train[epoch - 1] = loss(cache_train.logit, label_train)
        loss_valid[epoch - 1] = loss(cache_valid.logit, label_valid)
        if loss_valid[epoch - 1] < best_loss:
            best_loss = loss_valid[epoch - 1]
            best_epoch = epoch
            best_weight = weight
        if epoch == model.max_epoch:
            break
        weight = adam.update(weight, backward(model, param, cache_train, label_train, weight))
        if not np.all(np.isfinite(weight)):
            break
    if best_epoch == 0:
        raise TrainingFailedError(f"Validation loss was not finite in any epoch of model {split.model_index + 1}.")
    cache = forward(model, reshape_weight(model, best_weight), np.arange(model.num_subject))
    return Fit(loss_train, loss_valid, best_epoch, best_weight, cache.probability)
