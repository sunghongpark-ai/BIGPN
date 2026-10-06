from __future__ import annotations

import csv
import os
import time
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

try:
    from threadpoolctl import threadpool_limits
except ImportError:
    threadpool_limits = None

from .dataset import TARGETS, read_dataset
from .errors import InvalidOptionError
from .model import Model, Options, build_model
from .training import Fit, initialize_weight, split_indices, train

_WORKER: dict[str, Any] = {}


@dataclass(frozen=True)
class Result:
    model: Model
    risk: dict[str, np.ndarray]
    test_risk: dict[str, np.ndarray]
    weight: np.ndarray = field(repr=False)
    best_epoch: np.ndarray = field(repr=False)
    loss_train: np.ndarray = field(repr=False)
    loss_valid: np.ndarray = field(repr=False)
    elapsed_time: float
    options: Options

    def save(self, file: str | os.PathLike[str]) -> None:
        array = {
            "subject": np.array(self.model.subject),
            "target": np.array(self.model.target),
            "cv_data": self.model.cv_data,
            "cv_list": self.model.cv_list,
            "weight": self.weight,
            "best_epoch": self.best_epoch,
            "loss_train": self.loss_train,
            "loss_valid": self.loss_valid,
        }
        for name in self.model.target:
            array[f"risk_{name}"] = self.risk[name]
            array[f"test_risk_{name}"] = self.test_risk[name]
        np.savez_compressed(file, **array)


def run(
    cohort: str | os.PathLike[str],
    *,
    max_epoch: int,
    learn_rate: float,
    reg_gamma: float,
    network: str | os.PathLike[str] | None = None,
    pathway: str | os.PathLike[str] | None = None,
    num_iter: int = 100,
    num_fold: int = 5,
    seed: int = 1,
    targets: Sequence[str] | str = TARGETS,
    edge_threshold: float = 0.15,
    gradient: str = "exact",
    n_jobs: int | None = None,
    output: str | os.PathLike[str] | None = None,
    verbose: bool = False,
) -> Result:
    options = Options(max_epoch, learn_rate, reg_gamma, num_iter, num_fold, seed, edge_threshold, gradient)
    workers = worker_count(n_jobs)
    model = build_model(read_dataset(cohort, network, pathway, targets), options)
    clock = time.perf_counter()
    fits = fit_all(model, workers, report if verbose else None)
    elapsed = time.perf_counter() - clock
    probability = np.stack([fit.probability for fit in fits], axis=2)
    risk = {name: probability[index].T.copy() for index, name in enumerate(model.target)}
    result = Result(
        model=model,
        risk=risk,
        test_risk=summarize_test_risk(model, probability),
        weight=np.stack([fit.weight for fit in fits]),
        best_epoch=np.array([fit.best_epoch for fit in fits], dtype=np.int64),
        loss_train=np.stack([fit.loss_train for fit in fits]),
        loss_valid=np.stack([fit.loss_valid for fit in fits]),
        elapsed_time=elapsed,
        options=options,
    )
    if output is not None:
        write_risk(output, model, result.test_risk)
    return result


def worker_count(n_jobs: int | None) -> int:
    if n_jobs is None:
        return 1
    if isinstance(n_jobs, bool) or not isinstance(n_jobs, int) or n_jobs == 0 or n_jobs < -1:
        raise InvalidOptionError("n_jobs must be None, -1 or a positive integer.")
    return (os.cpu_count() or 1) if n_jobs == -1 else n_jobs


def fit_model(model: Model, model_index: int) -> Fit:
    split = split_indices(model, model_index)
    return train(model, split, initialize_weight(model, split.iteration))


def fit_all(model: Model, workers: int, callback: Callable[[Model, int, Fit], None] | None = None) -> list[Fit]:
    fits: dict[int, Fit] = {}
    if workers == 1:
        for model_index in range(model.num_model):
            fits[model_index] = fit = fit_model(model, model_index)
            if callback is not None:
                callback(model, model_index, fit)
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=worker_setup, initargs=(model,)) as pool:
            pending = {pool.submit(worker_fit, model_index): model_index for model_index in range(model.num_model)}
            try:
                for future in as_completed(pending):
                    model_index = pending[future]
                    fits[model_index] = fit = future.result()
                    if callback is not None:
                        callback(model, model_index, fit)
            except BaseException:
                pool.shutdown(wait=False, cancel_futures=True)
                raise
    return [fits[model_index] for model_index in range(model.num_model)]


def worker_setup(model: Model) -> None:
    _WORKER["model"] = model
    if threadpool_limits is not None:
        _WORKER["limit"] = threadpool_limits(limits=1)


def worker_fit(model_index: int) -> Fit:
    return fit_model(_WORKER["model"], model_index)


def report(model: Model, model_index: int, fit: Fit) -> None:
    iteration, test_fold, valid_fold = (int(value) for value in model.cv_list[model_index])
    print(
        f"Model {model_index + 1}/{model.num_model} (iteration {iteration}, test fold {test_fold}, "
        f"validation fold {valid_fold}): best epoch {fit.best_epoch}, "
        f"validation loss {fit.loss_valid[fit.best_epoch - 1]:.6f}",
        flush=True,
    )


def summarize_test_risk(model: Model, probability: np.ndarray) -> dict[str, np.ndarray]:
    total = np.zeros((model.num_target, model.num_subject, model.num_iter))
    count = np.zeros((1, model.num_subject, model.num_iter))
    for model_index, (iteration, test_fold, _) in enumerate(model.cv_list):
        test = model.cv_data[iteration - 1] == test_fold
        total[:, test, iteration - 1] = total[:, test, iteration - 1] + probability[:, test, model_index]
        count[0, test, iteration - 1] = count[0, test, iteration - 1] + 1
    average = total / count
    return {name: average[index].T.copy() for index, name in enumerate(model.target)}


def write_risk(file: str | os.PathLike[str], model: Model, test_risk: dict[str, np.ndarray]) -> None:
    value = np.column_stack([test_risk[name].mean(axis=0) for name in model.target])
    with open(Path(file), "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["ID", *(f"P{name}" for name in model.target)])
        for subject, row in zip(model.subject, value):
            writer.writerow([subject, *(format(number, ".15g") for number in row)])
