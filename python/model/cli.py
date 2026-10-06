from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from .dataset import TARGETS
from .errors import BIGPNError
from .runner import run


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        prog="bigpn",
        description="Train BIGPN with repeated stratified cross-validation on a cohort CSV file.",
    )
    command.add_argument(
        "cohort", nargs="?", default=Path(__file__).resolve().parents[2] / "dataset" / "sample.csv", type=Path, help="cohort CSV with ID, Y<target> label columns and one column per protein"
    )
    command.add_argument("--max-epoch", type=int, default=20, help="maximum number of training epochs")
    command.add_argument("--learn-rate", type=float, default=0.001, help="Adam learning rate")
    command.add_argument("--reg-gamma", type=float, default=0.0001, help="L2 regularization coefficient")
    command.add_argument("--network", type=Path, help="PPI network CSV (default: synthetic network scores)")
    command.add_argument("--pathway", type=Path, help="pathway hierarchy CSV (default: synthetic pathway topology)")
    command.add_argument(
        "--num-iter", type=int, default=1, help="number of cross-validation repetitions (default: 1)"
    )
    command.add_argument("--num-fold", type=int, default=5, help="number of folds (default: 5)")
    command.add_argument("--seed", type=int, default=1, help="seed of the first repetition (default: 1)")
    command.add_argument("--targets", nargs="+", default=list(TARGETS), help="target names (default: abt gfa nfl tau)")
    command.add_argument("--edge-threshold", type=float, default=0.15, help="lowest PPI score kept (default: 0.15)")
    command.add_argument(
        "--gradient", choices=("exact", "legacy"), default="exact", help="gradient of the propagation parameters"
    )
    command.add_argument("--jobs", type=int, default=None, help="worker processes; -1 uses every CPU (default: 1)")
    command.add_argument("--output", type=Path, help="write ID and mean test-fold risk per target to this CSV")
    command.add_argument("--save", type=Path, help="write all results to this .npz file")
    command.add_argument("--verbose", action="store_true", help="print one line per trained model")
    return command


def main(argv: Sequence[str] | None = None) -> int:
    command = parser()
    argument = command.parse_args(argv)
    try:
        result = run(
            argument.cohort,
            max_epoch=argument.max_epoch,
            learn_rate=argument.learn_rate,
            reg_gamma=argument.reg_gamma,
            network=argument.network,
            pathway=argument.pathway,
            num_iter=argument.num_iter,
            num_fold=argument.num_fold,
            seed=argument.seed,
            targets=argument.targets,
            edge_threshold=argument.edge_threshold,
            gradient=argument.gradient,
            n_jobs=argument.jobs,
            output=argument.output,
            verbose=argument.verbose,
        )
    except BIGPNError as error:
        command.exit(1, f"bigpn: error: {error}\n")
    if argument.save is not None:
        result.save(argument.save)
    model = result.model
    print(f"Trained {model.num_model} models on {model.num_subject} subjects in {result.elapsed_time:.1f} s.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
