from __future__ import annotations

import csv
import math
import os
import re
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .errors import InvalidDatasetError, InvalidOptionError

TARGETS = ("abt", "gfa", "nfl", "tau")
RESOURCE = Path(__file__).resolve().parents[2] / "dataset"
NETWORK_FILE = RESOURCE / "network.csv"
PATHWAY_FILE = RESOURCE / "pathway.csv"
NUMBER = re.compile(r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
IDENTIFIER = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,62}")
MATLAB_KEYWORDS = frozenset(
    {
        "break",
        "case",
        "catch",
        "classdef",
        "continue",
        "else",
        "elseif",
        "end",
        "for",
        "function",
        "global",
        "if",
        "otherwise",
        "parfor",
        "persistent",
        "return",
        "spmd",
        "switch",
        "try",
        "while",
    }
)


@dataclass(frozen=True)
class Dataset:
    target: tuple[str, ...]
    node: tuple[str, ...]
    depth: np.ndarray
    link: np.ndarray
    protein: tuple[str, ...]
    network: np.ndarray
    subject: tuple[str, ...]
    x: np.ndarray
    y: np.ndarray


def read_dataset(
    cohort: str | os.PathLike[str],
    network: str | os.PathLike[str] | None = None,
    pathway: str | os.PathLike[str] | None = None,
    targets: Sequence[str] | str = TARGETS,
) -> Dataset:
    target = check_targets(targets)
    node, depth, link = read_pathway(PATHWAY_FILE if pathway is None else pathway)
    protein = tuple(name for name, level in zip(node, depth) if level == 0)
    weight = read_network(NETWORK_FILE if network is None else network, protein)
    subject, x, y = read_cohort(cohort, protein, target)
    return Dataset(target, node, depth, link, protein, weight, subject, x, y)


def check_targets(targets: Sequence[str] | str) -> tuple[str, ...]:
    if isinstance(targets, str):
        targets = (targets,)
    try:
        target = tuple(targets)
    except TypeError:
        target = ()
    if not target or not all(isinstance(name, str) for name in target):
        raise InvalidOptionError("targets must be nonempty text.")
    valid = all(IDENTIFIER.fullmatch(name) and name not in MATLAB_KEYWORDS for name in target)
    if not valid or len(set(target)) != len(target):
        raise InvalidOptionError(
            "targets must be unique identifiers that start with a letter, contain only letters, digits and "
            "underscores, have at most 63 characters and are not MATLAB keywords."
        )
    return target


def read_pathway(file: str | os.PathLike[str]) -> tuple[tuple[str, ...], np.ndarray, np.ndarray]:
    name = os.fspath(file)
    header, body = read_csv(name)
    column = column_index(name, header, ("Node", "Level", "Parent"))
    node = [row[column[0]] for row in body]
    parent = [row[column[2]] for row in body]
    level = parse_numbers(name, [[row[column[1]]] for row in body], "Level", 1)[:, 0]
    if not node:
        raise InvalidDatasetError(f"{name} has no rows.")
    if "" in node:
        raise InvalidDatasetError(f"{name}: Node must not be empty (line {node.index('') + 2}).")
    invalid = (level != np.round(level)) | (level < 0)
    if invalid.any():
        raise InvalidDatasetError(f"{name}: Level must be a nonnegative integer (line {int(np.argmax(invalid)) + 2}).")
    declared: dict[str, set[int]] = {}
    for item, value in zip(node, level):
        declared.setdefault(item, set()).add(int(value))
    conflict = sorted(item for item, values in declared.items() if len(values) > 1)
    if conflict:
        raise InvalidDatasetError(f"{name}: nodes declared at more than one level: {name_list(conflict)}.")
    node_level = {item: values.pop() for item, values in declared.items()}
    top = max(node_level.values())
    if top < 1 or set(range(top + 1)) - set(node_level.values()):
        raise InvalidDatasetError(f"{name}: levels must run from 0 (proteins) to the top pathway level without gaps.")
    order = sorted(node_level, key=lambda item: (node_level[item], item))
    position = {item: index for index, item in enumerate(order)}
    depth = np.array([node_level[item] for item in order], dtype=np.int64)
    if len(set(zip(node, parent))) < len(node):
        raise InvalidDatasetError(f"{name} has duplicate Node-Parent rows.")
    unknown = sorted({item for item in parent if item and item not in position})
    if unknown:
        raise InvalidDatasetError(f"{name}: parents not declared as nodes: {name_list(unknown)}.")
    link = np.array(
        [(position[child], position[upper]) for child, upper in zip(node, parent) if upper], dtype=np.int64
    ).reshape(-1, 2)
    parent_depth = depth[link[:, 1]]
    child_depth = depth[link[:, 0]]
    if (parent_depth < 1).any():
        protein_parent = sorted({order[index] for index in link[parent_depth < 1, 1]})
        raise InvalidDatasetError(f"{name}: proteins (level 0) cannot be parents: {name_list(protein_parent)}.")
    skip = (child_depth != 0) & (child_depth != parent_depth - 1)
    if skip.any():
        raise InvalidDatasetError(
            f"{name} has {int(skip.sum())} links whose child is neither a protein nor a pathway one level "
            "below its parent."
        )
    childless = sorted(set(np.flatnonzero(depth >= 1).tolist()) - set(link[:, 1].tolist()))
    if childless:
        raise InvalidDatasetError(
            f"{name}: pathways without children: {name_list([order[index] for index in childless])}."
        )
    link = link[np.lexsort((link[:, 0], link[:, 1]))]
    return tuple(order), depth, link


def read_network(file: str | os.PathLike[str], protein: Sequence[str]) -> np.ndarray:
    name = os.fspath(file)
    header, body = read_csv(name)
    column = column_index(name, header, ("Protein1", "Protein2", "Score"))
    position = {item: index for index, item in enumerate(protein)}
    unknown = sorted({row[column[side]] for row in body for side in (0, 1) if row[column[side]] not in position})
    if unknown:
        raise InvalidDatasetError(f"{name}: proteins absent from the pathway file: {name_list(unknown)}.")
    score = parse_numbers(name, [[row[column[2]]] for row in body], "Score", 1)[:, 0]
    if (score < 0).any():
        raise InvalidDatasetError(f"{name}: Score must be nonnegative (line {int(np.argmax(score < 0)) + 2}).")
    first = np.array([position[row[column[0]]] for row in body], dtype=np.int64)
    second = np.array([position[row[column[1]]] for row in body], dtype=np.int64)
    if (first == second).any():
        raise InvalidDatasetError(
            f"{name}: self-interactions are not allowed (line {int(np.argmax(first == second)) + 2})."
        )
    pair = np.sort(np.column_stack((first, second)), axis=1)
    if len({tuple(item) for item in pair.tolist()}) < len(pair):
        raise InvalidDatasetError(f"{name} lists a protein pair more than once.")
    weight = np.zeros((len(protein), len(protein)))
    weight[pair[:, 0], pair[:, 1]] = score
    return weight + weight.T


def read_cohort(
    file: str | os.PathLike[str], protein: Sequence[str], target: Sequence[str]
) -> tuple[tuple[str, ...], np.ndarray, np.ndarray]:
    name = os.fspath(file)
    header, body = read_csv(name)
    label = tuple(f"Y{item}" for item in target)
    expected = ("ID", *label, *protein)
    if len(set(expected)) < len(expected):
        raise InvalidOptionError(f"The column names ID, {', '.join(label)} and the protein names must be distinct.")
    column = column_index(name, header, expected)
    subject = tuple(row[column[0]] for row in body)
    if not subject:
        raise InvalidDatasetError(f"{name} has no subjects.")
    if "" in subject:
        raise InvalidDatasetError(f"{name}: ID must not be empty (line {subject.index('') + 2}).")
    if len(set(subject)) < len(subject):
        repeat = sorted(item for item, count in Counter(subject).items() if count > 1)
        raise InvalidDatasetError(f"{name}: duplicate IDs: {name_list(repeat)}.")
    count = len(target)
    y = parse_numbers(name, [[row[index] for index in column[1 : count + 1]] for row in body], "label", count).T
    invalid = (y != 0) & (y != 1)
    if invalid.any():
        which, line = np.argwhere(invalid.T)[0][::-1]
        raise InvalidDatasetError(f"{name}: {label[which]} must be 0 or 1 (line {line + 2}).")
    x = parse_numbers(name, [[row[index] for index in column[count + 1 :]] for row in body], "protein", len(protein)).T
    return subject, np.ascontiguousarray(x), np.ascontiguousarray(y)


def column_index(name: str, header: Sequence[str], expected: Sequence[str]) -> list[int]:
    if len(set(header)) < len(header):
        repeat = sorted(item for item, count in Counter(header).items() if count > 1)
        raise InvalidDatasetError(f"{name}: duplicate column names: {name_list(repeat)}.")
    missing = [item for item in expected if item not in header]
    if missing:
        raise InvalidDatasetError(f"{name}: missing columns: {name_list(missing)}.")
    extra = sorted(set(header) - set(expected))
    if extra:
        raise InvalidDatasetError(f"{name}: unexpected columns: {name_list(extra)}.")
    return [header.index(item) for item in expected]


def read_csv(name: str) -> tuple[list[str], list[list[str]]]:
    try:
        with open(name, newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle, strict=True)
            try:
                rows = list(reader)
            except csv.Error as error:
                raise InvalidDatasetError(f"{name}: line {reader.line_num} is not valid CSV ({error}).") from None
    except OSError as error:
        raise InvalidDatasetError(f"Cannot open {name}.") from error
    except UnicodeDecodeError as error:
        raise InvalidDatasetError(f"{name} is not UTF-8 text.") from error
    while rows and not rows[-1]:
        rows.pop()
    if not rows:
        raise InvalidDatasetError(f"{name} is empty.")
    width = len(rows[0])
    for number, fields in enumerate(rows, 1):
        if not fields:
            raise InvalidDatasetError(f"{name}: line {number} is empty.")
        if any("\n" in value or "\r" in value for value in fields):
            raise InvalidDatasetError(f"{name}: line {number} has a line break inside a quoted field.")
        if len(fields) != width:
            raise InvalidDatasetError(f"{name}: line {number} has {len(fields)} fields but the header has {width}.")
    cell = [[value.strip(" \t") for value in fields] for fields in rows]
    if "" in cell[0]:
        raise InvalidDatasetError(f"{name} has an empty column name.")
    return cell[0], cell[1:]


def parse_numbers(name: str, cell: Sequence[Sequence[str]], what: str, width: int) -> np.ndarray:
    value = []
    for row, fields in enumerate(cell):
        number = [float(text) if NUMBER.fullmatch(text) else math.nan for text in fields]
        for text, item in zip(fields, number):
            if not math.isfinite(item):
                raise InvalidDatasetError(f"{name}: {what} value '{text}' on line {row + 2} is not a finite number.")
        value.append(number)
    return np.array(value, dtype=float).reshape(len(cell), width)


def name_list(name: Iterable[str]) -> str:
    item = list(name)
    text = ", ".join(item[:5])
    if len(item) > 5:
        text = f"{text} and {len(item) - 5} more"
    return text
