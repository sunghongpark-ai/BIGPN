import csv
from pathlib import Path

import numpy as np

from .dataset import PATHWAY_FILE, TARGETS
from .model import integer_option
from .stream import matlab_rand, matlab_randperm, matlab_stream


def generate_sample(folder=None, seed=20261006, num_subject=906):
    seed = integer_option(seed, "seed", 0)
    num_subject = integer_option(num_subject, "num_subject", 1)
    if seed >= 2**32:
        raise ValueError("seed must be smaller than 2**32")
    folder = Path(folder) if folder is not None else PATHWAY_FILE.parent
    counts = [113, 98, 78, 57, 18]
    protein = [f"SYN_P{index + 1:03d}" for index in range(counts[0])]
    levels = [protein] + [[f"SYN_L{level}_{index + 1:03d}" for index in range(count)]
                          for level, count in enumerate(counts[1:], 1)]
    stream = matlab_stream(seed)
    labels = (matlab_rand(stream, len(TARGETS), num_subject) < 0.5).astype(int)
    features = 0.05 + 0.9 * matlab_rand(stream, len(protein), num_subject)
    pairs = [(first, second) for first in range(len(protein)) for second in range(first + 1, len(protein))]
    selected = sorted((pairs[index - 1] for index in matlab_randperm(stream, len(pairs))[:1255]))
    scores = 0.15 + 0.849 * matlab_rand(stream, len(selected), 1)[:, 0]
    links = []
    parent_index = 0
    for level in range(1, len(levels)):
        candidates = [(node, 0) for node in protein]
        if level > 1:
            candidates += [(node, level - 1) for node in levels[level - 1]]
        for parent in levels[level]:
            count = 3 + int(parent_index < 121)
            chosen = matlab_randperm(stream, len(candidates))[:count] - 1
            links.extend((candidates[index][0], candidates[index][1], parent) for index in chosen)
            parent_index += 1
    present = {node for node, _, _ in links}
    links.extend((node, level, "") for level, nodes in enumerate(levels) for node in nodes if node not in present)
    links.sort(key=lambda item: (item[1], item[0], item[2]))
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "sample.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["ID", *["Y" + name for name in TARGETS], *protein])
        for index in range(num_subject):
            writer.writerow([f"SYN_BIGPN_{index + 1:04d}", *labels[:, index], *[f"{value:.6f}" for value in features[:, index]]])
    with (folder / "network.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Protein1", "Protein2", "Score"])
        writer.writerows([protein[first], protein[second], f"{score:.6f}"] for (first, second), score in zip(selected, scores))
    with (folder / "pathway.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Node", "Level", "Parent"])
        writer.writerows(links)
    return folder / "sample.csv"
