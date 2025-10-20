import json
import pickle
import tempfile
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd

from simble.cell import Cell
from simble.location import as_enum
from simble.parsing import get_parser, validate_and_process_args
from simble.settings import s
from simble.simble import logger, set_logger
from simble.simulation import simulate
from simble.target import TargetAminoPair
from simble.tree import Node


def do_simulation(i, seed, filename, naive_pickle):

    naive = pickle.loads(naive_pickle)
    with open(filename, "r", encoding="utf-8") as f:
        settings = json.load(f, object_hook=as_enum)
    s.update_from_dict(settings)
    s._x_RNG = np.random.default_rng(seed) # pylint: disable=protected-access
    set_logger()

    print(f"Starting simulation {i}")

    clone_id = i+1
    root = Node(naive, clone_id=clone_id)
    TARGET_PAIR = TargetAminoPair(
        naive.heavy_chain.get_gapped_sequence(),
        naive.light_chain.get_gapped_sequence(),
        naive.heavy_chain.cdr3_length,
        naive.light_chain.cdr3_length)
    TARGET_PAIR.mutate(s.TARGET_MUTATIONS_HEAVY, s.TARGET_MUTATIONS_LIGHT)

    sampled, pop_data, df = simulate(clone_id, TARGET_PAIR, [root], root)

    sampled_cells = [x.cell for x in sampled]

    sampled_pickled = pickle.dumps(sampled_cells)
    target_pickled = pickle.dumps(TARGET_PAIR)

    return {
        "sampled": sampled_pickled,
        "target_pair": target_pickled,
        "pop_data": pop_data,
        "df": df,
        "clone_id": clone_id
        }

def get_affinity_info(cell, clone_id):
    return {"clone_id": clone_id, 
            "cell_id": f"{clone_id}_{id(cell)}",
            "celltype": cell.cell_type.value,
            "affinity": cell.affinity,
            "cross_reactivity": cell.cross_reactivity}


def process_results(results):
    all_results = {}
    for key in results[0].keys():
        all_results[key] = [x[key] for x in results]

    all_results["sampled"] = [pickle.loads(x) for x in all_results["sampled"]]
    all_results["target_pair"] = [pickle.loads(x) for x in all_results["target_pair"]]

    affinity_data = []

    for i in range(len(all_results["clone_id"])):
        print(f"Clone ID: {all_results['clone_id'][i]}")
        sampled = all_results["sampled"][i]
        all_other_targets = [
            x
            for x in all_results["target_pair"]
            if x != all_results["target_pair"][i]
            ]
        if len(all_other_targets) == 0:
            all_other_targets = [all_results["target_pair"][i]]
            print("Only one target, using self for cross-reactivity")
        for cell in sampled:
            real_affinity = cell.affinity
            cross_reactivity = []
            for x in all_other_targets:
                cell.calculate_affinity(x)
                cross_reactivity.append(cell.affinity)
            cell.cross_reactivity = sum(cross_reactivity)/len(cross_reactivity)
            cell.affinity = real_affinity
        affinity_info = [get_affinity_info(x, all_results["clone_id"][i]) for x in sampled]
        affinity_data.extend(affinity_info)

    affinity_data = pd.DataFrame(affinity_data)
    affinity_data.to_csv(s.RESULTS_DIR + "/affinity_data.csv", index=False)


def main():
    parser = get_parser()

    args = parser.parse_args()
    warnings = validate_and_process_args(args)

    set_logger()

    for warning in warnings:
        logger.warning(warning)

    if args.seed is not None:
        seed = args.seed
        ss = np.random.SeedSequence(seed)
    else:
        ss = np.random.SeedSequence()
    seeds = ss.spawn(args.n)
    print(f"Seed: {ss.entropy}")

    naive = Cell(None, None, created_at=0)
    naive_pickle = pickle.dumps(naive)

    with tempfile.NamedTemporaryFile(mode="w") as tmpf:
        json.dump(s, tmpf, default=lambda o: o.encode(), indent=4)
        tmpf.flush()
        logger.info("Starting simulation")
        if args.processes > 1:
            with Pool(processes=args.processes) as pool:
                result = pool.starmap(partial(do_simulation, filename=tmpf.name,
                                              naive_pickle=naive_pickle), zip(range(args.n), seeds))

        else:
            result = []
            for i in range(args.n):
                result.append(do_simulation(i, seeds[i], tmpf.name, naive_pickle))

    process_results(result)

if __name__ == "__main__":
    main()
