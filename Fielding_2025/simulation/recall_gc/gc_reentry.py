#############################################
## Hunter J. Melton and Jessie J. Fielding ##
################## 8/28/25 ##################

# Script to simulate a recall GC reaction using simble
# Initial GC reaction runs from generation 1-100, one other B cell is randomly chosen to wait 1000 gens,
# then re-enters a new GC and that reaction goes from generation 1101-1200

import argparse
import json
import logging
import os
from time import time
from collections import namedtuple
from functools import partial
from multiprocessing import Pool
import tempfile

import numpy as np
import pandas as pd
from simble.cell import Cell, CellType
import simble.helper as helper
from simble.location import LocationName, as_enum
from simble.settings import s
from simble.simble import process_results, set_logger
from simble.simulation import simulate
from simble.target import TargetAminoPair
from simble.tree import Node, simplify_tree
from simble.parsing import get_parser, validate_and_process_args

logger = logging.getLogger(__package__)

def do_simulation(i, seed, filename):
    """Runs a single simulation with the given seed and settings."""
    with open(filename, "r", encoding="utf-8") as f:
        settings = json.load(f, object_hook=as_enum)
    s.update_from_dict(settings)
    s._x_RNG = np.random.default_rng(seed) # pylint: disable=protected-access
    set_logger()
    logger.info("Starting simulation %s", i)

    if not os.path.exists(s.RESULTS_DIR):
        os.mkdir(s.RESULTS_DIR)
    
    curr_results = f'{s.RESULTS_DIR}/results{i}/'
    if s.DEV and not os.path.exists(curr_results):
        os.mkdir(curr_results)
    start = time()

    # Do 100 gens in first GC sample at 50 and 100
    # 12 cells of each type at each time point, 4 total time points
    GC = [x for x in s.LOCATIONS if x.name == LocationName.GC][0]
    OTHER = [x for x in s.LOCATIONS if x.name == LocationName.OTHER][0]
    GC.sample_size = 12
    OTHER.sample_size = 20 
    GC.sample_times = [50, 100] 
    OTHER.sample_times = [50, 100]

    curr_time = 0
    clone_id = i+1
    naive = Cell(None, None, created_at=curr_time)
    root = Node(naive, clone_id=clone_id)
    airr = []
    TARGET_PAIR = TargetAminoPair(
        naive.heavy_chain.get_gapped_sequence(), 
        naive.light_chain.get_gapped_sequence(), 
        naive.heavy_chain.cdr3_length, 
        naive.light_chain.cdr3_length)
    TARGET_PAIR.mutate(s.TARGET_MUTATIONS_HEAVY, s.TARGET_MUTATIONS_LIGHT)

    sampled, pop_data, dev_df = simulate(clone_id, TARGET_PAIR, [root], root)

    sequenced = []
    to_sequence = s.RNG.choice([x for x in sampled if x.cell.location == LocationName.OTHER and x.generation == OTHER.sample_times[0]+1], size=12, replace=False)
    sequenced.extend(to_sequence)
    to_sequence = s.RNG.choice([x for x in sampled if x.cell.location == LocationName.OTHER and x.generation == OTHER.sample_times[1]+1], size=12, replace=False)
    sequenced.extend(to_sequence)
    to_sequence = s.RNG.choice([x for x in sampled if x.cell.location == LocationName.GC and x.generation == GC.sample_times[0]+1], size=12, replace=False)
    sequenced.extend(to_sequence)
    to_sequence = s.RNG.choice([x for x in sampled if x.cell.location == LocationName.GC and x.generation == GC.sample_times[1]+1], size=12, replace=False)
    sequenced.extend(to_sequence)

    to_reenter = s.RNG.choice([x for x in sampled if x not in sequenced and x.cell.location == LocationName.OTHER and x.generation == OTHER.sample_times[-1]+1], size=1, replace=False)

    first_end_time = max(GC.sample_times[-1], OTHER.sample_times[-1])
    new_gc_time = first_end_time + 1000

    def make_new_child(node, new_generation):
        child_cell = Cell(
            node.cell.heavy_chain.copy(),
            node.cell.light_chain.copy(),
            location=node.cell.location,
            created_at=new_generation)
        child_node = Node(
            child_cell,
            parent=node,
            heavy_mutations=0,
            light_mutations=0,
            generation=new_generation
            )
        child_cell.calculate_affinity(TARGET_PAIR)
        node.add_child(child_node)
        return child_node

    reentering_cells = [make_new_child(x, new_gc_time) for x in to_reenter]
    reentering_cells_gc = reentering_cells
    for x in reentering_cells_gc:
        x.cell.location = LocationName.GC
        x.cell.cell_type = CellType.DEFAULT

    # second gc reaction at time 1100
    TARGET_PAIR.mutate(s.TARGET_MUTATIONS_HEAVY, s.TARGET_MUTATIONS_LIGHT)

    # sample at 50 gens after first and second vaccine
    GC.sample_size = 12
    OTHER.sample_size = 12
    GC.sample_times = [new_gc_time+50, new_gc_time+100] 
    OTHER.sample_times = [new_gc_time+50, new_gc_time+100]

    # sample the recall reaction
    second_sampled, second_pop_data, _ = simulate(clone_id, TARGET_PAIR, reentering_cells_gc, root, new_gc_time)
    to_sequence = s.RNG.choice([x for x in second_sampled if x.cell.location == LocationName.OTHER and x.generation == OTHER.sample_times[0]+1], size=12, replace=False)
    sequenced.extend(to_sequence)
    to_sequence = s.RNG.choice([x for x in second_sampled if x.cell.location == LocationName.OTHER and x.generation == OTHER.sample_times[1]+1], size=12, replace=False)
    sequenced.extend(to_sequence)
    to_sequence = s.RNG.choice([x for x in second_sampled if x.cell.location == LocationName.GC and x.generation == GC.sample_times[0]+1], size=12, replace=False)
    sequenced.extend(to_sequence)
    to_sequence = s.RNG.choice([x for x in second_sampled if x.cell.location == LocationName.GC and x.generation == GC.sample_times[1]+1], size=12, replace=False)
    sequenced.extend(to_sequence)

    end_time = max(GC.sample_times[-1], OTHER.sample_times[-1])

    sampled = sequenced

    sampled_ids = [id(x.cell) for x in sampled]
    fasta_string = "".join([x.cell.as_fasta(x.sampled_time) for x in sampled])

    airr = [x for node in sampled for x in node.cell.as_AIRR(node.sampled_time)]
    airr = pd.DataFrame(airr)
    airr["sequence_id"] = airr["sequence_id"].apply(lambda x: f"{clone_id}_{x}")
    airr["cell_id"] = airr["cell_id"].apply(lambda x: f"{clone_id}_{x}")
    airr["clone_id"] = clone_id


    if s.MEMORY_SAVE:
        # in memory saving mode we don't keep the full tree
        newick = ""
        pruned_newick = ""
        pruned_time_tree = ""
        pruned = root
    elif s.KEEP_FULL_TREE:
        newick = f'{root.write_newick()};'
        pruned = root.prune_subtree(sampled_ids)
        pruned_newick = f'{pruned.write_newick()};'
        pruned_time_tree = f'{pruned.write_newick(time_tree=True)};'
    else:
        newick = ""
        pruned = root
        pruned_newick = f'{root.write_newick()};'
        pruned_time_tree = f'{root.write_newick(time_tree=True)};'

    pruned = pruned.prune_subtree(sampled_ids)
    simplified_tree = simplify_tree(pruned)
    simplified_tree_newick = f'{simplified_tree.write_newick()};'
    simplified_time_tree_newick = f'{simplified_tree.write_newick(time_tree=True)};'
    end = time()

    logger.debug("Time taken: %s", end - start)
    return {
        "airr": airr, 
        "fasta": fasta_string, 
        "full_tree": newick, 
        "pruned_tree": pruned_newick,
        "pruned_time_tree": pruned_time_tree, 
        "simplified_tree": simplified_tree_newick,
        "simplified_time_tree": simplified_time_tree_newick,
        "data": dev_df, 
        "clone_id": clone_id, 
        "pop_data": pop_data,
        "targets": {
            "clone_id": clone_id,
            "heavy": TARGET_PAIR.heavy.amino_acid_seq,
            "light": TARGET_PAIR.light.amino_acid_seq
            }
        }

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
        seed = 89243 #TyCHE
        ss = np.random.SeedSequence(seed)
    seeds = ss.spawn(args.n)
    print(f"Seed: {ss.entropy}")

    with tempfile.NamedTemporaryFile(mode="w") as tmpf:
        json.dump(s, tmpf, default=lambda o: o.encode(), indent=4)
        tmpf.flush()
        start = time()
        logger.info("Starting simulation")
        if args.processes > 1:
            with Pool(processes=args.processes) as pool:
                result = pool.starmap(
                    partial(do_simulation, filename=tmpf.name),
                    zip(range(args.n), seeds)
                    )
        else:
            result = []
            for i in range(args.n):
                result.append(
                    do_simulation(i, seeds[i], tmpf.name)
                    )

    process_results(result)

    end = time()
    logger.debug("Program finished! Total time taken: %s", end-start)

if __name__ == "__main__":
    main()