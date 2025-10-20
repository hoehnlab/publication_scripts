import os

import numpy as np
import pandas as pd

from simble.cell import Cell
from simble.parsing import get_parser, validate_and_process_args
from simble.settings import s
from simble.simble import logger, process_results, set_logger
from simble.simulation import run_simulation, simulate
from simble.target import TargetAminoPair
from simble.tree import Node, simplify_tree


def do_selection_simulation(i, naive):
    main_dir = s.RESULTS_DIR
    s.RESULTS_DIR = s.RESULTS_DIR + "/selection/"
    if not os.path.exists(s.RESULTS_DIR):
        os.mkdir(s.RESULTS_DIR)
    results = do_simulation(i, naive)
    process_results([results])
    s.RESULTS_DIR = main_dir

def do_neutral_simulation(i, naive):
    main_dir = s.RESULTS_DIR
    s.SELECTION = False
    s.RESULTS_DIR = s.RESULTS_DIR + "/neutral/"
    if not os.path.exists(s.RESULTS_DIR):
        os.mkdir(s.RESULTS_DIR)
    results = do_simulation(i, naive)
    process_results([results])
    s.RESULTS_DIR = main_dir

def do_uniform_simulation(i, naive):
    main_dir = s.RESULTS_DIR
    s.SELECTION = False
    s.UNIFORM = True
    s.RESULTS_DIR = s.RESULTS_DIR + "/uniform/"
    if not os.path.exists(s.RESULTS_DIR):
        os.mkdir(s.RESULTS_DIR)
    s.SEQUENCE_LENGTH = len(naive.heavy_chain.nucleotide_seq)
    results = run_simulation(i, s.RESULTS_DIR)
    process_results([results])
    s.RESULTS_DIR = main_dir


def do_simulation(i, naive):

    clone_id = i+1
    root = Node(naive, clone_id=clone_id)
    TARGET_PAIR = TargetAminoPair(
        naive.heavy_chain.get_gapped_sequence(),
        naive.light_chain.get_gapped_sequence(),
        naive.heavy_chain.cdr3_length,
        naive.light_chain.cdr3_length)
    TARGET_PAIR.mutate(s.TARGET_MUTATIONS_HEAVY, s.TARGET_MUTATIONS_LIGHT)

    sampled, pop_data, df = simulate(clone_id, TARGET_PAIR, [root], root)

    sampled_ids = [id(x.cell) for x in sampled]
    fasta_string = "".join([x.cell.as_fasta(x.sampled_time) for x in sampled])

    airr = [x for node in sampled for x in node.cell.as_AIRR(node.sampled_time)]
    airr = pd.DataFrame(airr)
    airr["sequence_id"] = airr["sequence_id"].apply(lambda x: f"{clone_id}_{x}")
    airr["cell_id"] = airr["cell_id"].apply(lambda x: f"{clone_id}_{x}")
    airr["clone_id"] = clone_id

    newick = ""
    pruned = root.prune_subtree(sampled_ids)
    pruned_newick = f'{pruned.write_newick()};'
    pruned_time_tree = f'{pruned.write_newick(time_tree=True)};'

    simplified_tree = simplify_tree(pruned)
    simplified_tree_newick = f'{simplified_tree.write_newick()};'
    simplified_time_tree_newick = f'{simplified_tree.write_newick(time_tree=True)};'

    return {
        "airr": airr, 
        "fasta": fasta_string, 
        "full_tree": newick, 
        "pruned_tree": pruned_newick,
        "pruned_time_tree": pruned_time_tree, 
        "simplified_tree": simplified_tree_newick,
        "simplified_time_tree": simplified_time_tree_newick,
        "data": df, 
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
        ss = np.random.SeedSequence()
    print(f"Seed: {ss.entropy}")

    s.UNIFORM = False
    s.SELECTION = True

    s._x_RNG = np.random.default_rng(seed) # pylint: disable=protected-access

    if not os.path.exists(s.RESULTS_DIR):
        os.mkdir(s.RESULTS_DIR)

    naive = Cell(None, None, created_at=0)


    logger.info("Starting simulations")
    do_selection_simulation(0, naive)
    do_neutral_simulation(0, naive)
    do_uniform_simulation(0, naive)


if __name__ == "__main__":
    main()