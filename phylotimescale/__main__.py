#!/usr/bin/env python

"""
Command line interface for weekdayfinder.
"""

import argparse
from phylotimescale import Simulator


def parse_command_line():
    "Parses arguments for Simulator."

    # Init parser and add arguments.
    parser = argparse.ArgumentParser()

    # Positional arguments.
    parser.add_argument(
        "tree",
        help = "Set true tree."
    )
    parser.add_argument(
        "reps",
        help = "Set number of replicates.",
        type = int
    )

    # Optional arguments.
    parser.add_argument(
        "--min_Ne",
        help = "Set lower bound of Ne variation.",
        type = int,
        nargs = "?",
        default = 100000
    )
    parser.add_argument(
        "--max_Ne",
        help = "Set upper bound of Ne variation.",
        type = int,
        nargs = "?",
        default = 1000000
    )
    parser.add_argument(
        "--min_g",
        help = "Set upper bound of g variation.",
        type = int,
        nargs = "?",
        default = 1
    )
    parser.add_argument(
        "--max_g",
        help = "Set upper bound of g variation.",
        type = int,
        nargs = "?",
        default = 1
    )
    parser.add_argument(
        "--outdir",
        help = "Set directory for output.  Default is current folder.",
        nargs = "?",
        const = "."
    )
    parser.add_argument(
        "--outprefix",
        help = "Set filename for output.",
        nargs = "?",
        const = "test"
    )
    parser.add_argument(
        "--ipcoal_kwargs",
        help = """Set params for ipcoal.  Enter as a dict.  Example:
        {
            "mut": 1e-8, 
            "recomb": 1e-9, 
            "nloci": 100, 
            "nsites": 1000,
        }
        """,
        nargs = "?",
        default = {
            "mut": 1e-8, 
            "recomb": 1e-9, 
            "nloci": 100, 
            "nsites": 1000,
        }
    )
    parser.add_argument(
        "--chronos_constraints",
        help = """Set params for chronos.  Enter as a dict.  Example:
        {
            ("r0", "r1"): (1e2, 1e4),
            ("r1", "r6"): (1e5, 1e5),        
        }
        """,
        nargs = "?",
        default = {
            ("r0", "r1"): (1e2, 1e4),
            ("r1", "r6"): (1e5, 1e5),        
        }
    )
    parser.add_argument(
        "--mb_params",
        help = """Set params for mrbayes.  Enter as a list of lists.  Example:
        [
        ["test1", "r0 r1 r2 r3", "uniform(1, 10)"],
        ["test2", "r4 r5", "uniform(1, 10)"],
        ["test3", "r6 r7 r8 r9", "uniform(1, 10)"]
        ],
        """,
        nargs = "?",
        default = [
        ["test1", "r0 r1 r2 r3", "uniform(1, 10)"],
        ["test2", "r4 r5", "uniform(1, 10)"],
        ["test3", "r6 r7 r8 r9", "uniform(1, 10)"]
        ],
    )
    # Look into action = append instead.  This can create a list: [constraint, tips, prior].  Then in Simulator, parse this
    # to build the constraint blocks in the Nexus file.  Use calib_string still, but with a list of 3-item lists or dicts.

    # The parse_args() call currently has no arguments, but I can add one to specify how to handle the append stuff.
    parser.add_argument(
        "--mb_treeagepr",
        help = "Set prior for mrbayes tree.",
        nargs = "?",
        const = "uniform(1, 10)"
    )
    parser.add_argument(
        "--seed",
        help = "Set seed for reproducibility.",
        type = int,
        nargs = "?",
        default = 12345
    )

    # Parse arguments.
    args = parser.parse_args()
    return args

def main():
    """
    Run main function on parsed arguments.
    """

    # Get arguments from command line as a dict-like object.
    args = parse_command_line()

    # Run Simulator if positional arguments are satisfied.
    if args.tree and args.reps:
        sim = Simulator(
        tree = args.tree,
        reps = args.reps,
        min_Ne = args.min_Ne,
        max_Ne = args.max_Ne,
        ipcoal_kwargs = args.ipcoal_kwargs,
        chronos_constraints = args.chronos_constraints,
        mb_params = args.mb_params,
        mb_treeagepr = args.mb_treeagepr,
        seed = args.seed
    )
        sim.run()

    # Save results to current directory.
        import os
        sim.data.to_csv(os.path.abspath("./sim.csv"))
    else:
        print("Tree and number of replicates required to run.")


if __name__ == "__main__":
    main()