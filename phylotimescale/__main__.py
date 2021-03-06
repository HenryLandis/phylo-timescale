#!/usr/bin/env python

"""
Command line interface for weekdayfinder.
"""

import argparse
from phylotimescale import Simulator


def parse_command_line():
    """
    Parses arguments for Simulator.
    """

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
        help = "Set params for ipcoal.  Use quotes for args with spaces.",
        metavar = "key=value",
        nargs = "*",
        default = {
            "mut": 1e-8, 
            "recomb": 1e-9, 
            "nloci": 100, 
            "nsites": 1000,
        }
    )
    parser.add_argument(
        "--chronos_constraints",
        help = "Set params for chronos.  Use quotes for args with spaces.",
        nargs = 4,
        #action = "append",
        default = None
    )
    parser.add_argument(
        "--mb_params",
        help = "Set params for mrbayes.",
        nargs = 3,
        action = "append",
        default = None
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

def parse_pair_for_ipcoal(s):
    """
    Parses a key-value pair, separated by '='.
    """

    items = s.split('=')
    key = items[0].strip()
    if len(items) > 1:
        value = '='.join(items[1:])
    return (key, value)


def parse_ipcoal_to_dict(items):
    """
    Parse a series of key-value pairs and return a dictionary of strings : flaots for ipcoal.
    """

    d = {}
    if items:
        for item in items:
            key, value = parse_pair_for_ipcoal(item)
            d[key] = float(value)
    return d


def parse_lists_for_chronos(s):
    """
    Parse a list of lists and return a dictionary of tuples : tuples for chronos.  Each individual call of --chronos_constraints
    is expected to include exactly two tips and exactly two age values, of the form '(tip1)' '(min.age)' '(tip2)' '(max.age)'.  
    Example:

    ("r0", "r19", "r32") : (8460000, 12310000, 10380000), ("r8", "r25", "r35") : (18460000, 22310000, 20380000)
    """
    d = {}
    print(s)
    tip1 = ()
    min_age = ()
    tip2 = ()
    max_age = ()
    for i in s[0].split(" "):
        tip1 = tip1 + (i ,)
    for i in s[1].split(" "):
        min_age = min_age + (i ,)
    for i in s[2].split(" "):
        tip2 = tip2 + (i ,)
    for i in s[3].split(" "):
        max_age = max_age + (i ,)
    d[tip1] = min_age
    d[tip2] = max_age
    print(d) # This dict is formatted exactly how I want for the chronos script.  I need to edit Simulator to build the arguments
    # directly from this setup.  Use d.get() to retrieve values.  Write a small function to enter values and get keys.

    # Alternative: do list of lists like for mrbayes.  Advantage is to input clades in a sensible manner.  Also requires Simulator
    # edits.
    print([str(i[0]) for i in d.keys()])
    print([str(i[1]) for i in d.keys()])
    print([str(i[0]) for i in d.values()])
    print([str(i[1]) for i in d.values()])
    return d


def main():
    """
    Run main function on parsed arguments.
    """

    # Get arguments from command line as a dict-like object.
    args = parse_command_line()

    # Build dicts for ipcoal and chronos params.
    ipcoal_dict = parse_ipcoal_to_dict(args.ipcoal_kwargs)
    # chronos_dict = parse_lists_for_chronos(args.chronos_constraints)

    # Run Simulator if positional arguments are satisfied.
    if args.tree and args.reps:
        sim = Simulator(
        tree = args.tree,
        reps = args.reps,
        min_Ne = args.min_Ne,
        max_Ne = args.max_Ne,
        ipcoal_kwargs = ipcoal_dict,
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