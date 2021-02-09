#!/usr/bin/env python

"""
Describe this entire script here.
"""

import os
import tempfile
import subprocess
import toytree
import ipcoal
import numpy as np
import pandas as pd
import ipyrad.analysis as ipa


# makes pandas DF with wide columns look nicer
pd.set_option("max_colwidth", 14)


# Rscript that can accept arbitrary number of calibrated nodes for chronos.
# entered as {(tip1, tip2): (min_age, max_age)} to apply the 
# age constraint to the mrca of the two tips.
RSTRING = """
library(ape)

# load the Python variables
btree <- read.tree(text="{raxml_tree}")
min_ages <- c({min_ages})
max_ages <- c({max_ages})
tips1 <- c({tips1})
tips2 <- c({tips2})
lamb <- {lamb}
model <- {model}

# run the R code 
nodes <- c()
for (i in 1:length(tips1)) {{
    mrca <- getMRCA(btree, c(tips1[i], tips2[i]))
    nodes <- append(nodes, mrca)
}}

calib <- data.frame(node=nodes, age.min=min_ages, age.max=max_ages)

ctree <- chronos(
    btree, 
    lambda=lamb,
    model=model,
    calibration=calib,
)
write.tree(ctree)
"""

# strings to define calibration settings in mrbayes.
treeage_string = "prset treeagepr={};"
calib_string = '''constraint {0}={1};
  calibrate {0}={2};'''
topology_string = '''prset topologypr=constraints({});
  prset nodeagepr=calibrated;'''



class Simulator:
    """
    Describe this class here...

    ... stores results in a Pandas DataFrame as an attribute
    of the Simulator class in .data.  

    Parameters
    -----------
    sptree (str):
        Enter a newick string or toytree of the true species tree, which
        is ultrametric with edge lengths in units of time (e.g., years),
        but not necessarily in generations (you can set the generation
        times for all species using a range set w/ min_Gen and max_Gen).
    reps (int):
        The number of simulations replicates to run.
    min_N (int):
        ...
    max_N (int):
        ...

    ipcoal_kwargs (dict):
        A dictionary of arguments to the ipcoal.Model object init. 
        Examples include {'mut': 1e-8, 'recomb': 1e-9}.
    """
    def __init__(
        self,
        tree, 
        reps,
        min_Ne=100000,
        max_Ne=1000000,
        min_Gen=1,
        max_Gen=1,
        outdir=".",
        outprefix="test",
        ipcoal_kwargs={
            "mut": 1e-8, 
            "recomb": 1e-9, 
            "nloci": 100, 
            "nsites": 1000,
        },
        chronos_constraints={},
        mb_names=[],
        mb_tips=[],
        mb_priors=[],
        mb_treeagepr=None,
        seed=None,
        ):
        
        # Store the species tree and its rooting
        self.sptree = toytree.tree(tree)
        self.root = self.sptree.get_tip_labels(self.sptree.treenode.children[0].idx)

        # number of sampled replicates to simulate
        self.reps = int(reps)

        # store ratevar params
        self.min_n = min_Ne
        self.max_n = max_Ne
        self.min_g = min_Gen
        self.max_g = max_Gen

        # store arguments to raxml and chronos, pull nsites and nloci out.
        self.ipcoal_kwargs = ipcoal_kwargs
        self.chronos_constraints = chronos_constraints
        if "nloci" in self.ipcoal_kwargs:
            self.nloci = self.ipcoal_kwargs.pop("nloci")
        else:
            self.nloci = 100
        if "nsites" in self.ipcoal_kwargs:
            self.nsites = self.ipcoal_kwargs.pop("nsites")
        else:
            self.nsites = 1000

        # path to output files
        self.outdir = os.path.realpath(os.path.expanduser(outdir))
        self.prefix = os.path.join(self.outdir, outprefix)

        # random number generator
        self.rng = np.random.default_rng(seed=seed)
        
        # the Ne and Gentime values used in simulations
        # sample arrays of random values to apply to edges
        self.samp_ns = self.rng.integers(
            low=self.min_n, 
            high=(self.min_n + 1 if self.min_n == self.max_n else self.max_n),
            size=(self.reps, self.sptree.nnodes),
        )

        # if no variation in G then do not sample  
        self.samp_gs = self.rng.uniform(
            low=self.min_g, 
            high=self.max_g, 
            size=(self.reps, self.sptree.nnodes),            
        )

        # store mrbayes params
        self.mb_names = mb_names
        self.mb_tips = mb_tips
        self.mb_priors = mb_priors
        self.mb_treeagepr = mb_treeagepr

        # setup dataframe for results. Will be eventually filled with a
        # mix of integers and string types.
        self.data = pd.DataFrame(
            columns=[
                "spp_tree",
                "phy_seqpath",
                "nex_seqpath",
                "nsnps",
                "raxml_tree",
                "chronos_correlated",
                "chc_simple_error",
                "chc_relative_error",
                "chronos_relaxed",
                "chr_simple error",
                "chr_relative_error",
                "mrbayes_tree",
                "mb_simple_error",
                "mb_relative_error"
            ],
            index=range(self.reps),
            data=0,
            dtype=int,
        )



    def run(self):
        """
        Runs a series of functions to:
            1. sample sets of rate-var parameters to apply to spp. tree
            2. simulate genealogies & seqs on transformed spp. trees.
            3. infer gene trees on simulated sequences.
            4. infer chronogram from raxml gene trees.           
        """
        self.transform_sptree_to_gen_units()
        self.simulate_geneal_and_seqs()
        self.batch_raxml()
        self.batch_chronos()
        self.batch_mrbayes()



    def transform_sptree_to_gen_units(self):
        """
        Apply rate variation to the edges of the species tree by 
        transforming edge lengths into units of generations.
        """ 
        # iterate over sampled sets of values to create transformed trees.
        for idx in self.data.index:

            # get a tree copy and random values
            tre = self.sptree.copy()

            # Divide edge lengths (absolute time) by generation time to 
            # convert tree to units of number of generations. Tree is 
            # now likely to be non-ultrametric.
            gdict = dict(zip(range(self.sptree.nnodes), self.samp_gs[idx]))
            tre = tre.set_node_values("g", values=gdict)
            tre = tre.set_node_values(
                "dist", 
                {i: j.dist / j.g for (i, j) in tre.idx_dict.items()}
            )

            # save tree to newick with Ne values as "names" in format=1
            self.data.loc[idx, "spp_tree"] = tre.write(tree_format=0)

        # report to user
        print(f"applied rate variation to {self.reps} spp. trees")



    def simulate_geneal_and_seqs(self):
        """
        Setup ipcoal simualtion using sptree in units of generations
        and apply Ne values from the .samp_ns array. Simulate 
        genealogies and sequence data on each tree.
        """
        for idx in self.data.index:
            
            # load the transformed sptree
            tre = toytree.tree(self.data.at[idx, "spp_tree"])

            # set Ne values on the tree, which ipcoal expects
            tre = tre.set_node_values(
                "Ne", 
                dict(zip(range(tre.nnodes), self.samp_ns[idx])),
            )

            # simulate genealogies on this species tree
            model = ipcoal.Model(
                tree=tre, 
                nsamples=2, 
                seed=self.rng.integers(0, 1e9),
                **self.ipcoal_kwargs,
            )
            model.sim_loci(self.nloci, self.nsites)
            
            # Write a diploid phylip file.
            model.write_concat_to_phylip(
                name=self.prefix + "_{}".format(idx),
                outdir=self.outdir,
                diploid=True,
            )

            # Write a diploid nexus file.
            model.write_concat_to_nexus(
                name=self.prefix + "_{}".format(idx),
                outdir=self.outdir,
                diploid=True,
            )

            # store the number of snps
            self.data.loc[idx, "nsnps"] = model.df.nsnps.sum()

            # store the path to the sequence alignment
            self.data.loc[idx, "phy_seqpath"] = os.path.join(
                self.outdir, self.prefix + "_{}.phy".format(idx)
            )

            self.data.loc[idx, "nex_seqpath"] = os.path.join(
                self.outdir, self.prefix + "_{}.nex".format(idx)
            )
        print("simulated sequences on {} species trees.".format(self.reps))


           
    def batch_raxml(self):
        """
        Infer raxml trees from sequence data.
        """       
        for idx in self.data.index:
        
            # Define and run raxml object.
            rax = ipa.raxml(
                name="tmp",
                data=self.data.at[idx, "phy_seqpath"],
                workdir=tempfile.gettempdir(),
                N=100,
                T=8,
            )
            rax.run(force=True, block=False, quiet=True)
            
            # save the newick string to file
            raxtree = toytree.tree(rax.trees.bipartitions)
            raxtree = raxtree.root(self.root)
            self.data.loc[idx, "raxml_tree"] = raxtree.write(tree_format=0)

            # until/unless we add mrbayes, we can remove sequence files here
            # os.remove(self.data.at[idx, "seqpath"])


    def batch_chronos(self):
        """
        Run chronos on raxml trees to infer an ultrametric tree with 
        both rates=gamma and rates=correlated.
        """
        for idx in self.data.index:

            # build the chronos argument dict
            chronos_kwargs = {
                "raxml_tree": self.data.at[idx, "raxml_tree"],
                "min_ages": ",".join(
                    [str(i[0]) for i in self.chronos_constraints.values()]
                ),
                "max_ages": ",".join(
                    [str(i[1]) for i in self.chronos_constraints.values()]
                ),
                "tips1": ",".join(
                    ["'{}'".format(i[0]) for i in self.chronos_constraints.keys()]
                ),
                "tips2": ",".join(
                    ["'{}'".format(i[1]) for i in self.chronos_constraints.keys()]
                ),
                "lamb": str(1.0),
            }

            # run on both models:
            for model in ["relaxed", "correlated"]:

                # store the model as an arg to chronos
                chronos_kwargs["model"] = "'{}'".format(model)

                # build Rstring with chronos information.
                rstring = RSTRING.format(**chronos_kwargs)

                # Write the R script to a file.
                tmp_rscript = os.path.join(tempfile.gettempdir(), "tmp.R")
                with open(tmp_rscript, 'w') as out:
                    out.write(rstring)
            
                # call the Rscript file to get relaxed chronos tree
                out = subprocess.run(
                    ["Rscript", tmp_rscript], 
                    check=True,
                    stdout=subprocess.PIPE,
                )

                # load the output
                ctre = out.stdout.decode().strip().split("[1] ")[-1].strip('"')
                key = "chronos_{}".format(model)
                self.data.loc[idx, key] = ctre

    def batch_mrbayes(self):
        '''Run mrbayes on sequence data from ipcoal to infer trees.'''

        for idx in self.data.index:

            # Define and run mrbayes object.
            mb = ipa.mrbayes(
                data = self.data.at[idx, "nex_seqpath"],
                name = "tmp",
                workdir = tempfile.gettempdir(),
                clock_model = 2,
                constraints = self.sptree,
                ngen = int(1e6),
                nruns = 2,
            )

            # Add priors for tree age & clade names, plus topology fixing.
            with open(mb.nexus, 'r+') as fd:
                contents = fd.readlines()
                fd.seek(0)
                contents.insert(24, "  " + treeage_string.format(self.mb_treeagepr) + "\n\n")
                contents.insert(27, "  " + topology_string.format(", ".join(self.mb_names)) + "\n")
                contents.insert(29, '''  startvals Tau = fixedtree;
            propset ExtSprClock(Tau,V)$prob=0;
            propset NNIClock(Tau,V)$prob=0;
            propset ParsSPRClock(Tau,V)$prob=0;\n\n''')
                fd.seek(0)

                # Add calibrations for individual clades.
                for i in range(len(self.mb_names)):
                    if i == len(self.mb_names) - 1:
                        contents.insert((25 + i), "  " + calib_string.format(self.mb_names[i], self.mb_tips[i], self.mb_priors[i]))
                        fd.seek(0)
                    else: 
                        contents.insert((25 + i), 
                        "  " + calib_string.format(self.mb_names[i], self.mb_tips[i], self.mb_priors[i]) + "\n\n")
                        fd.seek(0)

                # Remove unnecessary line.
                for line in contents:
                    if line.startswith("  prset topologypr=fixed(fixedtree);"):
                        pass
                    else:
                        fd.write(line)

                # if above fails, split read/write as two with-loops, use contents variable to transfer

            # Run mrbayes object.
            subprocess.run(
                    ["mb", mb.nexus], 
                    check=True,
                    stdout=subprocess.PIPE,
                )

            # Store output.
            self.data.loc[idx, "mrbayes_tree"] = toytree.tree(mb.trees.constre, tree_format=10)









if __name__ == "__main__":

    # Run this test with ctrl+shift+B in editor

    sptree = toytree.rtree.unittree(ntips=10, treeheight=1e6, seed=123)
    print(sptree)

    sim = Simulator(
        tree=sptree, 
        reps=1, 
        min_Ne=1e4,
        max_Ne=1e4,
        seed=123,
        ipcoal_kwargs={
            "mut": 1e-8, 
            "recomb": 1e-9, 
            "nloci": 20, 
            "nsites": 1000,
        },
        chronos_constraints={
            ("r0", "r1"): (1e2, 1e4),
            ("r1", "r6"): (1e5, 1e5),        
        },
        mb_names=["test1", "test2", "test3"],
        mb_tips=["r0 r1 r2 r3", "r4 r5", "r6 r7 r8 r9"],
        mb_priors=["uniform(1, 100)", "uniform(1, 100)", "uniform(1, 100)"],
        mb_treeagepr="uniform(1, 100)"
    )
    sim.run()

    # show results
    print(sim.data.T)
    chtree = toytree.tree(sim.data.at[0, 'chronos_relaxed'])
    print(chtree)
    mbtree = toytree.tree(sim.data.at[0, 'mrbayes_tree'])
    print(mbtree)

    # save results to csv
    sim.data.to_csv("~/phylo-timescale/sim.csv")