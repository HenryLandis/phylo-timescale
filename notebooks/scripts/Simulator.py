#!/usr/bin/env python

"""
Describe this entire script here.
"""


import os
import toytree
import ipcoal
import numpy as np
import pandas as pd
from loguru import logger



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

    """    
    def __init__(
        self, 
        sptree, 
        reps,
        min_Ne=100000,
        max_Ne=1000000,
        min_Gen=1,
        max_Gen=1,
        outdir=".",
        outprefix="test",
        seed=None,
        ):
        
        # Store initial arguments.
        self.sptree = toytree.tree(sptree)
        self.reps = int(reps)

        # store ratevar params
        self.min_n = min_Ne
        self.max_n = max_Ne
        self.min_g = min_Gen
        self.max_g = max_Gen

        # path to output files
        self.outdir = os.path.realpath(os.path.expanduser(outdir))
        self.prefix = os.path.join(self.outdir, outprefix)

        # random number generator
        self.rng = np.random.default_rng(seed=seed)
        
        # the Ne and Gentime values used in simulations
        # sample arrays of random values to apply to edges
        self.samp_ns = self.rng.integers(
            low=self.min_n, 
            high=(self.max_n + 1 if self.max_n == self.max_n else self.min_n),
            size=(self.reps, self.sptree.nnodes),
        )

        # if no variation in G then do not sample  
        self.samp_gs = self.rng.uniform(
            low=self.min_g, 
            high=self.max_g, 
            size=(self.reps, self.sptree.nnodes),            
        )

        # setup dataframe for results. Will be eventually filled with a
        # mix of integers and string types.
        self.data = pd.DataFrame(
            columns=[
                "spp_tree",
                "seqalign",
                "nsnps",
                "raxml_tree",
                "chronos_relaxed",
                "chronos_strict",
                "error",
            ],
            index=range(self.reps),
            data=0,
            dtype=int,
        )

     
        
    # , outdir, file_marker, path, min_ages, max_ages, tips, lamb):
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
        # self.batch_raxml()
        # self.batch_chronos()



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
            self.data.loc[idx, "spp_trees"] = tre.write(tree_format=0)

        # report to user
        logger.info(f"applied rate variation to {self.reps} spp. trees")



    def simulate_geneal_and_seqs(self, outdir, file_marker):
        """
        Setup ipcoal simualtion using sptree in units of generations
        and apply Ne values from the .samp_ns array. Simulate 
        genealogies and sequence data on each tree.
        """       
        for idx in self.data.index:
            
            # load the transformed sptree
            tre = self.data.at[idx, "spp_trees"]

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
                mut=1e-8,
                recomb=1e-9,
            )
            model.sim_loci(1000, 100) # (loci, bp) 
            
            # Write a diploid phylip file.
            model.write_concat_to_phylip(
                outdir=self.outdir, 
                diploid=True,
                name=self.prefix,
            )

            # store genealogy to species tree
            self.data.loc[idx, ""]
            
        # Save the list of concatenated phylip files to instance variable.   
        self.seqs = seqs


            
    def batch_raxml(self, outdir, file_marker, outdir_for_chronos):
        '''
        Infer raxml trees from sequence data.
        '''
        
        np.random.seed(self.seed)
        raxtrees = []
        counter = 0
        for i in self.seqs:
        
            # Increment counter.
            counter += 1
        
            # Define and run raxml object.
            rax = ipa.raxml(
                name = file_marker + '{0:03}'.format(counter),
                data = i,
                workdir = outdir,
                N = 100,
                T = 10 # Core assignment appropriate for pinky or other external server.
            )
            rax.run(force = True, block = False)
            
            # Take the raxml result and save as a newick string (required format for chronos).
            rax_result = outdir + "RAxML_bipartitions." + file_marker + '{0:03}'.format(counter)
            rax_toytree = toytree.tree(rax_result).root(["r7", "r8", "r9", "r10", "r11"]) # Rooting midpoint of true tree.
            rax_toytree.write(outdir_for_chronos + file_marker + '{0:03}'.format(counter) + ".tre")
            
            # Add raxtree to list.
            # raxtrees.append(outdir_for_chronos + file_marker + '{0:03}'.format(counter) + ".tre")
            raxtrees.append(rax_toytree)
            
        # Save to instance variable and dataframe.
        self.raxtrees = raxtrees
        for i in self.df.index:
            self.df.loc[i, "rax_trees"] = self.raxtrees[i]