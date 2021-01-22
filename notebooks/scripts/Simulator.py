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
        min_Ne=5e5,
        max_Ne=5e6,
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
        
        # setup dataframe for results. Will be eventually filled with a
        # mix of integers and string types.
        self.data = pd.DataFrame(
            columns=[
                "spp_trees",
                "raxml_trees",
                "nsnps",
                "chr_trees_relaxed",
                "chr_trees_strict",
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
        #self.sample_sptree_rate_var()
        # self.batch_ipcoal()
        # self.batch_raxml()
        # self.batch_chronos()



    def sample_sptree_rate_var(self):
        """
        Apply rate variation to the edges of the species tree.
        """ 

        # sample arrays of random values to apply to edges
        samp_nes = self.rng.randint(
            low=self.min_n, 
            high=self.max_n, 
            size=(self.reps, self.sptree.ntips),
        )
        samp_gs = self.rng.randint(
            low=self.min_g, 
            high=self.max_g, 
            size=(self.reps, self.sptree.ntips),            
        )

        # iterate over sampled sets of values to create transformed trees.
        for idx in self.data.index:
                  
            # get a tree copy and random values
            tre = self.sptree.copy()
            vals_n = samp_nes[idx]
            vals_g = samp_gs[idx]

            # apply Nes to the tree
            tre = tre.set_node_values(
                "Ne", dict(zip(range(self.sptree.ntips), vals_n))
            )

            # apply Gs to the tree
            tre = tre.set_node_values(
                "g", dict(zip(range(self.sptree.ntips), vals_g))
            )
    
            # Divide edge lengths (absolute time) by generation time to 
            # convert tree to units of generations.
            tre = tre.set_node_values(
                "dist", 
                {i: j.dist / j.g for (i, j) in tre.idx_dict.items()}
            )
    
            # write tree to newick with Ne values on nodes
            self.data.loc[idx, "spp_trees"] = tre.write()

        logger.info(f"applied rate variation to {self.reps} spp. trees")

        #     reps.append(i)
        #     savefile = outdir + file_marker + '{0:03}'.format(counter) + ".tre"
        #     i.write(savefile)
            
        # # Save to instance variable and dataframe.
        # self.reps = reps
        # for i in range(self.ntrees):
        #     self.df.loc[i, "rate_trees"] = self.reps[i]
    
        # return "Saved " + str(self.ntrees) + " trees."
    



    def batch_ipcoal(self, outdir, file_marker):
        """
        Generate sequence data on tree variants.
        """       
        counter = 0
        seqs = []
        for i in self.reps:
        
            # Increment counter.
            counter += 1
        
            # Run ipcoal.
            model = ipcoal.Model(i, nsamples = 2, seed = self.seed) 
            model.sim_loci(1000, 100) # (loci, bp) 
            
            # get phy output name
            outname = file_marker + "_diploid" + '{0:03}'.format(counter) + ".phy"

            # Write a diploid phylip file.
            model.write_concat_to_phylip(
                outdir=outdir, 
                diploid=True,
                name=outname,
            )
            seqs.append(outdir + file_marker + "_diploid" + '{0:03}'.format(counter) + ".phy")
            
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