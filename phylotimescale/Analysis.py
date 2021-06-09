#!/usr/bin/env python

"""
Describe this entire script here.
"""

import toytree
import numpy as np
import pandas as pd
from statistics import mean

# makes pandas DF with wide columns look nicer
pd.set_option("max_colwidth", 14)

class Analysis:
    """
    Calculate simple and relative errors of trees from Simulator object compared to a given true tree.
    """    
    def __init__(
        self,
        tree,
        df,
        seed=None
        ):

        # Store the species tree.
        self.sptree = toytree.tree(tree)
        
        # Store copy of pd dataframe from Simulator object.
        self.data = df.copy()
    


    def run(self, calc = "mean"):
        """
        Run simple and relative error.  
        """
        self.simple_error(calc = calc)
        self.relative_error(calc = calc)



    def simple_error(self, calc = "mean"):
        "Calculate simple error between true tree and chronos or mrbayes output."
        
        # Get edge lengths of true tree.
        true_edge_lengths = self.sptree.get_edge_values(feature="height")
        
        for idx in self.data.index:

            # Calculate error for three tree types.
            for model in ["chronos_correlated", "chronos_relaxed", "mrbayes_tree"]:

                # Get edge lengths to subtract from true edge lengths.
                dtree = toytree.tree(self.data.at[idx, model])
                dtree_edge_lengths = dtree.get_edge_values(feature="height")

                # Average of absolute difference between node neights.
                if calc =="mean":

                # mrbayes: correct for output and rotations, then get average per-branch deviation from true node height.
                    # if model == "mrbayes_tree":
                        # dtree_edge_lengths = dtree_edge_lengths * 1e6
                    dists = []
                    for nidx in self.sptree.idx_dict:
                        tips = self.sptree.get_tip_labels(nidx)
                        dtree_nidx = dtree.get_mrca_idx_from_tip_labels(tips)
                        dtree_node = dtree.idx_dict[dtree_nidx]
                        true_node = self.sptree.idx_dict[nidx]
                        # print(true_node.height, dtree_node.height, dtree_node.height * 1e6)
                        if true_node.is_leaf() == False:
                            diff = abs(true_node.height - (dtree_node.height * 1e6)) / true_node.height
                        dists.append(diff)
                        # print(dists)
                    if model == "chronos_correlated":
                        self.data.loc[idx, "chc_simple_error"] = mean(dists)
                    elif model == "chronos_relaxed":
                        self.data.loc[idx, "chr_simple_error"] = mean(dists)
                    elif model == "mrbayes_tree":
                        self.data.loc[idx, "mb_simple_error"] = mean(dists)

                # chronos: calculate array.
                    #else:
                        #subtract_array = abs(true_edge_lengths - dtree_edge_lengths) / true_edge_lengths
                        #if model == "chronos_correlated":
                            #self.data.loc[idx, "chc_simple_error"] = subtract_array.mean()
                        #elif model == "chronos_relaxed":
                            #self.data.loc[idx, "chr_simple_error"] = subtract_array.mean()

                # Sum of squares of absolute differences between node neights.
                elif calc == "ss":
                    if model == "mrbayes_tree":
                        dists = []
                        for nidx in self.sptree.idx_dict:
                            tips = self.sptree.get_tip_labels(nidx)
                            dtree_nidx = dtree.get_mrca_idx_from_tip_labels(tips)
                            dtree_node = dtree.idx_dict[dtree_nidx]
                            true_node = self.sptree.idx_dict[nidx]
                            # print(true_node.height, dtree_node.height, dtree_node.height * 1e6)
                            if true_node.height == 0:
                                sqdiff = np.square(abs(1 - (dtree_node.height * 1e6)))
                            else:
                                sqdiff = np.square(abs(true_node.height - (dtree_node.height * 1e6)) / true_node.height)
                            dists.append(sqdiff)
                        # print(dists)
                        self.data.loc[idx, "mb_simple_error"] = sum(dists)
                    else:
                        subtract_array = abs(true_edge_lengths - dtree_edge_lengths) / true_edge_lengths
                        squared_array = np.square(subtract_array)
                        sum_squares = np.sum(squared_array)
                        if model == "chronos_correlated":
                            self.data.loc[idx, "chc_simple_error"] = sum_squares
                        elif model == "chronos_relaxed":
                            self.data.loc[idx, "chr_simple_error"] = sum_squares

                # Raise error for unexpected string.
                else:
                    print("Please enter a valid calc type.")



    def relative_error(self, calc = "mean"):
        "Calculate relative error between true tree and chronos or mrbayes output."

        # Get edge lengths of true tree.
        true_edge_lengths = self.sptree.get_edge_values(feature="height")
        
        for idx in self.data.index:

            # Calculate error for three tree types.
            for model in ["chronos_correlated", "chronos_relaxed", "mrbayes_tree"]:

                # Scale tree to match height of true tree, then get edge lengths to subtract from true edge lengths.
                dtree = toytree.tree(self.data.at[idx, model])
                dtree = dtree.mod.node_scale_root_height(treeheight=list(self.sptree.get_feature_dict("height").keys())[0])
                dtree_edge_lengths = dtree.get_edge_values(feature="height")

               # Average of absolute difference between node neights.
                if calc =="mean":

                # mrbayes: correct for output and rotations, then get average per-branch deviation from true node height.
                    #if model == "mrbayes_tree":
                        # dtree_edge_lengths = dtree_edge_lengths * 1e6
                    dists = []
                    for nidx in self.sptree.idx_dict:
                        tips = self.sptree.get_tip_labels(nidx)
                        dtree_nidx = dtree.get_mrca_idx_from_tip_labels(tips)
                        dtree_node = dtree.idx_dict[dtree_nidx]
                        true_node = self.sptree.idx_dict[nidx]
                        # print(true_node.height, dtree_node.height, dtree_node.height * 1e6)
                        if true_node.is_leaf() == False:
                            diff = abs(true_node.height - dtree_node.height) / true_node.height
                        dists.append(diff)
                        # print(dists)
                    if model == "chronos_correlated":
                        self.data.loc[idx, "chc_relative_error"] = mean(dists)
                    elif model == "chronos_relaxed":
                        self.data.loc[idx, "chr_relative_error"] = mean(dists)
                    elif model == "mrbayes_tree":
                        self.data.loc[idx, "mb_relative_error"] = mean(dists)
                        # print(list(dtree.get_feature_dict("height").keys())[0])

                # chronos: calculate array.
                    #else:
                        #subtract_array = abs(true_edge_lengths - dtree_edge_lengths) / true_edge_lengths
                        #if model == "chronos_correlated":
                            #self.data.loc[idx, "chc_relative_error"] = subtract_array.mean()
                        #elif model == "chronos_relaxed":
                            #self.data.loc[idx, "chr_relative_error"] = subtract_array.mean()

                # Sum of squares of absolute differences between node neights.
                elif calc == "ss":
                    if model == "mrbayes_tree":
                        dists = []
                        for nidx in self.sptree.idx_dict:
                            tips = self.sptree.get_tip_labels(nidx)
                            dtree_nidx = dtree.get_mrca_idx_from_tip_labels(tips)
                            dtree_node = dtree.idx_dict[dtree_nidx]
                            true_node = self.sptree.idx_dict[nidx]
                            # print(true_node.height, dtree_node.height, dtree_node.height * 1e6)
                            if true_node.height == 0:
                                sqdiff = np.square(abs(1 - dtree_node.height))
                            else:
                                sqdiff = np.square(abs(true_node.height - dtree_node.height) / true_node.height)
                            dists.append(sqdiff)
                        # print(dists)
                        self.data.loc[idx, "mb_relative_error"] = sum(dists)
                    else:
                        subtract_array = abs(true_edge_lengths - dtree_edge_lengths) / true_edge_lengths
                        squared_array = np.square(subtract_array)
                        sum_squares = np.sum(squared_array)
                        if model == "chronos_correlated":
                            self.data.loc[idx, "chc_relative_error"] = sum_squares
                        elif model == "chronos_relaxed":
                            self.data.loc[idx, "chr_relative_error"] = sum_squares

                # Raise error for unexpected string.
                else:
                    print("Please enter a valid calc type.")









if __name__ == "__main__":

    # Run this test with ctrl+shift+B in editor

    from phylotimescale import Simulator

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
        mb_params=[
        ["test1", "r0 r1 r2 r3", "uniform(1, 10)"],
        ["test2", "r4 r5", "uniform(1, 10)"],
        ["test3", "r6 r7 r8 r9", "uniform(1, 10)"]
        ],
        mb_treeagepr="uniform(1, 10)"
    )
    sim.run()

    ana=Analysis(
        tree=sim.sptree,
        df=sim.data
        )
    ana.run()

    # show results
    print(ana.data.T)
    chtree = toytree.tree(ana.data.at[0, 'chronos_relaxed'])
    print(chtree)
    mbtree = toytree.tree(ana.data.at[0, 'mrbayes_tree'])
    print(mbtree)

    # save results to csv
    ana.data.to_csv("~/phylo-timescale/ana.csv")