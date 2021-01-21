#!/usr/bin/env python

"""
Describe this entire script here.
"""

import subprocess
import toytree
import numpy as np
import pandas as pd



class Analysis:
    """
    Describe this class here...
    """    
    def __init__(self, seed, sptree, df):
        
        # Store initial arguments.
        self.seed = seed
        self.sptree = sptree
        self.df = df
        
        # Objects to be created.
        self.chrtrees = []
        self.errors = []
        

    def batch_chronos(self, path_with_marker, min_ages, max_ages, tip1, tip2, lamb):
        """
        Run chronos on raxml trees.  Format min_ages, max_ages, tip1 
        and tip2 as tuples of the same length, formatted with
        double quotes enclosing parentheses, with constituent elements 
        in single quotes.
        Examples: min_ages = "('5000000', '10000000')"; tip1 = "('r0', 'r4')"
        """   
        np.random.seed(self.seed)
        chrtrees = []
        counter = 0
        for i in self.df["rax_trees"]:
        
            # Increment counter.
            counter += 1
        
            # Rstring with chronos information.
            rstring = f"""library(ape)
btree <- read.tree("{i}")
min_ages <- c{min_ages}
max_ages <- c{max_ages}
tip1 <- c{tip1}
tip2 <- c{tip2}
nodes <- c()
for (i in 1:length(tip1)) {{
    mrca <- getMRCA(btree, c(tip1[i], tip2[i]))
    nodes <- append(nodes, mrca)
}}
calib <- data.frame(node = nodes, age.min = as.numeric(min_ages), age.max = as.numeric(max_ages))
ctree <- chronos(btree, lambda = {lamb}, model = "relaxed", calibration = calib)
write.tree(ctree)"""
            
            # Write the R script to a file.
            with open(path_with_marker + '{0:03}'.format(counter) + ".R", 'w') as out:
                out.write(rstring)
            
            # byte = rstring.encode()
            # temp = tempfile.NamedTemporaryFile()
            # temp.write(byte)
            # temp.seek(0)
            # print(temp.read())
    
            cmd = ["Rscript", path_with_marker + '{0:03}'.format(counter) + ".R"] # Runs R script saved to path.
            out = subprocess.check_output(cmd).decode()
            results, tree = [i.strip().strip('"') for i in out.split("[1]")]

            # Add chronos tree to list.
            chrtrees.append(tree)
            
            # temp.close()
    
        # Save to instance variable and dataframe.
        self.chrtrees = chrtrees
        for i in self.df.index:
            self.df.loc[i, "chr_trees_relax"] = self.chrtrees[i]
    


    def calculate_error(self):
        "Calculate error between true tree and chronos trees."
        
        errors = []
        
        # Get edge lengths of true tree.
        true_edge_lengths = self.sptree.get_edge_values(feature = "dist")
        
        # For each chrtree, get edge lengths to subtract from true edge lengths.
        for i in self.chrtrees:
            chrtree = toytree.tree(i)
            chr_edge_lengths = chrtree.get_edge_values(feature = "dist")
            subtract_array = true_edge_lengths - chr_edge_lengths
            
            # Square each element in the array.
            squared_array = np.square(subtract_array)
            
            # Sum all elements in the array (sum of squares).
            sum_squares = np.sum(squared_array)
            
            # Add error to list.
            errors.append(sum_squares)
            
        # Save to instance variable and dataframe.
        self.errors = errors
        for i in self.df.index:
            self.df.loc[i, "error"] = self.errors[i]
    


    def batch_mrbayes(self):
        "call mrbayes (maybe using ipa) to infer trees"
        pass