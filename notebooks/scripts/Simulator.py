class Simulator:
    
    def __init__(self, seed, sptree, ntrees):
        
        # Store initial arguments.
        self.seed = seed # Seed for all instances of the class.  
        self.sptree = sptree # True species tree.
        self.ntrees = ntrees # Number of replicate trees.
        # self.ntips = ntips  Number of tips.
        
        # Objects to be created.
        self.reps = []
        self.seqs = []
        self.raxtrees = []
        
        # Organize results into a dataframe.
        self.df = pd.DataFrame({
            
            # Trees with rate variation applied to edges.
            "rate_trees": ["" for i in range(ntrees)],
            
            # Gene trees with variable edge lengths inferred by RAxML.
            "rax_trees": ["" for i in range(ntrees)],
            
            # Number of SNPs from RAxML trees.
            "nsnps": 0,

            # Ultrametric trees inferred from chronos under relaxed model.
            "chr_trees_relax": ["" for i in range(ntrees)],
            
            # # Ultrametric trees inferred from chronos under strict model.
            "chr_trees_strict": ["" for i in range(ntrees)],
            
            # Error between tree from chronos and true tree.
            "error" : ["" for i in range(ntrees)]
        })
        
    # Execute all the major functions in a row.
    def run(self, outdir, file_marker, path, min_ages, max_ages, tips, lamb):
        self.batch_treedata()
        self.batch_ipcoal()
        self.batch_raxml()
        self.batch_chronos()

    def batch_treedata(self, outdir, file_marker):
        '''
        Apply rate variation to the edges of the species tree.
        '''
        
        np.random.seed(self.seed)
        self.reps = [self.sptree for i in range(self.ntrees)]
        reps = []
        counter = 0
        for i in self.reps:
        
            # Increment counter.
            counter +=1
    
            # Set Ne values from an interval onto tree.
            dict_ne = {j.name : np.random.randint(5e5, 5e6) for j in i.get_feature_dict()}
            i = i.set_node_values("Ne", dict_ne)
    
            # Set g values from a normal distribution onto tree.
            dict_g = {j.name : np.random.normal(1, 0.2) for j in i.get_feature_dict()}
            i = i.set_node_values("g", dict_g)
    
            # Divide edge lengths (absolute time) by generation time to convert tree to units of generations.
            i = i.set_node_values(
            "dist",
            {j.name: j.dist / j.g for j in i.get_feature_dict()}
            )
    
            # Save tree to a list and as a separate newick file.
            reps.append(i)
            savefile = outdir + file_marker + '{0:03}'.format(counter) + ".tre"
            i.write(savefile)
            
        # Save to instance variable and dataframe.
        self.reps = reps
        for i in self.df.index:
            self.df.loc[i, "rate_trees"] = self.reps[i]
    
        return "Saved " + str(self.ntrees) + " trees."
    
    def batch_ipcoal(self, outdir, file_marker):
        '''
        Generate sequence data on tree variants.
        '''
        
        counter = 0
        seqs = []
        for i in self.reps:
        
            # Increment counter.
            counter += 1
        
            # Run ipcoal.
            model = ipcoal.Model(i, nsamples = 2, seed = self.seed) 
            model.sim_loci(1000, 100) # (loci, bp) 
            
            # Write a diploid phylip file.
            file = model.write_concat_to_phylip(outdir = outdir, diploid = True,
                                      name = file_marker + "_diploid" + '{0:03}'.format(counter) + ".phy")
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