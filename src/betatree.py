#!/ebio/ag-neher/share/programs/EPD/bin/python
'''
author:     Taylor Kessinger & Richard Neher
date:       10/07/2014
content:    generate beta coalescent trees and calculate their SFS
'''
import numpy as np
import random as rand
import scipy.special as sf
from Bio import Phylo

class betatree(object):    
    '''
    class that simulates a beta coalescent tree
    parameters:
    sample_size -- number of leaves
    alpha       -- parameter of the merger distribution 2 for Kingman, 1 for BSC
    '''
    def __init__(self, sample_size, alpha=2):
        self.alpha=alpha
        self.n = sample_size
        self.k = np.arange(self.n+1)
        self.inv_k = 1.0/np.arange(1,self.n+2)
        self.inv_kkp1 = self.inv_k/np.arange(2,self.n+3)
        self.cum_sum_inv_kkp1 = np.array([0]+np.cumsum(self.inv_kkp1).tolist())

    def init_tree(self):
        '''
        constructs the blocks that are to be merged, each leave corresponds
        to a BioPython clade object
        '''
        self.blocks = [ Phylo.BaseTree.Clade(name=str(i), branch_length=0) 
                      for i in range(self.n)]

    def coalescence_event(self):
        '''
        choose the time of the next merger,the number of blocks to merge,
        and perform the merger
        '''
        merger_size = self.whichp(len(self.blocks))
        waiting_time = self.waiting_time()
        # branch length is first accumulated. once the tree is done, take
        # differentials of parent and child branchlength to get the actual
        # branchlength
        for clade in self.blocks:
            clade.branch_length+=waiting_time

        #randomly pick some blocks
        merging_blocks = rand.sample(self.k[:len(self.blocks)], merger_size) 
        self.merge_clades(merging_blocks)

    def merge_clades(self, merging_blocks):
        '''
        creates a new clade whose children are the merging blocks
        '''
        new_clade = Phylo.BaseTree.Clade(clades = [self.blocks[i] 
                                                   for i in merging_blocks])
        new_clade.branch_length = self.blocks[merging_blocks[0]].branch_length
        # remove the merging blocks from the active blocks
        for i in sorted(merging_blocks, reverse=True):
            self.blocks.pop(i)
        self.blocks.append(new_clade)

    def clean_up_subtree(self, clade):
        '''
        calculate the branch length and number of children for each node
        '''
        if clade.is_terminal():
            clade.weight = 1
            return
        else:
            clade.weight=0
            clade.branch_length-=clade.clades[0].branch_length
            for child in clade.clades:
                self.clean_up_subtree(child)
                clade.weight+=child.weight
            return

    def waiting_time(self):
        '''
        returns the waiting time to the next merger. 
        '''
        b=len(self.blocks)
        if self.alpha==1:   #the BSC merger rate
            dt = rand.expovariate((b-1)) 
        elif self.alpha==2: #the Kingman merger rate
            dt = rand.expovariate(b*(b-1)/2.0) 
        else: # the general beta coalescent case
            exptime1 = sf.gamma(2)/(sf.gamma(2-self.alpha)*sf.gamma(self.alpha))
            exptime2 = sum([1.0*b*np.exp(-sf.gammaln(b-k+1)-sf.gammaln(k+1)+sf.gammaln(k-self.alpha)+sf.gammaln(b-k+self.alpha)) for k in range(b+1)[2:]])
            dt = rand.expovariate(exptime1*exptime2) #this product gives the Beta coalescent merger rate
        return dt

    def whichp(self,b):
        '''
        generates the merger size distribution, then samples from it.
        parameters:
            b: the number of extant lineages.
        '''
        if self.alpha == 1: #BSC case
            lambs = np.zeros(b)
            lambs[1:] = 1.0/np.arange(1,b)/np.arange(2,b+1)
            lambs = np.cumsum(lambs)
            rand = np.random.uniform(0,self.cum_sum_inv_kkp1[b-1])
            return np.where(self.cum_sum_inv_kkp1[:b] > rand)[0][0]+1
        elif self.alpha==2: #Kingman case
            return 2
        else: #other Beta coalescents
            lambs = np.zeros(b)
            for k in range(1,b): #the ugly below expression is the full probability for a given Beta coalescent with parameter self.alpha
                lambs[k] = lambs[k-1]+1.0*np.exp(-sf.gammaln(b-(k+1)+1)-sf.gammaln((k+1)+1)+sf.gammaln((k+1)-self.alpha)+sf.gammaln(b-(k+1)+self.alpha))
            rand = np.random.uniform(0,lambs[-1])
            return np.where(lambs > rand)[0][0]+1

    def coalesce(self):
        '''
        simulates the Beta coalescent process for arbitrary alpha.
        parameters:
            K0: the initial population size.
            alpha: parameter for the Beta coalescent. set to 2 for Kingman and 1 for Bolthausen-Sznitman.
        '''

        self.init_tree()
        #while the whole tree is not merged yet
        while len(self.blocks) != 1:
            self.coalescence_event()
        
        self.clean_up_subtree(self.blocks[0])
        self.BioTree = Phylo.BaseTree.Tree(root=self.blocks[0])


class SFS(betatree):
    '''
    class the generates many trees and accumulates an SFS
    '''
    def __init__(self, sample_size, alpha=2):
        betatree.__init__(self,sample_size, alpha)
        self.alleles=[]

    def glob_trees(self, ntrees=10):
        '''
        generate many trees, accumulate the SFS
        parameters:
        ntrees -- number of trees to generate
        '''
        for ti in xrange(ntrees):
            self.coalesce()
            self.alleles.append([(clade.weight,clade.branch_length) 
                                 for clade in self.BioTree.get_terminals()
                                 +self.BioTree.get_nonterminals()])

    def getSFS(self, ntrees = 10):
        '''
        calculate an SFS based on ntrees trees
        '''
        self.glob_trees(ntrees)
        self.sfs = np.zeros(self.n+1)
        for aset in self.alleles:
            for w,l in aset:
                self.sfs[w]+=l

    def binSFS(self, bins=10):
        '''
        use the precalcutated SFS and bin it.
        '''
        pass

if __name__=='__main__':
    import matplotlib.pyplot as plt
    myT = betatree(100,2)
    myT.coalesce()
    Phylo.draw(myT.BioTree)

    myT = betatree(100,1)
    myT.coalesce()
    Phylo.draw(myT.BioTree)

    myT = betatree(100,1.5)
    myT.coalesce()
    Phylo.draw(myT.BioTree)

    mySFS = SFS(100,alpha=1)
    mySFS.getSFS(ntrees=1000)
    plt.figure()
    plt.plot(np.linspace(0,1,mySFS.n+3)[1:-1], mySFS.sfs)
    plt.yscale('log')


