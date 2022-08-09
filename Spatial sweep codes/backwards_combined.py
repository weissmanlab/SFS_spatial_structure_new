#!/usr/bin/env python3
'''
Coalescent simulation with recombination.
Takes the output of the forward time simulations as input. 
Read the frequency of beneficial mutation and go backwards in time until n individuals all coalesce.
As time goes back, record the number of leaves for each node. 
The neutral mutations arise as a Poisson process with rate Un = 1 under our assumptions. 
The expected SFS is the histogram of the leaf counts 

Pre-sweep coalescent process can account for more than one merger in a single generation
when there are many individuals are left (large compared to sqrt(2 * N))

Code can handle 2-D populations by converting the 2-D array into a 1-D array with 
migration patterns calculated using a single index. 


'''

import numpy as np
from multiprocessing import Pool
from multiprocessing import Process
import sys
from numpy import random
from functions_combined import *
import time

start = time.time()



L = int(sys.argv[1]) # number of demes
rho = int(sys.argv[2]) # deme capacity
s = float(sys.argv[3]) # selection coefficient
m = float(sys.argv[4]) # migration rate
l0 = int(sys.argv[5]) # location of the origin of the sweep
nbase = int(sys.argv[6]) # sample size for the coalescent simulation
N_SFS = int(sys.argv[7]) # number of coalescent simulation we run for ensemble average.
T_after_fix = int(sys.argv[8]) # number of generations between the end of the sweep and sampling
Nforw = int(sys.argv[9]) # number of forward simulations, used for parallel jobs
dimension = int(sys.argv[10]) # To set spatial dimension of population (1D or 2D) 
r = float(sys.argv[11]) # recombination rate

extra_gen_cutoff = sys.argv[12] if len(sys.argv) >= 13 else int(L ** 2 / m / rho)

## python backwards_combined.py L N s m l0 n_base N_SFS T_after_fix Nforw dimensions
## python backwards_combined.py 500 20000 0.05 0.250 1 10000 4 3000 1 1


if (dimension == 1):
    print ('1_D')
    fname = 'L={}_N={}_s={:.6f}_m={:.6f}_{}.txt'.format(L, rho, s, m, Nforw)
    print(fname)
    lines = np.loadtxt(fname, dtype = np.int64)

elif(dimension == 2):
    print('2-D')
    ##Read input file and flatten lines[i] into 1-D array here for use within existing framework
    fname = 'L={}_N={}_s={:.3f}_m={:.2f}_l0={}_Nforw={}.txt'.format(L, rho, s, m, l0, Nforw)
    #print(fname)
    lines = np.loadtxt(fname, dtype = np.int64)
    tfinal = int(len(lines) / L)
    lines = lines.reshape(tfinal, L, L)  ###lines = Data from forward simulation
    lines = [lines[i].flatten() for i in range(len(lines))]
    print(tfinal)    
else:
    print("Dimension not recognised")
    quit()

def runner(idx):
    '''
    idx: argumemt to spawn multiple jobs. Not relevant for the body of the simulation. 
    
     Reads a generation, starting from the last one. Tracks the sampled indviduals backwards in time 
    till all mergers have happened and then calculates the branch length to get the SFS.
       '''
    rho_e = lines[-1]    ##Mutants in a deme. We pick the absolute last generation we simulated forward in time till to start simulating backwards. 
    n = nbase
    SFS = np.zeros(n)
    rho_e = (rho_e).astype(np.int64)

    '''Creating the primary data structure by sampling'''
     # data structure('individuals') format will be [mut_type, deme index, individual index(inside the deme)] in a row and as many rows as the individuals we sample
    
    individuals = sample_data(rho_e, n, rho)
    unique, leaf_counts = np.unique(individuals, axis = 0, return_counts = True) #To avoid overcounting of repeated locations and ensure we actually sample what we want
 
    # The SFS is the histogram of these leaf counts counted over all generations
    hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))


    '''Starting simulations backward in time now'''    
 
    t_after_fix = 0    
    '''Going backwards in time for user given generations after sweep has fixed'''
    #Since the population doesnt change once sweep has fixed, we don't need to change the input data for this
    
    while (t_after_fix < T_after_fix) and (len(individuals) > 1):
        SFS += hist
        t_after_fix += 1
        # Since every individual has the mutant allele after fixation, recombination step is not needed
        individuals_post_migration = migration(rho_e, individuals, rho, m, L, dimension)
        individuals, leaf_counts = coalescent(individuals_post_migration, leaf_counts)
        hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))

    '''Going back in time during the sweep, all the way back to the time when mutation is first seeded'''
    line_num = -1
    
    while (len(individuals) > 1) and (line_num > -len(lines)):  ##The maximum we can go backward in time is till T2
        #print('sweep')
        SFS += hist
        line_num -= 1

        rho_e_parent = (lines[line_num]).astype(np.int64) ##Getting the parent generation from each time step of our forward simulation
        individuals_post_migration = migration(rho_e_parent, individuals, rho, m, L, dimension)
        individuals_post_recombination = recombination(rho_e_parent, individuals_post_migration, rho, r)
        individuals, leaf_counts = coalescent(individuals_post_recombination, leaf_counts)
        hist, bin_edges = np.histogram(leaf_counts,bins = np.arange(1, n + 2))
        rho_e = rho_e_parent
    

    '''
    The first round of coalescence simulation ends when we get to the time
    when the beneficial mmutation arose, and all the left over individuals
    are wild type (WT). 
    If there is no recombination, the simulation should end here since the MRCA 
    will be the benficial mutation that gave rise to the sweep. 
    
    From this point till we reach the MRCA, the coalescence will be extremely slow. The MRCA can be different
    from the beneficial mutation because of recombination. 
    We run the same kind of simulation until the individuals
    disperse for L^2 / m / N. We speed up the simulation by recording the number
    of generations between the coalescence events.
    '''
    

    left_individuals = len(individuals)  ##Individuals left to still coalesce in time before the mutation arose
    # Make sure that all individuals are WT (mutation type = 0)
    mut_types, deme_arr, ind_in_deme_arr = individuals.T
    mut_types = np.zeros_like(mut_types)
    mut_types = (mut_types).astype(np.int64)
    individuals = np.vstack((mut_types, deme_arr, ind_in_deme_arr)).T

    branch_len = 0 # number of generations until first merging event
    extra_gen = 0 # extra run time before stopping the coalescent simulation.


    while left_individuals > 1 and extra_gen < extra_gen_cutoff:  ###Will it be same for 2-D
        #print('extra gens')
        # Only update call np.histogram when there is nonzero number of coalescent events.
        SFS += hist
        extra_gen += 1
        rho_e_parent_pre_sweep = np.zeros_like(rho_e)
        rho_e_parent_pre_sweep = (rho_e_parent_pre_sweep).astype(np.int64)

        # Before the beneficial mutation appeared, everyone is WT. Thus, rho_e_parent is just an array of zeros.
        # Also, we can skip the recombination step, since the genotype will always be zero for everyone.
        individuals_post_migration = migration(rho_e_parent_pre_sweep, individuals, rho, m, L, dimension)
        individuals, leaf_counts = coalescent(individuals_post_migration, leaf_counts)
        current_individuals_counts = len(individuals)
        if current_individuals_counts < left_individuals:
            hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))
        left_individuals = current_individuals_counts
        
        
    '''    
    After the individuals disperse, we may ignore the spatial structure.
    If there are more than sqrt(2 * N * L) individuals left, it is likely
    to have more than 1 merging event in a single generation. Thus, we manually
    draw a parent for every individual and find if some of them are the same.
    Once there are much smaller number of individuals left, it will take a long
    time to coalesce. Thus, we draw T2 (time until the first coalescent) from
    geometric probablity distribution. (This is found from probablity of choosing a different
    parent for every individual left.)
    '''

    
    while left_individuals > 1:
        T2 = random.exponential(rho * L / ((left_individuals - 1) * left_individuals / 2)) ###Drawing the T2 from geometric distribution
        SFS += hist * T2 ##That many branches will exist


        # Now choose two random branches that will coalesce in T2 generations. 
        coal_inds = random.choice(range(left_individuals), size = 2, replace = False)
        # Like in the coalescent function, as two branches merge, two leaf counts are summed up and the number of remaining individuals decreases by 1.
        leaf_counts_copy = leaf_counts
        leaf_counts_copy[coal_inds[0]] += leaf_counts[coal_inds[1]]
        leaf_counts_copy[coal_inds[1]] = 0
        leaf_counts = leaf_counts_copy[leaf_counts_copy > 0]
        hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))

        left_individuals -= 1


    f = np.arange(1, n + 1) / n  ##Calculating frquency of mutant
    H = np.sum(2 * f * (1 - f) * SFS) / np.sum(SFS)  ##Calculating Heterezygosity
    return SFS, H


#############################################################
'''Main body to create multiple backward trees and average'''   
#############################################################


if __name__ == '__main__':
    p = Pool(20)  # 20 is the nunmber of cores we used to spawn jobs parallely on, adjust as needed
    ret = p.map(runner, range(N_SFS))  ###Number of simulations you want to average over    
    SFS_items = [r[0] for r in ret]
    H_items = [r[1] for r in ret]
    # Averaging all SFSs 
    SFS = np.sum(SFS_items, axis=0)
    SFS /= N_SFS
    np.savetxt('expected_SFS_L={}_N={}_s={:.3f}_m={:.2f}_r={:.2e}_nsample={}_t_after_fix={}_Nback={}_Nforw={}.txt'.format(L,
                rho, s, m, r, nbase, T_after_fix, N_SFS, Nforw), SFS)


    end = time.time()
    print(start-end) 
    
