from multiprocessing import Pool
from tqdm import tqdm
import algGenetic
import greedySolver2
import hill_climb
import multiprocessing
import numpy as np
import toolbox as tb
import copy

def neighbor_nodf(mtx, MT, F, deg_min, neg_del, sums, o_pos, z_pos):
    """
    Computes the neighbor of a matrix (i.e. a matrix with a hemming
    distance of 1) along with the associated nodf value.
    Output is the new nodf value.
    Note that the input parameters will be modified in this method.
    """
    # Remove a link:
    nodf0, sums = tb.nodf_one_link_removed(mtx, MT, F, deg_min, neg_del, sums, o_pos)
    # Add a different link:
    nodf1, sums = tb.nodf_one_link_added(mtx, MT, F, deg_min, neg_del, sums, z_pos)
    # one link has been moved from o_idx to z_idx and the new nodf value
    # can be returned
    return nodf1, sums

def cost(nodf):
    """
    A very simple cost function
    """
    return 1.0 - nodf

def sim_anneal_step(mtx, temp, iters, support_data):
    """
    One step of the simmulated annealing algorithm.
    Supply a initial solution and nodf meta-data.
    Output will be the optimal matrix found in this process
    along with it's cost
    """
    # Do the reqired dice rolls for this algorithm vectorised and in advance:
    decision = np.random.rand(iters)
    MT, F, DM, ND, S = support_data

    opt_mtx = np.copy(mtx)
    opt_cost = cost(tb.nodf(mtx))
    old_cost = opt_cost
    ones = tb.get_valid_ones(mtx)
    zeros = tb.get_promising_zeros(mtx)
    opt_time = -1
    for i in range(iters):
        o_idx = np.random.randint(len(ones))
        z_idx = np.random.randint(len(zeros))

        o_pos = ones[o_idx]
        z_pos = zeros[z_idx]
        # compute the new
        new_nodf, S = neighbor_nodf(mtx, MT, F, DM, ND, S, o_pos, z_pos)
        new_cost = cost(new_nodf)
        if(new_cost < opt_cost):
            opt_cost = new_cost
            opt_mtx = np.copy(mtx)
            opt_time = i
        acc_prob = tb.acceptProb(old_cost, new_cost, temp)
        if(decision[i] <= acc_prob):
            # accept the new solution by not changing it back
            old_cost = new_cost
            # modify the zero and ones list accordingly:
            # ones = tb.get_valid_ones(mtx)
            ones = tb.get_valid_ones(mtx)
            zeros[z_idx] = o_pos
        else:
            # reject the new solution by reverting back
            a, S = neighbor_nodf(mtx, MT, F, DM, ND, S, z_pos, o_pos)
            # old_cost will not be updated!
    support_data = [MT, F, DM, ND, S]
    return [opt_mtx, opt_cost, support_data]

def sim_anneal_opt(mtx, alpha = 0.9985, iters = 6, init_temp = 0.25, min_temp = 10**(-4)):
    """
    A Simulated annealing algorithm. The default values are copied from my
    supervisors implementation. Given these defaults the user only needs to
    supply a initial solution [mtx_0] to get a well performing simulated
    annealing algorithm.
    Other parameters:
        alpha: cooling factor
        iters: the number of iterations at each temp level
        bar: show a progress bar
    Output will be the optimal solution encountered in the whole algorithm.

    Rest is clear.
    """
    NodesA, NodesB = mtx.shape
    tot_steps = int(np.ceil(np.log(min_temp / init_temp) / np.log(alpha))) * iters
    cool_steps = int(np.ceil(np.log(min_temp / init_temp) / np.log(alpha)))

    # Initialise optimal parameters:
    opt_mtx = np.copy(mtx)
    opt_nodf = tb.nodf(opt_mtx)
    opt_cost = cost(opt_nodf)
    temp = init_temp
    eps = 10**(-12)

    # Initialise support data:
    support_data = tb.init_nodf(mtx)
    opt_support_data = copy.deepcopy(support_data)
    for i in tqdm(range(cool_steps)):
        [new_mtx, new_cost, support_data] = sim_anneal_step(mtx,temp,iters,support_data = support_data)
        # print(1.0 - opt_cost, temp, tb.nodf(mtx))
        if(new_cost + eps < opt_cost):
            # Found a new optimal matrix.
            opt_mtx = np.copy(new_mtx)
            # Try hillclimbing to improve it!
            opt_mtx = hill_climb.full_hill_climb(opt_mtx,R = 3, multithread = True)
            # Update optimal parameters
            opt_nodf = tb.nodf(opt_mtx)
            opt_cost = cost(opt_nodf)
            opt_support_data = tb.init_nodf(opt_mtx)
            print("New best NODF value found: {}".format(opt_nodf))
        # Test to see if we should go back to the opt. solution
        acc_prob = tb.acceptProb(opt_cost, new_cost, temp)
        if(np.random.random() > acc_prob):
            # Go back to the optimal solution:
            mtx = np.copy(opt_mtx)
            support_data = copy.deepcopy(opt_support_data)
        temp = temp * alpha
    return opt_mtx

if(__name__ == "__main__"):
    NodesA = 14
    NodesB = 13
    Edges = 52
    mtx = greedySolver2.greedySolve(NodesA, NodesB, Edges)
    mtx = hill_climb.full_hill_climb(mtx)
    mtx = sim_anneal_opt(mtx)
    nodf = tb.nodf(mtx)
    print(mtx.sum())
    print("NODF = {}".format(nodf))
