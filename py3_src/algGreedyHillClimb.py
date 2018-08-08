from tqdm import tqdm
import numpy as np
import time
import toolbox as tb
import greedySolver2
import hill_climb

def optimise(NodesA, NodesB, Edges, verbose = True):
    mtx = greedySolver2.greedySolve(NodesA, NodesB, Edges)
    mtx = hill_climb.full_hill_climb(mtx, R = 3)
    return mtx
