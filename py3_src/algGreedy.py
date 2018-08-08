from tqdm import tqdm
import numpy as np
import time
import toolbox as tb
import greedySolver2
import hill_climb
import sys

def optimise(NodesA, NodesB, Edges, verbose = True):
    mtx = greedySolver2.greedySolve(NodesA, NodesB, Edges)
    return mtx

