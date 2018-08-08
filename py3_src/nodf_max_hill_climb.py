from bipartite_graph2 import BipartiteGraph
from multiprocessing import Pool
import bipartite_graph2
import bisect
import collections
import matplotlib.pyplot as plt
import numpy as np
import random
import time
import toolbox
from tqdm import tqdm

def isValidMatrix(Mtx):
    NodesA, NodesB = Mtx.shape
    M1 = np.matmul(Mtx, np.ones((NodesB, 1)))
    M2 = np.matmul(np.transpose(Mtx), np.ones((NodesA, 1)))

    Valid1 = (M1 > 0.0)
    Valid2 = (M2 > 0.0)

    Valid1_1 = np.all(Valid1)
    Valid2_2 = np.all(Valid2)
    return (Valid1_1 and Valid2_2)

# Where do we need to add a node to maximise the
# NODF measure?
def max_nodf_add_one_interaction(mtx):
    my_zeros = np.where(mtx == 0.0)
    bestNODF = -100.0
    bestMTX = mtx
    for i in range(len(my_zeros[0])):
        xpos = my_zeros[0][i]
        ypos = my_zeros[1][i]
        newMTX = np.copy(mtx)
        newMTX[xpos, ypos] = 1.0
        if(isValidMatrix(newMTX)):
            newNODF = toolbox.nodf(newMTX)
            if(newNODF > bestNODF):
                bestNODF = newNODF
                bestMTX = newMTX
    return [bestMTX, bestNODF]

# swap a zero and a one entry in the matrix and compute the
# new nodf measure
def swapEntriesNODF(params):
    mtx = params[0]
    zpos = params[1]
    opos = params[2]
    mymtx = np.copy(mtx)
    mymtx[zpos[0], zpos[1]] = 1.0
    mymtx[opos[0], opos[1]] = 0.0
    if(isValidMatrix(mymtx)):
        newnodf = toolbox.nodf(mymtx)
        return newnodf
    else:
        return -100.0

# Get all the zero and one entry positons in the matrix
# and put them into a usable format
def getZerosAndOnes(mtx):
    zList = np.array(np.where(mtx == 0.0)).T
    oList = np.array(np.where(mtx == 1.0)).T
    return [zList, oList]

def hill_climb(mtx):
    mtx_backup = np.copy(mtx)
    oldNODF = toolbox.nodf(mtx)
    NodesA, NodesB = mtx.shape
    mtx_has_changed = False
    [zList, oList] = getZerosAndOnes(mtx)
    params = []
    for one_elem in oList:
        for xshift, yshift in list(zip([-1, 1], [-1, 1])):
            shift = np.array([xshift, yshift])
            new_pos = one_elem + shift
            if(np.all(new_pos >= 0) and new_pos[0] < NodesA and new_pos[1] < NodesB):
                if(mtx[new_pos[0], new_pos[1]] == 0.0):
                    params.append([mtx, new_pos, one_elem])
    with Pool(4) as pool:
        newNODFList = list(tqdm(pool.imap(swapEntriesNODF, params),total=len(params)))

    maxidx = newNODFList.index(max(newNODFList))
    zpos = params[maxidx][1]
    opos = params[maxidx][2]
    mtx[zpos[0], zpos[1]] = 1.0
    mtx[opos[0], opos[1]] = 0.0
    newNODF = newNODFList[maxidx]
    if(newNODF > oldNODF):
        return [mtx, True]
    else:
        return [mtx_backup, False]

