from multiprocessing import Pool
from tqdm import tqdm
from timeit import timeit
import time
import itertools
import numpy as np
import random
import toolbox as tb
import nodf_max_hill_climb
import greedySolver2
import hill_climb

def secondEntry(lst):
    return lst[1]

def genRandGraph(NodesA, NodesB, Edges):
    graph = np.zeros((NodesA, NodesB))
    #Initial fill:
    graph[:, 0] = 1.0
    graph[0, :] = 1.0
    filled = NodesA + NodesB - 1
    # Fill the rest:
    rows = list(range(1, NodesA))
    cols = list(range(1, NodesB))
    idxes = list(itertools.product(rows, cols))
    probs = []
    for idx1, idx2 in idxes:
        #probs.append(1.0 / (idx1 + idx2))
        probs.append(1.0)
    probs = np.array(probs) / np.sum(probs)
    myIdxes = np.random.choice(range(len(idxes)), Edges - filled, p = probs)
    for i in myIdxes:
        idxX, idxY = idxes[i]
        graph[idxX, idxY] = 1.0
    return graph

def crossover(parentA, parentB):
    Edges = np.sum(parentA)
    NodesA, NodesB = parentA.shape
    #decision = (np.random.random((NodesA, NodesB)) < 0.5).astype(float)
    decision = np.random.random((NodesA, NodesB))
    child = parentA * decision + parentB * (1.0 - decision)
    while(np.sum(child) > Edges):
        i = np.random.randint(1, NodesA)
        j = np.random.randint(1, NodesB)
        child[i,j] = 0
    while(np.sum(child) < Edges):
        i = np.random.randint(1, NodesA)
        j = np.random.randint(1, NodesB)
        child[i,j] = 1
    return child

def mutate(graph):
    graph = tb.findNeighborMtx(graph)
    return graph

def local_optimisation(graph):
    graph = hill_climb.full_hill_climb(graph)
    return graph

def simulateGeneration(population, fitness, elite=0.05, rmut=0.08, hcmut=0.08):
    # I assume that the population is sorted by fitness
    population = list(population)
    fitness = np.array(fitness)
    fitness = fitness / np.sum(fitness)
    Nelite = round(len(population) * elite)
    Nmutants = round(len(population) * rmut)

    # Implementing elites
    children = population[0:Nelite]

    # Compute the random mutants
    with Pool(4) as pool:
        mutants = pool.map(mutate, population[0:Nmutants])
    children.extend(mutants)

    # Compute the children.
    Nchild = len(population) -Nelite -Nmutants
    parentsA = population[0:Nchild]
    parentsB = parentsA
    np.random.shuffle(parentsA)
    np.random.shuffle(parentsB)
    #parentsA = np.random.choice(range(len(population)), Left, p = fitness)
    #parentsB = np.random.choice(range(len(population)), Left, p = fitness)
    params = list(zip(parentsA, parentsB))
    #for (Aidx, Bidx) in list(zip(parentsA, parentsB)):
    #    params.append((population[Aidx], population[Bidx]))
    # Computing the children:
    with Pool(4) as pool:
        children2 = pool.starmap(crossover, params)
    children.extend(children2)

    # Hill climb on the population:
    Noptim = int(round(len(children) * hcmut))
    if(Noptim > 0):
        idxes = np.random.choice(len(children), Noptim, replace = False)
        params = []
        for idx in idxes:
            params.append(children[idx])
        with Pool(4) as pool:
            opt_children = pool.map(local_optimisation, params)
            opt_children = pool.map(tb.optimize_mtx_order, opt_children)
        for idx, child in list(zip(idxes, opt_children)):
            children[idx] = child


    return children

def optimise(NodesA, NodesB, Edges, verbose = True):
    popSize = 100
    params = []
    for i in range(popSize):
        params.append([NodesA, NodesB, Edges])

    if(verbose):
        print("Generating graphs")
    with Pool(4) as pool:
        population = list(tqdm(pool.imap(greedySolver2.greedy_gen, params), total = popSize))

    if(verbose):
        print("Starting genetic algorithm")
    Gens = 500
    optgraph = genRandGraph(NodesA, NodesB, Edges)
    optnodf = tb.nodf(optgraph)
    for i in tqdm(range(Gens)):
        with Pool(4) as pool:
            fitness = pool.map(tb.nodf, population)
        popfitness = list(zip(population, fitness))
        popfitness.sort(key = secondEntry, reverse = True)
        population, fitness = zip(*popfitness)
        if(fitness[0] > optnodf):
            optnodf = fitness[0]
            optgraph = population[0]
            if(verbose):
                print("Optimal nodf value is {}".format(optnodf))
        rmut = 0.2
        hcmut = 0.1
        population = simulateGeneration(population, fitness, rmut = rmut, hcmut = hcmut)
    return optgraph

def benchCreation(NodesA, NodesB, Edges, N):
    # Benchmark the creation of graphs:
    t1 = time.time()
    params = []
    for i in range(N):
        params.append((NodesA, NodesB, Edges))
    with Pool() as pool:
        graphs = pool.starmap(genRandGraph, params)
    t2 = time.time()
    tdiff = t2 - t1
    print("Generation: \t{}s".format(tdiff / N))

def benchCrossOver(NodesA, NodesB, Edges, N):
    my_graphs = []
    for i in range(32):
        my_graphs.append(genRandGraph(NodesA, NodesB, Edges))
    t1 = time.time()
    params = []
    for i in range(N):
        idx1 = np.random.randint(32)
        idx2 = np.random.randint(32)
        g1 = my_graphs[idx1]
        g2 = my_graphs[idx2]
        params.append((g1, g2))
    with Pool() as pool:
        children = pool.starmap(crossover, params)
    t2 = time.time()
    tdiff = t2 - t1
    print("Crossover: \t{}s".format(tdiff / N))

def benchMutation(NodesA, NodesB, Edges, N):
    my_graphs = []
    for i in range(32):
        my_graphs.append(genRandGraph(NodesA, NodesB, Edges))
    t1 = time.time()
    params = []
    for i in range(N):
        idx1 = np.random.randint(32)
        g1 = my_graphs[idx1]
        params.append(g1)
    with Pool() as pool:
        mutants = pool.map(mutate, params)
    t2 = time.time()
    tdiff = t2 - t1
    print("Mutation: \t{}s".format(tdiff / N))

def benchFitnessEval(NodesA, NodesB, Edges, N):
    my_graphs = []
    for i in range(32):
        my_graphs.append(genRandGraph(NodesA, NodesB, Edges))
    t1 = time.time()
    params = []
    for i in range(N):
        idx1 = np.random.randint(32)
        g1 = my_graphs[idx1]
        params.append(g1)
    with Pool(1) as pool:
        nodfs = pool.map(tb.nodf, params)
    t2 = time.time()
    tdiff = t2 - t1
    print("Fitness: \t{}s".format(tdiff / N))

if(__name__ == "__main__"):
    NodesA = 16
    NodesB = 44
    Edges = 278

    #benchCreation(NodesA, NodesB, Edges, N)
    #benchCrossOver(NodesA, NodesB, Edges, N)
    #benchMutation(NodesA, NodesB, Edges, N)
    #benchFitnessEval(NodesA, NodesB, Edges, N)
    optimise(NodesA, NodesB, Edges)
