import algGreedy
import algGreedyHillClimb
import algSimulatedAnneal
import csv
import numpy as np
import sys
import toolbox as tb


if __name__ == "__main__":
    if(len(sys.argv) != 5):
        print("Invalid function call!")
    else:
        NodesA = int(sys.argv[1])
        NodesB = int(sys.argv[2])
        Edges = int(sys.argv[3])
        quality = int(sys.argv[4])
        mtx = np.zeros((NodesA, NodesB))
        if(quality == 0):
            # Use simple greedy algorithm:
            mtx = algGreedy.optimise(NodesA, NodesB, Edges)
        elif(quality ==1):
            mtx = algGreedyHillClimb.optimise(NodesA, NodesB, Edges)
        elif(quality == 2):
            mtx = algSimulatedAnneal.optimise(NodesA, NodesB, Edges)
        else:
            print("Invalid quality parameter")
        # Output matrix to csv file
        mtx = mtx.astype(int)
        np.savetxt("outputmtx.csv",mtx, delimiter = ",", fmt = "%i")
