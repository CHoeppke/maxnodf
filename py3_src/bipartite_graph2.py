import numpy as np
import matplotlib.pyplot as plt
from nestedness_calculator import NestednessCalculator

class BipartiteGraph(object):

    Nnodes = 0
    NnodesA = 0
    NnodesB = 0
    Edges = []
    Mtx = False
    MtxUpdated = True

    def __init__(self, NinA, NinB):
        self.NnodesA = NinA
        self.NnodesB = NinB
        self.Nnodes = NinA + NinB
        self.Edges = []

    def getEdges(self):
        return self.Edges

    def getNodesA(self):
        return self.NnodesA

    def getNodesB(self):
        return self.NnodesB

    def hasEdge(self, NodeA, NodeB):
        hasEdge = False
        if([NodeA, NodeB] in self.Edges):
            hasEdge = True
        return hasEdge

    def getMatrix(self):
        if(self.MtxUpdated == True):
            Mtx = np.array([[0 for k in range(self.NnodesB)] for j in
                range(self.NnodesA)])
            # Now fill the graph with the appropriate entiries
            for [Aidx, Bidx] in self.Edges:
                # Debug message
                # print("Has edge at position ({}, {})".format(Aidx, Bidx))
                # The matrix has the be symmetrical.
                Mtx[Aidx, Bidx] = 1.0
            self.Mtx = Mtx
            self.MtxUpdated = False
        # Speedup if the matrix has already been computed
        return self.Mtx

    def addEdge(self, node1, node2):
        if(self.hasEdge(node1, node2) == False):
            self.Edges.append([node1, node2])
            self.MtxUpdated = True

    def copy(self):
        newGraph = BipartiteGraph(self.NnodesA, self.NnodesB)
        for [n1,n2] in self.Edges:
            newGraph.addEdge(n1,n2)
        return newGraph

    def isValidNetwork(self):
        Mtx = self.getMatrix()
        RowSum = np.matmul(Mtx, np.ones((self.NnodesB, 1)))
        for rsum in RowSum:
            if(rsum == 0.0):
                return False
        ColSum = np.matmul(np.transpose(Mtx), np.ones((self.NnodesA, 1)))
        for csum in ColSum:
            if(csum == 0.0):
                return False
        # Both the RowSum and ColSum test is okay.
        # the means that the graph represents a valid network
        return True

    # Compute the nestedness
    def NODF(self):
        # First compute all row pairs and their DF, PO values
        # Then use the transpose of the matrix to do the same of the
        # column pairs
        print("Implement NODF function")

    def get_row_marginal_totals(self):
        Mtx = self.getMatrix()
        MT = np.zeros((self.NnodesA, 1))
        ONES = np.ones((self.NnodesB,1))
        MT = np.matmul(Mtx, ONES)
        return MT

    def get_col_marginal_totals(self):
        Mtx = self.getMatrix()
        ONES = np.ones((self.NnodesA,1))
        MT = np.matmul(np.transpose(Mtx), ONES)
        return MT

    def get_row_DF_and_PO(self):
        # First select all row pairs with i < j
        MT = self.get_row_marginal_totals()
        DF = np.zeros((self.NnodesA, self.NnodesB))
        for j in range(1, self.NnodesA):
            for i in range(0, self.NnodesB):
                # Now compute the Marginal totals for i and j
                if(MT[j] < MT[i]):
                    DF[i,j] = 100.0
                else:
                    DF[i,j] = 0.0

        # Now compute the N and PO values:
        print("Implement get_row_DF_and_PO")

# Computes a Graph from a given matrix.
def getGraphFromMatrix(mtx):
    NodesA, NodesB = mtx.shape
    newGraph = BipartiteGraph(NodesA, NodesB)
    for i in range(NodesA):
        for j in range(NodesB):
            if( mtx[i,j] >= 0.5):
                newGraph.addEdge(i, j)
    return newGraph

def graphInList(graphList, graph):
    for graphB in graphList:
        if(np.array_equal(graphB.getMatrix(), graph.getMatrix())):
            return True

    return False

def generateAllGraphs(NvertA, NvertB, nEdge):
    print("Generating graphs with {}, {}, {}".format(NvertA, NvertB, nEdge))
    if(nEdge == 0):
        myGraph = BipartiteGraph(NvertA, NvertB)
        return [myGraph]
    elif(nEdge > 0):
        #Add get list of smaller graphs and add an edge:
        smallGraphs = generateAllGraphs(NvertA, NvertB, nEdge - 1)
        myGraphs = []
        for myGraph in smallGraphs:
            for Aidx in range(NvertA):
                for Bidx in range(NvertB):
                    # Create a copy of the current graph
                    newGraph = myGraph.copy()
                    # print("Created a copy graph")
                    if(newGraph.hasEdge(Aidx, Bidx) == False):
                        newGraph.addEdge(Aidx, Bidx)
                        # The graph I wanted to create is valid. So I will
                        # add the new graph to the list myGraphs
                        # print("Added an edge")
                        if(graphInList(myGraphs, newGraph) == False):
                            myGraphs.append(newGraph)
                            # print("Added the graph")
        return myGraphs
    return NULL

if(__name__ == '__main__'):
    # Try and create every BipartieteGraph with
    # Nvert nodes and Nedge edges. This should help me get a first overview
    # of the bigger problem

    # Only give 4x3 graphs with 4 edges
    NnodesA = 4
    NnodesB = 3
    graphList = generateAllGraphs(NnodesA, NnodesB, 5)
    x = str(len(graphList))
    print("There are {} graphs with the given spec.".format(x))
    # For 4x3 graphs this should give 220 graphs.
    # So far this works!

    # Filter the graphList for valid networks:
    validGraphs = []
    for graph in graphList:
        if(graph.isValidNetwork()):
            validGraphs.append(graph)

    x = str(len(validGraphs))
    print("There are {} valid networks with the given spec.".format(x))

    # Now let me compute the nestedness metric for all of those graphs.
    # I have found a NODF computing script on GitHub that I will use
    # for this task.
    maxNest = -100.0
    optGraph = np.zeros((NnodesA, NnodesB))
    for myGraph in validGraphs:
        graphMtx = myGraph.getMatrix()
        ncalc = NestednessCalculator(graphMtx)
        nodf = ncalc.nodf(graphMtx)
        if(nodf > maxNest):
            maxNest = nodf
            optGraph = graphMtx

    print("The maximal nestedness is {}.".format(maxNest))
    print("The maximal nestedness is attained using the graph:")
    print(optGraph)

