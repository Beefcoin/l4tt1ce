import numpy as np

""" This file contains all necessary Lattice Descriptors with all needed information such
as Directions, Discretized Velocity Directions, Weights and opposite directions for the
velocity set, to make calculations easier """

class D2Q9:
    #this is just a basic D2Q9 Lattice Descriptor
    Q = 9 #directions
    e = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
                  [1, 1], [-1, 1], [-1, -1], [1, -1]]) #discretized velocity set
    w = np.array([4/9] + [1/9]*4 + [1/36]*4) #weights (is this correct? I hope so)
    opp = [0, 3, 4, 1, 2, 7, 8, 5, 6] #opposite directions for easier calculations
