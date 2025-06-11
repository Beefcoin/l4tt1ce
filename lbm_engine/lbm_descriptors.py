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

import numpy as np

class D3Q19:
    # This is a basic D3Q19 Lattice Descriptor (thanks ChatGPT, this would have been ass to make by hand)
    Q = 19  # number of discrete velocity directions
    
    #discret velocity set
    e = np.array([
        [ 0,  0,  0], 
        [ 1,  0,  0], 
        [-1,  0,  0],  
        [ 0,  1,  0],  
        [ 0, -1,  0],  
        [ 0,  0,  1],  
        [ 0,  0, -1], 
        [ 1,  1,  0],  
        [-1,  1,  0],  
        [-1, -1,  0],  
        [ 1, -1,  0],  
        [ 1,  0,  1],  
        [-1,  0,  1],  
        [-1,  0, -1],  
        [ 1,  0, -1],  
        [ 0,  1,  1],  
        [ 0, -1,  1],  
        [ 0, -1, -1],  
        [ 0,  1, -1],  
    ])

    #  weights (from standard D3Q19 lattice)
    w = np.array([
        1/3,               
        1/18, 1/18,         
        1/18, 1/18,
        1/18, 1/18,
        1/36, 1/36,         
        1/36, 1/36,
        1/36, 1/36,
        1/36, 1/36,
        1/36, 1/36,
        1/36, 1/36
    ])

    #opposite directions for easier calculations
    opp = [
         0, 
         2, 1, 4, 3, 6, 5,
         9, 10, 7, 8,
        13, 14,11,12,
        17, 18,15,16
    ]
