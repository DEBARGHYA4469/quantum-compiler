#   ----------------------------------------  Quantum Gates --------------------------------------------------------------------- 

from phase_balance import *
import numpy as np
import math 

PI = math.acos(-1.0)

# Library Gates ............................
H = np.array([[1 , 1],[1 , -1]])/math.sqrt(2) # Hadamard  
T = np.array([[math.cos(PI/8)-1j*math.sin(PI/8) , 0],[0 , math.cos(PI/8)+1j*math.sin(PI/8)]]) # Phase Gate 

X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]]) 
Z = np.array([[1,0],[0,-1]])

