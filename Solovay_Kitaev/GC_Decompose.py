# --------------------------------------Group Commutator Decompose ---------------------------------------------------------

import math  
import numpy as np
from Gates import *
from util import *
  

def solve_phi(theta):
	a=math.sin(theta/2.0)
	return 2*math.asin(math.sqrt(math.sqrt(0.5*(1 + math.sqrt(1-a**2))))) 

def commutator(V,W):
	return np.matmul(V,np.matmul(W,np.matmul(dagger(V),dagger(W))))

def gc_decompose(U): # decompose U = SVWV*W*S* and return V',W' where V = SVS* ,  W' = SWS*
	theta = 2*math.acos(U[0][0].real) # calculate 
	phi = solve_phi(theta) 
	V = SU2(rotateX(phi))
	W = SU2(rotateY(phi))
	w = dagger(W)
	v = dagger(V)
	K = np.matmul(V,np.matmul(W,np.matmul(v,w)))
	w1,v1=diagonalize(U)
	w2,v2=diagonalize(K)
	
	S = SU2(np.matmul(v1,dagger(v2)))
	Vtilde = np.matmul(S,np.matmul(V,dagger(S)))
	Wtilde = np.matmul(S,np.matmul(W,dagger(S)))
	return Vtilde,Wtilde

def rotateX(phi): # V 
	X_U = np.array([[0,1],[1,0]])
	V = I*math.cos(phi/2)  - X_U*1j*math.sin(phi/2) 	
	return V

def rotateY(phi):  # W
	Y_U = np.array([[0,-1j],[1j,0]])
	W = math.cos(phi/2)*I  - 1j*math.sin(phi/2)*Y_U 
	return W 


#checked ............................................................................


